import math
import multiprocessing as mp
import pickle

from interpro7dw import metacyc, uniprot
from interpro7dw.interpro import oracle
from interpro7dw.interpro.utils import copy_dict
from interpro7dw.utils import logger
from interpro7dw.utils.store import BasicStore, Directory, KVStore
from .utils import dump_to_tmp, unpack_pdb_matches


MIN_SIMILARITY = 0.75
MAIN_RANKS = [
    "superkingdom",
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species"
]


def export_sim_entries(matches_file: str, output: str):
    logger.info("starting")

    num_proteins = {}  # number of proteins matched by entry
    num_overlaps = {}  # number of proteins where two entries overlap >= 50%
    entry2type = {}    # type (in lower case) of each InterPro entry

    with KVStore(matches_file) as matches:
        for i, (signatures, entries) in enumerate(matches.values()):
            entry2locations = {}

            for entry_acc, entry in entries.items():
                if entry_acc not in entry2type:
                    entry2type[entry_acc] = entry["type"].lower()

                # Flatten the entry's matches (only one frag / location)
                entry2locations[entry_acc] = []
                for loc in entry["locations"]:
                    entry2locations[entry_acc].append((
                        loc["fragments"][0]["start"],
                        loc["fragments"][0]["end"],
                    ))

            # Evaluate how InterPro entries overlap
            for entry_acc, locations in entry2locations.items():
                try:
                    num_proteins[entry_acc] += 1
                except KeyError:
                    num_proteins[entry_acc] = 1

                for other_acc, other_locations in entry2locations.items():
                    if other_acc >= entry_acc:
                        continue

                    try:
                        entry_overlaps = num_overlaps[entry_acc]
                    except KeyError:
                        entry_overlaps = num_overlaps[entry_acc] = {}

                    try:
                        overlaps = entry_overlaps[other_acc]
                    except KeyError:
                        overlaps = entry_overlaps[other_acc] = [0, 0]

                    flag = 0
                    for start1, end1 in locations:
                        length1 = end1 - start1 + 1

                        for start2, end2 in other_locations:
                            length2 = end2 - start2 + 1
                            overlap = min(end1, end2) - max(start1, start2) + 1

                            if not flag & 1 and overlap >= length1 * 0.5:
                                # 1st time fragments overlap >= 50% of entry 1
                                flag |= 1
                                overlaps[0] += 1

                            if not flag & 2 and overlap >= length2 * 0.5:
                                # 1st time fragments overlap >= 50% of entry 2
                                flag |= 2
                                overlaps[1] += 1

                        if flag == 3:
                            # Both cases already happened
                            break

            if (i + 1) % 1e7 == 0:
                logger.info(f"{i + 1:>15,}")

        logger.info(f"{i + 1:>15,}")

    supfam = "homologous_superfamily"
    types = (supfam, "domain", "family", "repeat")
    overlapping_entries = []
    for entry_acc, overlaps in num_overlaps.items():
        entry_cnt = num_proteins[entry_acc]

        for other_acc, (cnt1, cnt2) in overlaps.items():
            other_cnt = num_proteins[other_acc]

            # Independent coefficients
            coef1 = cnt1 / (entry_cnt + other_cnt - cnt1)
            coef2 = cnt2 / (entry_cnt + other_cnt - cnt2)

            # Final coefficient (average of independent coefficients)
            coef = (coef1 + coef2) * 0.5

            # Containment indices
            cont1 = cnt1 / entry_cnt
            cont2 = cnt2 / other_cnt

            if all(e < MIN_SIMILARITY for e in (coef, cont1, cont2)):
                continue

            # Entries are deemed similar
            type1 = entry2type[entry_acc]
            type2 = entry2type[other_acc]
            if ((type1 == supfam and type2 in types)
                    or (type2 == supfam and type1 in types)):
                overlapping_entries.append((entry_acc, other_acc))

    with open(output, "wb") as fh:
        pickle.dump(overlapping_entries, fh)

    logger.info("done")


def _init_entry_xrefs() -> dict:
    return {
        "dom_orgs": set(),
        "enzymes": set(),
        "genes": set(),
        "matches": 0,
        "reactome": set(),
        "proteins": [],
        "proteomes": set(),
        "structures": set(),
        "struct_models": {
            "alphafold": 0
        },
        "taxa": {}
    }


def _process_entries(proteins_file: str, matches_file: str,
                     alphafold_file: str, proteomes_file: str,
                     domorgs_file: str, pdb2matches_file: str,
                     evidences_file: str, protein2enzymes: dict,
                     protein2reactome: dict, start: str, stop: str | None,
                     workdir: Directory, queue: mp.Queue):
    entry2structures = unpack_pdb_matches(pdb2matches_file)
    proteins_store = KVStore(proteins_file)
    matches_store = KVStore(matches_file)
    alphafold_store = KVStore(alphafold_file)
    proteomes_store = KVStore(proteomes_file)
    domorgs_store = KVStore(domorgs_file)
    evidences_store = KVStore(evidences_file)

    i = 0
    tmp_stores = {}
    xrefs = {}
    for protein_acc, (signatures, entries) in matches_store.range(start, stop):
        protein = proteins_store[protein_acc]
        protein_id = protein["identifier"]
        taxon_id = protein["taxid"]
        proteome_id = proteomes_store.get(protein_acc)
        evidence, gene_name = evidences_store[protein_acc]

        try:
            domain = domorgs_store[protein_acc]
        except KeyError:
            domain_id = None
            domain_members = []
        else:
            domain_id = domain["id"]
            domain_members = domain["members"]

        in_alphafold = len(alphafold_store.get(protein_acc, [])) > 0

        for is_interpro, obj in [(False, signatures), (True, entries)]:
            for entry_acc, match in obj.items():
                if entry_acc in xrefs:
                    entry = xrefs[entry_acc]
                else:
                    entry = xrefs[entry_acc] = _init_entry_xrefs()

                entry["matches"] += len(match["locations"])
                entry["proteins"].append((protein_acc, protein_id,
                                          in_alphafold))

                if taxon_id in entry["taxa"]:
                    entry["taxa"][taxon_id] += 1
                else:
                    entry["taxa"][taxon_id] = 1

                if entry_acc in domain_members:
                    entry["dom_orgs"].add(domain_id)

                if proteome_id:
                    entry["proteomes"].add(proteome_id)

                """
                Use `pop()` instead of `get()` so we add structures once
                per signature 
                """
                for pdb_id, v in entry2structures.pop(entry_acc, {}).items():
                    ratio = v["coverage"] / v["length"]
                    entry["structures"].add((pdb_id, ratio))

                if gene_name:
                    entry["genes"].add(gene_name)

                if in_alphafold:
                    entry["struct_models"]["alphafold"] += 1

                if is_interpro:
                    for ecno in protein2enzymes.get(protein_acc, []):
                        entry["enzymes"].add(ecno)

                    pathways = protein2reactome.get(protein_acc, [])
                    for pathway_id, pathway_name in pathways:
                        entry["reactome"].add((pathway_id, pathway_name))

        i += 1
        if i == 1e5:
            dump_to_tmp(xrefs, tmp_stores, workdir)
            queue.put((False, i))
            i = 0

    dump_to_tmp(xrefs, tmp_stores, workdir)
    proteins_store.close()
    matches_store.close()
    alphafold_store.close()
    proteomes_store.close()
    domorgs_store.close()
    evidences_store.close()

    queue.put((False, i))
    queue.put((True, tmp_stores))


def export_xrefs(uniprot_uri: str, proteins_file: str, matches_file: str,
                 alphafold_file: str, proteomes_file: str, domorgs_file: str,
                 pdb2matches_file: str, evidences_file: str, taxa_file: str,
                 metacyc_file: str, output: str,
                 interpro_uri: str | None = None, processes: int = 8,
                 tempdir: str | None = None):
    """Export InterPro entries and member database signatures cross-references.
    For each entry or signature, the following information is saved:
        - proteins matched (and number of matches)
        - proteomes
        - PDBe structures
        - taxa
        - domain organisations
        - ENZYME numbers
        - MetaCyc and Reactome pathways

    :param uniprot_uri: UniProt Oracle connection string
    :param proteins_file: KVStore file of protein info
    :param matches_file: KVStore file of protein matches
    :param alphafold_file: KVStore file of proteins with AlphaFold models
    :param proteomes_file: KVStore file of protein-proteome mapping
    :param domorgs_file: KVStore file of domain organisations
    :param pdb2matches_file: File of PDB matches
    :param evidences_file: KVStore file of protein evidences/genes
    :param taxa_file: File of taxonomic information
    :param metacyc_file: MetaCyc tar archive
    :param output: Output BasicStore file
    :param interpro_uri: InterPro Oracle connection string. If provided,
    update database with pathways cross-references
    :param processes: Number of workers
    :param tempdir: Temporary directory
    """
    logger.info("loading Swiss-Prot data")
    protein2enzymes = uniprot.proteins.get_swissprot2enzyme(uniprot_uri)
    protein2reactome = uniprot.proteins.get_swissprot2reactome(uniprot_uri)

    logger.info("iterating proteins")
    processes = max(1, processes - 1)
    with KVStore(proteins_file) as store:
        keys = store.get_keys()

    chunksize = math.ceil(len(keys) / processes)
    queue = mp.Queue()
    workers = []
    for i in range(processes):
        start = keys[i*chunksize]
        try:
            stop = keys[(i+1)*chunksize]
        except IndexError:
            stop = None

        workdir = Directory(tempdir=tempdir)
        p = mp.Process(target=_process_entries,
                       args=(proteins_file, matches_file, alphafold_file,
                             proteomes_file, domorgs_file, pdb2matches_file,
                             evidences_file, protein2enzymes,
                             protein2reactome, start, stop, workdir, queue))
        p.start()
        workers.append((p, workdir))

    entry2stores = {}
    progress = 0
    milestone = step = 1e7
    work_done = 0
    while work_done < len(workers):
        is_done, obj = queue.get()
        if is_done:
            work_done += 1
            for entry_acc, entry_store in obj.items():
                if entry_acc in entry2stores:
                    entry2stores[entry_acc].append(entry_store)
                else:
                    entry2stores[entry_acc] = [entry_store]
        else:
            progress += obj
            if progress >= milestone:
                logger.info(f"{progress:>15,.0f}")
                milestone += step

    logger.info(f"{progress:>15,.0f}")
    for p, workdir in workers:
        p.join()

    logger.info("loading MetaCyc pathways")
    ec2metacyc = metacyc.get_ec2pathways(metacyc_file)

    logger.info("loading taxonomy")
    with open(taxa_file, "rb") as fh:
        taxa = pickle.load(fh)

    """
    Define lineage but with major ranks only.
    Some ranks will stay empty (None) because not all clades 
    exist (e.g. no family between an order and a genus).
    """
    logger.info("create lineage for main taxonomic ranks")
    for info in taxa.values():
        lineage = [None] * len(MAIN_RANKS)

        for node_id in info["lineage"]:
            node = taxa[node_id]

            try:
                i = MAIN_RANKS.index(node["rank"])
            except ValueError:
                pass
            else:
                lineage[i] = node_id

        info["main_ranks"] = lineage

    logger.info("writing final file")
    entry2pathways = {}
    with BasicStore(output, mode="w") as store:
        progress = 0
        total = len(entry2stores)
        milestone = step = math.ceil(0.1 * total)
        for entry_acc in sorted(entry2stores):
            # Merge cross-references
            entry_xrefs = {}
            for entry_store in entry2stores[entry_acc]:
                for xrefs in entry_store:
                    copy_dict(xrefs, entry_xrefs, concat_or_incr=True)

            """
            Propagate protein count to ancestors and build a tree of
            taxonomy distribution
            """
            all_taxa = {}
            tree = {}
            for taxon_id, num_proteins in entry_xrefs["taxa"].items():
                is_species = False
                node_id = taxon_id
                while node_id:
                    node = taxa[node_id]

                    try:
                        all_taxa[node_id] += num_proteins
                    except KeyError:
                        all_taxa[node_id] = num_proteins

                    if node["rank"] == "species":
                        """
                        The taxon (with mapped sequences) 
                        or one of its ancestors is a species .
                        """
                        is_species = True

                    node_id = node["parent"]

                # Add lineage of major ranks in tree
                lineage = taxa[taxon_id]["main_ranks"]
                obj = tree
                unique_id = "1"  # default to root
                for i, (rank, node_id) in enumerate(zip(MAIN_RANKS, lineage)):
                    """
                    Since several nodes may have node_id set to None,
                    we need to create a unique identifier
                    """
                    if node_id:
                        unique_id = node_id
                    else:
                        unique_id += f"-{i}"

                    try:
                        node = obj[unique_id]
                    except KeyError:
                        if node_id:
                            sci_name = taxa[node_id]["sci_name"]
                        else:
                            sci_name = None

                        node = obj[unique_id] = {
                            "id": unique_id,
                            "rank": rank,
                            "name": sci_name,
                            "proteins": 0,
                            "species": 0,
                            "children": {}
                        }

                    node["proteins"] += num_proteins
                    if is_species:
                        node["species"] += 1

                    obj = node["children"]  # descends into children

            # Wraps superkingdoms in a "root" node
            num_proteins = 0
            num_species = 0
            children = []
            for node in tree.values():
                num_proteins += node["proteins"]
                num_species += node["species"]
                children.append(_format_node(node))

            entry_xrefs["taxa"] = {
                "all": all_taxa,
                "hit": entry_xrefs["taxa"],
                "tree": {
                    "id": "1",
                    "rank": None,
                    "name": "root",
                    "proteins": num_proteins,
                    "species": num_species,
                    "children": children
                }
            }

            # Adds MetaCyc pathways
            entry_xrefs["metacyc"] = set()
            for ecno in entry_xrefs["enzymes"]:
                for pathway_id, pathway_name in ec2metacyc.get(ecno, []):
                    entry_xrefs["metacyc"].add((pathway_id, pathway_name))

            entry2pathways[entry_acc] = [
                ("metacyc", entry_xrefs["metacyc"]),
                ("reactome", entry_xrefs["reactome"])
            ]

            store.write((entry_acc, entry_xrefs))

            progress += 1
            if progress == milestone:
                logger.info(f"{progress:>15,.0f} / {total:,}")
                milestone += step

        logger.info(f"{progress:>15,.0f} / {total:,}")

    size = 0
    for p, workdir, in workers:
        size += workdir.get_size()
        workdir.remove()

    logger.info(f"temporary files: {size / 1024 ** 2:.0f} MB")

    if interpro_uri:
        logger.info("updating ENTRY2PATHWAY")
        oracle.entries.update_pathways(interpro_uri, entry2pathways)

    logger.info("done")


def _format_node(node: dict) -> dict:
    children = []

    while node["children"]:
        child_id, child = node["children"].popitem()
        children.append(_format_node(child))

    node["children"] = children

    return node
