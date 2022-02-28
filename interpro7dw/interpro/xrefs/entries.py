import math
import multiprocessing as mp
import pickle
from typing import Optional

from interpro7dw import metacyc, uniprot
from interpro7dw.interpro import oracle
from interpro7dw.interpro.utils import copy_dict, overlaps_pdb_chain
from interpro7dw.utils import logger
from interpro7dw.utils.store import BasicStore, Directory, KVStore
from .utils import dump_to_tmp


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


def _process_entries(proteins_file: str, matches_file: str,
                     alphafold_file: str, proteomes_file: str,
                     domorgs_file: str, structures_file: str,
                     protein2enzymes: dict, protein2reactome: dict,
                     start: str, stop: Optional[str], workdir: Directory,
                     queue: mp.Queue):
    with open(structures_file, "rb") as fh:
        protein2structures = pickle.load(fh)

    proteins_store = KVStore(proteins_file)
    matches_store = KVStore(matches_file)
    alphafold_store = KVStore(alphafold_file)
    proteomes_store = KVStore(proteomes_file)
    domorgs_store = KVStore(domorgs_file)

    i = 0
    tmp_stores = {}
    xrefs = {}
    for protein_acc, (signatures, entries) in matches_store.range(start, stop):
        protein = proteins_store[protein_acc]
        protein_id = protein["identifier"]
        taxon_id = protein["taxid"]
        proteome_id = proteomes_store.get(protein_acc)
        structures = protein2structures.get(protein_acc, {})

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
            for entry_acc, entry in obj.items():
                if entry_acc in xrefs:
                    entry_xrefs = xrefs[entry_acc]
                else:
                    entry_xrefs = xrefs[entry_acc] = {
                        "dom_orgs": set(),
                        "enzymes": set(),
                        "matches": 0,
                        "reactome": set(),
                        "proteins": [],
                        "proteomes": set(),
                        "structures": set(),
                        "struct_models": {
                            "AlphaFold": 0
                        },
                        "taxa": {}
                    }

                entry_xrefs["matches"] += len(entry["locations"])
                entry_xrefs["proteins"].append((protein_acc, protein_id))

                if taxon_id in entry_xrefs["taxa"]:
                    entry_xrefs["taxa"][taxon_id] += 1
                else:
                    entry_xrefs["taxa"][taxon_id] = 1

                if entry_acc in domain_members:
                    entry_xrefs["dom_orgs"].add(domain_id)

                if proteome_id:
                    entry_xrefs["proteomes"].add(proteome_id)

                for pdbe_id, chains in structures.items():
                    for chain_id, segments in chains.items():
                        if overlaps_pdb_chain(entry["locations"], segments):
                            entry_xrefs["structures"].add(pdbe_id)
                            break  # Skip other chains

                if is_interpro:
                    if in_alphafold:
                        entry_xrefs["struct_models"]["AlphaFold"] += 1

                    for ecno in protein2enzymes.get(protein_acc, []):
                        entry_xrefs["enzymes"].add(ecno)

                    pathways = protein2reactome.get(protein_acc, [])
                    for pathway_id, pathway_name in pathways:
                        entry_xrefs["reactome"].add((pathway_id, pathway_name))

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

    queue.put((False, i))
    queue.put((True, tmp_stores))


def export_xrefs(uniprot_uri: str, proteins_file: str, matches_file: str,
                 alphafold_file: str, proteomes_file: str, domorgs_file: str,
                 struct_models_file: str, structures_file: str,
                 taxa_file: str, metacyc_file: str, output: str,
                 interpro_uri: Optional[str] = None, processes: int = 8,
                 tempdir: Optional[str] = None):
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
    :param struct_models_file: BasicStore file of structural models
    :param structures_file: File of protein-structures mapping
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
                             proteomes_file, domorgs_file, structures_file,
                             protein2enzymes, protein2reactome, start, stop,
                             workdir, queue))
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

    logger.info("loading structural models")
    struct_models = {}
    with BasicStore(struct_models_file, mode="r") as models:
        for model in models:
            signature_acc, entry_acc, algorithm = model[:3]

            try:
                obj = struct_models[algorithm]
            except KeyError:
                obj = struct_models[algorithm] = {}

            if signature_acc in obj:
                obj[signature_acc] += 1
            else:
                obj[signature_acc] = 1

            if entry_acc:
                if entry_acc in obj:
                    obj[entry_acc] += 1
                else:
                    obj[entry_acc] = 1

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
        i = 0
        n = len(entry2stores)
        for entry_acc in sorted(entry2stores):
            # Merge cross-references
            entry_xrefs = {}
            for entry_store in entry2stores[entry_acc]:
                for xrefs in entry_store:
                    copy_dict(xrefs, entry_xrefs, concat_or_incr=True)

            """
            Propagates number of proteins matched to ancestors,
            and build tree of taxonomic distribution.
            """
            entry_taxa = {}
            tree = {}
            while entry_xrefs["taxa"]:
                taxon_id, num_proteins = entry_xrefs["taxa"].popitem()

                # Propagates for all ancestors
                is_species = False
                node_id = taxon_id
                while node_id:
                    node = taxa[node_id]

                    if node["rank"] == "species":
                        """
                        The taxon (with mapped sequences) is a species 
                        or has a species in its lineage.
                        """
                        is_species = True

                    try:
                        entry_taxa[node_id] += num_proteins
                    except KeyError:
                        entry_taxa[node_id] = num_proteins

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
                "all": entry_taxa,
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

            # Add structural models
            models_xrefs = entry_xrefs["struct_models"]
            for algorithm, counts in struct_models.items():
                models_xrefs[algorithm] = counts.get(entry_acc, 0)

            store.write((entry_acc, entry_xrefs))

            i += 1
            if i % 1e4 == 0:
                logger.info(f"{i:>15,.0f} / {n}")

        logger.info(f"{i:>15,.0f} / {n}")

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


def _process_clans(member2clan: dict, proteins_file: str, matches_file: str,
                   proteomes_file: str, domorgs_file: str,
                   structures_file: str, start: str, stop: Optional[str],
                   workdir: Directory, queue: mp.Queue):
    with open(structures_file, "rb") as fh:
        protein2structures = pickle.load(fh)

    proteins_store = KVStore(proteins_file)
    matches_store = KVStore(matches_file)
    proteomes_store = KVStore(proteomes_file)
    domorgs_store = KVStore(domorgs_file)

    i = 0
    tmp_stores = {}
    xrefs = {}
    for protein_acc, (signatures, entries) in matches_store.range(start, stop):
        protein = proteins_store[protein_acc]
        taxon_id = protein["taxid"]
        proteome_id = proteomes_store.get(protein_acc)
        structures = protein2structures.get(protein_acc, {})

        try:
            domain = domorgs_store[protein_acc]
        except KeyError:
            domain_id = None
            domain_members = []
        else:
            domain_id = domain["id"]
            domain_members = domain["members"]

        for signature_acc, signature in signatures.items():
            if signature_acc not in member2clan:
                continue

            clan_acc, database = member2clan[signature_acc]
            if clan_acc in xrefs:
                clan_xrefs = xrefs[clan_acc]
            else:
                clan_xrefs = xrefs[clan_acc] = {
                    "dom_orgs": set(),
                    "entries": {
                        "all": set()
                    },
                    "proteins": [],
                    "proteomes": set(),
                    "structures": set(),
                    "taxa": set()
                }

            if signature_acc in domain_members:
                clan_xrefs["dom_orgs"].add(domain_id)

            clan_xrefs["entries"]["all"].add(signature_acc)
            if database in clan_xrefs["entries"]:
                clan_xrefs["entries"][database].add(signature_acc)
            else:
                clan_xrefs["entries"][database] = {signature_acc}

            clan_xrefs["proteins"].append(protein_acc)

            if proteome_id:
                clan_xrefs["proteomes"].add(proteome_id)

            for pdbe_id, chains in structures.items():
                for chain_id, segments in chains.items():
                    if overlaps_pdb_chain(signature["locations"], segments):
                        clan_xrefs["structures"].add(pdbe_id)
                        break  # Skip other chains

            clan_xrefs["taxa"].add(taxon_id)

        i += 1
        if i == 1e5:
            dump_to_tmp(xrefs, tmp_stores, workdir)
            queue.put((False, i))
            i = 0

    dump_to_tmp(xrefs, tmp_stores, workdir)
    proteins_store.close()
    matches_store.close()
    proteomes_store.close()
    domorgs_store.close()

    queue.put((False, i))
    queue.put((True, tmp_stores))


def export_clan_xrefs(clans_file: str, proteins_file: str, matches_file: str,
                      proteomes_file: str, domorgs_file: str,
                      structures_file: str, output: str, processes: int = 8,
                      tempdir: Optional[str] = None):
    logger.info("loading clan members")
    clans = {}
    member2clan = {}
    with open(clans_file, "rb") as fh:
        for clan_acc, clan in pickle.load(fh).items():
            database = clan["database"]
            members = set()

            for entry_acc, _, _ in clan["members"]:
                member2clan[entry_acc] = (clan_acc, database)
                members.add(entry_acc)

            clans[clan_acc] = (database, members)

    logger.info("iterating proteins")
    processes = max(1, processes - 1)
    with KVStore(proteins_file) as store:
        keys = store.get_keys()

    chunksize = math.ceil(len(keys) / processes)
    queue = mp.Queue()
    workers = []
    for i in range(processes):
        start = keys[i * chunksize]
        try:
            stop = keys[(i + 1) * chunksize]
        except IndexError:
            stop = None

        workdir = Directory(tempdir=tempdir)
        p = mp.Process(target=_process_clans,
                       args=(member2clan, proteins_file, matches_file,
                             proteomes_file, domorgs_file, structures_file,
                             start, stop, workdir, queue))
        p.start()
        workers.append((p, workdir))

    clan2stores = {}
    progress = 0
    milestone = step = 1e7
    work_done = 0
    while work_done < len(workers):
        is_done, obj = queue.get()
        if is_done:
            work_done += 1
            for clan_acc, clan_store in obj.items():
                if clan_acc in clan2stores:
                    clan2stores[clan_acc].append(clan_store)
                else:
                    clan2stores[clan_acc] = [clan_store]
        else:
            progress += obj
            if progress >= milestone:
                logger.info(f"{progress:>15,.0f}")
                milestone += step

    logger.info(f"{progress:>15,.0f}")
    for p, workdir in workers:
        p.join()

    logger.info("writing final file")
    with BasicStore(output, mode="w") as store:
        i = 0
        n = len(clan2stores)
        while clan2stores:
            clan_acc, clan_stores = clan2stores.popitem()
            clans.pop(clan_acc)

            # Merge cross-references
            clan_xrefs = {}
            for clan_store in clan_stores:
                for xrefs in clan_store:
                    copy_dict(xrefs, clan_xrefs, concat_or_incr=True)

            store.write((clan_acc, clan_xrefs))

            i += 1
            if i % 1e3 == 0:
                logger.info(f"{i:>15,.0f} / {n}")

        logger.info(f"{i:>15,.0f} / {n}")

        logger.info(f"{len(clans)} clans without cross-references")
        for clan_acc, (database, members) in clans.items():
            store.write((clan_acc, {
                "dom_orgs": set(),
                "entries": {
                    "all": members,
                    database: members
                },
                "proteins": [],
                "proteomes": set(),
                "structures": set(),
                "taxa": set()
            }))

    size = 0
    for p, workdir, in workers:
        size += workdir.get_size()
        workdir.remove()

    logger.info(f"temporary files: {size / 1024 ** 2:.0f} MB")
    logger.info("done")
