import re
from typing import Dict

import cx_Oracle

from interpro7dw import metacyc, uniprot
from interpro7dw.interpro.utils import overlaps_pdb_chain
from interpro7dw.utils import logger
from interpro7dw.utils.store import SimpleStore, Store
from interpro7dw.utils.store import copy_dict, loadobj
from interpro7dw.utils.tempdir import TemporaryDirectory


RANKS = [
    "superkingdom",
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species"
]


def _dump(refs: Dict, stores: Dict[str, SimpleStore], dst: TemporaryDirectory):
    while refs:
        entry_acc, entry_xrefs = refs.popitem()
        try:
            store = stores[entry_acc]
        except KeyError:
            file = dst.mktemp()
            store = stores[entry_acc] = SimpleStore(file)

        store.add(entry_xrefs, keep_open=False)


def _format_node(node: dict) -> dict:
    children = []

    while node["children"]:
        child_id, child = node["children"].popitem()
        children.append(_format_node(child))

    node["children"] = children

    return node


def dump_entries(ipr_url: str, unp_url: str, proteins_file: str,
                 matches_file: str, proteomes_file: str, domorgs_file: str,
                 structures_file: str, taxa_file: str, metacyc_file: str,
                 alphafold_file: str, xrefs_file: str, **kwargs):
    """Export InterPro entries and member database signatures cross-references.
    For each entry or signature, the following information is saved:
        - proteins matched (and number of matches)
        - proteomes
        - PDBe structures
        - taxa
        - domain organisations
        - ENZYME numbers
        - MetaCyc and Reactome pathways

    :param ipr_url: InterPro Oracle connection string.
    :param unp_url: UniProt Oracle connection string.
    :param proteins_file: Store file of protein info.
    :param matches_file: Store file of protein matches.
    :param proteomes_file: Store file of protein-proteome mapping.
    :param domorgs_file: Store file of domain organisations.
    :param structures_file: File of PDBe structures.
    :param taxa_file: File of taxonomic information.
    :param metacyc_file: MetaCyc tar archive.
    :param alphafold_file: CSV file listing AlphaFold predictions.
    :param xrefs_file: Output SimpleStore file.
    """
    tempdir = kwargs.get("tempdir")

    logger.info("loading proteins with AlphaFold predictions")
    alphafold = {}
    with open(alphafold_file, "rt") as fh:
        for line in fh:
            protein_acc, start, end, model_id = line.rstrip().split(',')
            try:
                alphafold[protein_acc] += 1
            except KeyError:
                alphafold[protein_acc] = 1

    alphafold = {protein_acc
                 for protein_acc, cnt in alphafold.items()
                 if cnt == 1}

    logger.info("loading Swiss-Prot data")
    protein2enzymes = uniprot.misc.get_swissprot2enzyme(unp_url)
    protein2reactome = uniprot.misc.get_swissprot2reactome(unp_url)

    # Creates mapping protein -> structure -> chain -> locations
    logger.info("loading PDBe structures")
    protein2structures = {}
    for pdbe_id, entry in loadobj(structures_file).items():
        for protein_acc, chains in entry["proteins"].items():
            try:
                protein2structures[protein_acc][pdbe_id] = chains
            except KeyError:
                protein2structures[protein_acc] = {pdbe_id: chains}

    logger.info("iterating proteins")
    proteins = Store(proteins_file, "r")
    matches = Store(matches_file, "r")
    proteomes = Store(proteomes_file, "r")
    domorgs = Store(domorgs_file, "r")

    i = 0
    entry2store = {}
    xrefs = {}
    tempdir = TemporaryDirectory(root=tempdir)
    for i, (protein_acc, protein_matches) in enumerate(matches.items()):
        protein = proteins[protein_acc]
        protein_id = protein["identifier"]
        taxon_id = protein["taxid"]
        proteome_id = proteomes.get(protein_acc)
        structures = protein2structures.get(protein_acc, {})
        try:
            _, dom_id, dom_members, _ = domorgs[protein_acc]
        except KeyError:
            dom_id = None
            dom_members = []

        in_alphafold = protein_acc in alphafold

        for entry_acc, locations in protein_matches.items():
            try:
                entry_xrefs = xrefs[entry_acc]
            except KeyError:
                entry_xrefs = xrefs[entry_acc] = {
                    "dom_orgs": set(),
                    "enzymes": set(),
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

            entry_xrefs["matches"] += len(locations)
            entry_xrefs["proteins"].append((protein_acc, protein_id))

            try:
                entry_xrefs["taxa"][taxon_id] += 1
            except KeyError:
                entry_xrefs["taxa"][taxon_id] = 1

            if entry_acc in dom_members:
                entry_xrefs["dom_orgs"].add(dom_id)

            if proteome_id:
                entry_xrefs["proteomes"].add(proteome_id)

            for pdbe_id, chains in structures.items():
                for chain_id, segments in chains.items():
                    if overlaps_pdb_chain(locations, segments):
                        entry_xrefs["structures"].add(pdbe_id)
                        break  # Skip other chains

            if re.fullmatch(r"IPR\d{6}", entry_acc):
                if in_alphafold:
                    entry_xrefs["struct_models"]["alphafold"] += 1

                for ecno in protein2enzymes.get(protein_acc, []):
                    entry_xrefs["enzymes"].add(ecno)

                pathways = protein2reactome.get(protein_acc, [])
                for pathway_id, pathway_name in pathways:
                    entry_xrefs["reactome"].add((pathway_id, pathway_name))

        if (i + 1) % 1e4 == 0:
            _dump(xrefs, entry2store, tempdir)

            if (i + 1) % 10e6 == 0:
                logger.info(f"{i + 1:>15,}")

    _dump(xrefs, entry2store, tempdir)
    logger.info(f"{i + 1:>15,}")

    proteins.close()
    matches.close()
    proteomes.close()
    domorgs.close()

    # Frees memory
    alphafold.clear()
    protein2enzymes.clear()
    protein2reactome.clear()
    protein2structures.clear()

    logger.info("loading MetaCyc pathways")
    ec2metacyc = metacyc.get_ec2pathways(metacyc_file)

    logger.info("loading taxa")
    taxa = loadobj(taxa_file)

    logger.info("loading structural models")
    con = cx_Oracle.connect(ipr_url)
    cur = con.cursor()
    cur.execute(
        """
        SELECT LOWER(ALGORITHM), METHOD_AC, COUNT(*)
        FROM INTERPRO.STRUCT_MODEL
        GROUP BY METHOD_AC, ALGORITHM
        """
    )
    struct_models = {}
    for algorithm, entry_acc, cnt in cur:
        try:
            struct_models[algorithm][entry_acc] = cnt
        except KeyError:
            struct_models[algorithm] = {entry_acc: cnt}

    cur.close()
    con.close()

    logger.info("writing final file")
    with SimpleStore(xrefs_file) as store:
        # Case-insensitive order
        for entry_acc in sorted(entry2store, key=lambda x: x.lower()):
            entry_xrefs = {}

            # Merge cross-references
            entry_store = entry2store[entry_acc]
            for xrefs in entry_store:
                copy_dict(xrefs, entry_xrefs, concat_or_incr=True)

            """
            Defines lineages for all taxa with direct matches, but with only
            major ranks.
            Some ranks will stay empty (None, None) because not all clades
            exist (e.g. no family between an order and a genus).
            """
            lineages = {}
            for taxon_id in entry_xrefs["taxa"]:
                lineage = lineages[taxon_id] = [(None, None)] * len(RANKS)
                taxon = taxa[taxon_id]

                for node_id in taxon["lineage"]:
                    node = taxa[node_id]
                    try:
                        i = RANKS.index(node["rank"])
                    except ValueError:
                        pass
                    else:
                        lineage[i] = (node_id, node["sci_name"])

            """
            Propagates number of proteins matched to ancestors,
            and build tree of taxonomic distribution.
            """
            entry_taxa = {}
            tree = {}
            while entry_xrefs["taxa"]:
                taxon_id, num_proteins = entry_xrefs["taxa"].popitem()

                # Propagates for all clades
                taxon = taxa[taxon_id]
                for node_id in taxon["lineage"]:
                    try:
                        entry_taxa[node_id] += num_proteins
                    except KeyError:
                        entry_taxa[node_id] = num_proteins

                # Add lineage of major ranks in tree
                lineage = lineages[taxon_id]
                obj = tree
                for i, rank in enumerate(RANKS):
                    node_id, node_name = lineage[i]

                    try:
                        node = obj[node_id]
                    except KeyError:
                        node = obj[node_id] = {
                            "id": node_id,
                            "rank": rank,
                            "name": node_name,
                            "proteins": 0,
                            "children": {}
                        }

                    node["proteins"] += num_proteins
                    obj = node["children"]  # descends into children

            entry_xrefs["taxa"] = {
                "all": entry_taxa,
                "tree": [_format_node(node) for node in tree.values()]
            }

            # Adds MetaCyc pathways
            pathways = set()
            for ecno in entry_xrefs["enzymes"]:
                for pathway_id, pathway_name in ec2metacyc.get(ecno, []):
                    pathways.add((pathway_id, pathway_name))

            entry_xrefs["metacyc"] = pathways

            # Add structural models
            models_xrefs = entry_xrefs["struct_models"]
            for algorithm, counts in struct_models.items():
                models_xrefs[algorithm] = counts.get(entry_acc, 0)

            store.add((entry_acc, entry_xrefs))

    logger.info(f"temporary files: {tempdir.size / 1024 / 1024:.0f} MB")
    tempdir.remove()

    logger.info("done")


def dump_proteomes(proteins_file: str, matches_file: str, proteomes_file: str,
                   domorgs_file: str, structures_file: str, entries_file: str,
                   xrefs_file: str, **kwargs):
    """Export proteome cross-references, that is:
        - proteins
        - taxa
        - PDBe structures
        - InterPro entries and member database signatures
        - domain organisations
        - clans

    :param proteins_file: Store file of protein info.
    :param matches_file: Store file of protein matches.
    :param proteomes_file: Store file of protein-proteome mapping.
    :param domorgs_file: Store file of domain organisations.
    :param structures_file: File of PDBe structures.
    :param entries_file: File of InterPro entries
        and member database signatures.
    :param xrefs_file: Output SimpleStore.
    """
    tempdir = kwargs.get("tempdir")

    logger.info("loading entries")
    entries = loadobj(entries_file)

    # Creates mapping protein -> structures
    logger.info("loading PDBe structures")
    protein2structures = {}
    for pdbe_id, entry in loadobj(structures_file).items():
        for protein_acc in entry["proteins"]:
            try:
                protein2structures[protein_acc].add(pdbe_id)
            except KeyError:
                protein2structures[protein_acc] = {pdbe_id}

    logger.info("iterating proteins")
    proteins = Store(proteins_file, "r")
    matches = Store(matches_file, "r")
    proteomes = Store(proteomes_file, "r")
    domorgs = Store(domorgs_file, "r")

    i = 0
    stores = {}
    xrefs = {}
    tempdir = TemporaryDirectory(root=tempdir)
    for i, (protein_acc, proteome_id) in enumerate(proteomes.items()):
        protein = proteins[protein_acc]
        taxon_id = protein["taxid"]

        try:
            proteome_xrefs = xrefs[proteome_id]
        except KeyError:
            proteome_xrefs = xrefs[proteome_id] = {
                "domain_architectures": set(),
                "entries": {},
                "proteins": 0,
                "sets": set(),
                "structures": set(),
                "taxa": set()
            }

        proteome_xrefs["proteins"] += 1

        try:
            _, dom_id, _, _ = domorgs[protein_acc]
        except KeyError:
            pass
        else:
            proteome_xrefs["domain_architectures"].add(dom_id)

        for entry_acc in matches.get(protein_acc, []):
            entry = entries[entry_acc]
            try:
                proteome_xrefs["entries"][entry.database].add(entry_acc)
            except KeyError:
                proteome_xrefs["entries"][entry.database] = {entry_acc}

            if entry.clan:
                proteome_xrefs["sets"].add(entry.clan["accession"])

        try:
            pdbe_ids = protein2structures[protein_acc]
        except KeyError:
            pass
        else:
            proteome_xrefs["structures"] |= pdbe_ids

        proteome_xrefs["taxa"].add(taxon_id)

        if (i + 1) % 1e4 == 0:
            _dump(xrefs, stores, tempdir)

            if (i + 1) % 100e6 == 0:
                logger.info(f"{i + 1:>15,}")

    _dump(xrefs, stores, tempdir)
    logger.info(f"{i + 1:>15,}")

    proteins.close()
    matches.close()
    proteomes.close()
    domorgs.close()

    logger.info("writing final file")
    with SimpleStore(xrefs_file) as store:
        for proteome_id in sorted(stores):
            proteome_xrefs = {}

            # Merge cross-references
            proteome_store = stores[proteome_id]
            for xrefs in proteome_store:
                copy_dict(xrefs, proteome_xrefs, concat_or_incr=True)

            store.add((proteome_id, proteome_xrefs))

    logger.info(f"temporary files: {tempdir.size / 1024 / 1024:.0f} MB")
    tempdir.remove()

    logger.info("done")


def _propagate(taxa: Dict[str, list], base_xrefs: Dict, xrefs: Dict) -> Dict:
    all_xrefs = {}
    for taxon_id in taxa:
        # Init empty xrefs
        taxon_xrefs = all_xrefs[taxon_id] = {}
        copy_dict(base_xrefs, taxon_xrefs)

    while xrefs:
        taxon_id, taxon_xrefs = xrefs.popitem()
        lineage = taxa[taxon_id]

        for taxon_id in lineage:
            copy_dict(taxon_xrefs, all_xrefs[taxon_id], concat_or_incr=True)

    return all_xrefs


def dump_taxa(proteins_file: str, matches_file: str, proteomes_file: str,
              structures_file: str, entries_file: str, taxa_file: str,
              xrefs_file: str, **kwargs):
    """Export cross-references for taxa:
        - proteins
        - proteomes
        - InterPro entries and member database signatures
        - PDBe structures

    :param proteins_file: Store file of protein info.
    :param matches_file: Store file of protein matches.
    :param proteomes_file: Store file of protein-proteome mapping.
    :param structures_file: File of PDBe structures.
    :param entries_file: File of InterPro entries
    :param taxa_file: File of taxonomic information.
    :param xrefs_file: Output SimpleStore
    """
    tempdir = kwargs.get("tempdir")

    logger.info("loading entries")
    entries = loadobj(entries_file)

    logger.info("loading taxa")
    taxa = {}
    for taxon_id, taxon in loadobj(taxa_file).items():
        taxa[taxon_id] = taxon["lineage"]

    # Creates mapping protein -> structure -> chain -> locations
    logger.info("loading PDBe structures")
    structures = {}
    for pdbe_id, entry in loadobj(structures_file).items():
        for protein_acc, chains in entry["proteins"].items():
            try:
                structures[protein_acc][pdbe_id] = chains
            except KeyError:
                structures[protein_acc] = {pdbe_id: chains}

    logger.info("iterating proteins")
    base_xrefs = {
        "databases": {},
        "proteins": {
            "all": 0,
            "databases": {},  # grouped by entry database
            "entries": {}     # grouped by entry
        },
        "proteomes": set(),
        "structures": {
            "all": set(),
            "entries": {}  # overlapping with these entries
        }
    }
    proteins = Store(proteins_file, "r")
    matches = Store(matches_file, "r")
    proteomes = Store(proteomes_file, "r")

    i = 0
    stores = {}
    xrefs = {}
    tempdir = TemporaryDirectory(root=tempdir)
    for i, (protein_acc, protein) in enumerate(proteins.items()):
        taxon_id = protein["taxid"]
        proteome_id = proteomes.get(protein_acc)
        protein_structures = structures.get(protein_acc, {})
        protein_matches = matches.get(protein_acc, {})

        try:
            taxon_xrefs = xrefs[taxon_id]
        except KeyError:
            taxon_xrefs = xrefs[taxon_id] = {}
            copy_dict(base_xrefs, taxon_xrefs)

        taxon_xrefs["proteins"]["all"] += 1
        taxon_xrefs["proteomes"].add(proteome_id)

        # Add structures, regardless of entry matches
        taxon_xrefs["structures"]["all"] |= set(protein_structures.keys())

        databases = set()
        for entry_acc, locations in protein_matches.items():
            entry = entries[entry_acc]
            entry_db = entry.database

            try:
                taxon_xrefs["databases"][entry_db].add(entry_acc)
            except KeyError:
                taxon_xrefs["databases"][entry_db] = {entry_acc}

            if entry_db not in databases:
                # Counts the protein once per database
                try:
                    taxon_xrefs["proteins"]["databases"][entry_db] += 1
                except KeyError:
                    taxon_xrefs["proteins"]["databases"][entry_db] = 1
                finally:
                    databases.add(entry_db)

            try:
                taxon_xrefs["proteins"]["entries"][entry_acc] += 1
            except KeyError:
                taxon_xrefs["proteins"]["entries"][entry_acc] = 1

            tax2structs = taxon_xrefs["structures"]["entries"]
            for pdbe_id, chains in protein_structures.items():
                for chain_id, segments in chains.items():
                    if overlaps_pdb_chain(locations, segments):
                        try:
                            tax2structs[entry_acc].add(pdbe_id)
                        except KeyError:
                            tax2structs[entry_acc] = {pdbe_id}
                        finally:
                            break  # Skip other chains

        if (i + 1) % 1e4 == 0:
            xrefs = _propagate(taxa, base_xrefs, xrefs)
            _dump(xrefs, stores, tempdir)

            if (i + 1) % 100e6 == 0:
                logger.info(f"{i + 1:>15,}")

    xrefs = _propagate(taxa, base_xrefs, xrefs)
    _dump(xrefs, stores, tempdir)
    logger.info(f"{i + 1:>15,}")

    proteins.close()
    matches.close()
    proteomes.close()

    logger.info("writing final file")
    with SimpleStore(xrefs_file) as store:
        for taxon_id in sorted(stores):
            taxon_xrefs = {}

            # Merge cross-references
            taxon_store = stores[taxon_id]
            for xrefs in taxon_store:
                copy_dict(xrefs, taxon_xrefs, concat_or_incr=True)

            store.add((taxon_id, taxon_xrefs))

    logger.info(f"temporary files: {tempdir.size / 1024 / 1024:.0f} MB")
    tempdir.remove()

    logger.info("done")
