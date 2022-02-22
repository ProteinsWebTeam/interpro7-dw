import re
from typing import Dict

from interpro7dw import alphafold, metacyc, uniprot
from interpro7dw.interpro.utils import overlaps_pdb_chain
from interpro7dw.utils import logger
from interpro7dw.utils.store import Directory, SimpleStore, Store
from interpro7dw.utils.store import copy_dict, loadobj


def _dump(refs: Dict, stores: Dict[str, SimpleStore], dst: Directory):
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


def dump_entries(url: str, proteins_file: str, matches_file: str,
                 proteomes_file: str, domorgs_file: str, structures_file: str,
                 taxa_file: str, metacyc_file: str, alphafold_file: str,
                 struct_models_file: str, xrefs_file: str, **kwargs):
    """Export InterPro entries and member database signatures cross-references.
    For each entry or signature, the following information is saved:
        - proteins matched (and number of matches)
        - proteomes
        - PDBe structures
        - taxa
        - domain organisations
        - ENZYME numbers
        - MetaCyc and Reactome pathways

    :param url: UniProt Oracle connection string.
    :param proteins_file: Store file of protein info.
    :param matches_file: Store file of protein matches.
    :param proteomes_file: Store file of protein-proteome mapping.
    :param domorgs_file: Store file of domain organisations.
    :param structures_file: File of PDBe structures.
    :param taxa_file: File of taxonomic information.
    :param metacyc_file: MetaCyc tar archive.
    :param alphafold_file: CSV file listing AlphaFold predictions.
    :param struct_models_file: SimpleStore file of structural models.
    :param xrefs_file: Output SimpleStore file.
    """
    tempdir = kwargs.get("tempdir")

    logger.info("loading proteins with AlphaFold predictions")
    af_proteins = alphafold.get_proteins(alphafold_file, keep_fragments=False)

    logger.info("loading Swiss-Prot data")
    protein2enzymes = uniprot.misc.get_swissprot2enzyme(url)
    protein2reactome = uniprot.misc.get_swissprot2reactome(url)

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
    tempdir = Directory(tempdir=tempdir)
    for i, (protein_acc, protein_matches) in enumerate(matches.items()):
        protein = proteins[protein_acc]
        protein_id = protein["identifier"]
        taxon_id = protein["taxid"]
        proteome_id = proteomes.get(protein_acc)
        structures = protein2structures.get(protein_acc, {})
        try:
            domain = domorgs[protein_acc]
        except KeyError:
            dom_id = None
            dom_members = []
        else:
            dom_id = domain["id"]
            dom_members = domain["members"]

        in_alphafold = protein_acc in af_proteins

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
    af_proteins.clear()
    protein2enzymes.clear()
    protein2reactome.clear()
    protein2structures.clear()

    logger.info("loading MetaCyc pathways")
    ec2metacyc = metacyc.get_ec2pathways(metacyc_file)

    logger.info("loading taxa")
    taxa = loadobj(taxa_file)

    logger.info("loading structural models")
    struct_models = {}
    with SimpleStore(struct_models_file) as models:
        for model in models:
            entry_acc = model[0]
            algorithm = model[1]

            try:
                struct_models[algorithm][entry_acc] += 1
            except KeyError:
                struct_models[algorithm] = {entry_acc: 1}

    logger.info("writing final file")
    main_ranks = [
        "superkingdom",
        "kingdom",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species"
    ]
    with SimpleStore(xrefs_file) as store:
        for entry_acc, entry_store in entry2store.items():
            # Merge cross-references
            entry_xrefs = {}
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
                lineage = lineages[taxon_id] = [(None, None)] * len(main_ranks)
                taxon = taxa[taxon_id]

                for node_id in taxon["lineage"]:
                    node = taxa[node_id]
                    try:
                        i = main_ranks.index(node["rank"])
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
                lineage = lineages[taxon_id]
                obj = tree
                unique_id = "1"  # default to root
                for i, rank in enumerate(main_ranks):
                    node_id, node_name = lineage[i]

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
                        node = obj[unique_id] = {
                            "id": unique_id,
                            "rank": rank,
                            "name": node_name,
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




def dump_structures(proteins_file: str, matches_file: str, proteomes_file: str,
                    domorgs_file: str, structures_file: str, entries_file: str,
                    xrefs_file: str):

    logger.info("loading entries")
    entries = {}
    for entry in loadobj(entries_file).values():
        entries[entry.accession] = (entry.database, entry.clan)

    logger.info("loading PDBe structures")
    protein2structures = {}
    xrefs = {}
    for pdbe_id, entry in loadobj(structures_file).items():
        for protein_acc, chains in entry["proteins"].items():
            try:
                protein2structures[protein_acc][pdbe_id] = chains
            except KeyError:
                protein2structures[protein_acc] = {pdbe_id: chains}

        xrefs[pdbe_id] = {
            "domain_architectures": set(),
            "entries": {},
            "proteomes": set(),
            "proteins": 0,
            "sets": set(),
            "taxa": set()
        }

    logger.info("iterating proteins")
    proteins = Store(proteins_file, "r")
    matches = Store(matches_file, "r")
    proteomes = Store(proteomes_file, "r")
    domorgs = Store(domorgs_file, "r")

    i = 0
    for i, protein_acc in enumerate(sorted(protein2structures)):
        protein_structures = protein2structures[protein_acc]

        protein = proteins[protein_acc]
        taxon_id = protein["taxid"]
        proteome_id = proteomes.get(protein_acc)
        protein_matches = matches.get(protein_acc, {})

        try:
            dom_id = domorgs[protein_acc]["id"]
        except KeyError:
            dom_id = None

        for pdbe_id, chains in protein_structures.items():
            struct_xrefs = xrefs[pdbe_id]

            if dom_id:
                struct_xrefs["domain_architectures"].add(dom_id)

            if proteome_id:
                struct_xrefs["proteomes"].add(proteome_id)

            struct_xrefs["proteins"] += 1
            struct_xrefs["taxa"].add(taxon_id)

            for entry_acc, locations in protein_matches.items():
                database, clan = entries[entry_acc]

                for chain_id, segments in chains.items():
                    if overlaps_pdb_chain(locations, segments):
                        try:
                            struct_xrefs["entries"][database].add(entry_acc)
                        except KeyError:
                            struct_xrefs["entries"][database] = {entry_acc}

                        if clan:
                            struct_xrefs["sets"].add(clan["accession"])

                        break

        if (i + 1) % 1e4 == 0:
            logger.info(f"{i+1:>15,}")

    logger.info(f"{i + 1:>15,}")

    proteins.close()
    matches.close()
    proteomes.close()
    domorgs.close()

    logger.info("writing final file")
    with SimpleStore(xrefs_file) as store:
        for pdbe_id, struct_xrefs in xrefs.items():
            store.add((pdbe_id, struct_xrefs))

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
    tempdir = Directory(tempdir=tempdir)
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

        if (i + 1) % 1e6 == 0:
            _dump(xrefs, stores, tempdir)

            if (i + 1) % 10e6 == 0:
                logger.info(f"{i + 1:>15,}")

    _dump(xrefs, stores, tempdir)
    logger.info(f"{i + 1:>15,}")

    proteins.close()
    matches.close()
    proteomes.close()

    logger.info("loading taxa")
    taxa = loadobj(taxa_file)
    all_stores = {}
    for taxon_id, taxon in taxa.items():
        taxa[taxon_id] = set(taxon["lineage"])

        # Init final store
        all_stores[taxon_id] = SimpleStore(tempdir.mktemp())

    logger.info("propagate to ancestors")
    while stores:
        taxon_id, taxon_store = stores.popitem()

        # Merge cross-references
        taxon_xrefs = {}
        for xrefs in taxon_store:
            copy_dict(xrefs, taxon_xrefs, concat_or_incr=True)

        # Propagate to lineage (including current taxon!)
        lineage = taxa[taxon_id]
        for taxon_id in lineage:
            store = all_stores[taxon_id]
            store.add(taxon_xrefs, keep_open=False)

    logger.info("writing final file")
    with SimpleStore(xrefs_file) as final_store:
        while all_stores:
            taxon_id, taxon_store = all_stores.popitem()

            # Default cross-references
            taxon_xrefs = {}
            copy_dict(base_xrefs, taxon_xrefs)

            # Merge cross-references from descendants
            for xrefs in taxon_store:
                copy_dict(xrefs, taxon_xrefs, concat_or_incr=True)

            final_store.add((taxon_id, taxon_xrefs))

    logger.info(f"temporary files: {tempdir.size / 1024 / 1024:.0f} MB")
    tempdir.remove()

    logger.info("done")



