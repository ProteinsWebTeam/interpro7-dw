import os
import pickle
import shutil

from interpro7dw.interpro.utils import overlaps_pdb_chain
from interpro7dw.utils import logger
from interpro7dw.utils.store import Directory, KVStore
from . import config


def init_rel_doc() -> dict:
    return {key: None for key in config.REL_BODY["mappings"]["properties"]}


def join(*args, separator: str = " ") -> str:
    items = []

    for item in args:
        if item is None:
            continue
        elif isinstance(item, (int, float)):
            item = str(item)
        elif isinstance(item, (list, set, tuple)):
            item = separator.join(map(str, item))
        elif isinstance(item, dict):
            item = separator.join(map(str, item.values()))
        elif not isinstance(item, str):
            continue

        items.append(item)

    return separator.join(items)


def get_rel_doc_id(doc: dict) -> str:
    if doc["protein_acc"]:
        return join(doc["protein_acc"],
                    doc["proteome_acc"],
                    doc["entry_acc"],
                    doc["set_acc"],
                    doc["structure_acc"],
                    doc["structure_chain_acc"],
                    separator='-')
    elif doc["entry_acc"]:
        return join(doc["entry_acc"], doc["set_acc"], separator='-')

    return doc["tax_id"]


def export_documents(proteins_file: str, matches_file: str, domorgs_file: str,
                     proteomes_file: str, structures_file: str,
                     alphafold_file: str, proteomeinfo_file: str,
                     structureinfo_file: str, clans_file: str,
                     entries_file: str, taxa_file: str, outdirs: list[str],
                     version: str, cachesize: int = 100000):
    directories = []
    for path in outdirs:
        try:
            shutil.rmtree(path)
        except FileNotFoundError:
            pass

        os.makedirs(path, mode=0o775)
        directories.append(Directory(root=path))
        open(os.path.join(path, f"{version}{config.LOAD_SUFFIX}"), "w").close()

    logger.info("loading PDBe data")
    with open(structures_file, "rb") as fh:
        protein2structures = pickle.load(fh)

    with open(structureinfo_file, "rb") as fh:
        structures = pickle.load(fh)["entries"]

    logger.info("loading proteomes")
    with open(proteomeinfo_file, "rb") as fh:
        proteomes = pickle.load(fh)

    logger.info("loading taxonomy")
    with open(taxa_file, "rb") as fh:
        taxa = pickle.load(fh)

    logger.info("loading entries and clans")
    with open(entries_file, "rb") as fh:
        entries = pickle.load(fh)

    member2clan = {}
    with open(clans_file, "rb") as fh:
        for clan in pickle.load(fh).values():
            for entry_acc, _, _, _, _ in clan["members"]:
                member2clan[entry_acc] = (clan["accession"], clan["name"])

    logger.info("writing documents")
    proteins_store = KVStore(proteins_file)
    matches_store = KVStore(matches_file)
    proteomes_store = KVStore(proteomes_file)
    alphafold_store = KVStore(alphafold_file)
    domorgs_store = KVStore(domorgs_file)

    i = 0
    documents = []
    num_documents = 0
    seen_domains = set()
    seen_entries = set()
    seen_taxa = set()
    for i, (protein_acc, protein) in enumerate(proteins_store.items()):
        taxon_id = protein["taxid"]
        taxon = taxa[taxon_id]
        seen_taxa.add(taxon_id)

        try:
            domain = domorgs_store[protein_acc]
        except KeyError:
            domain_id = domain_str = None
            domain_members = set()
        else:
            domain_id = domain["id"]
            domain_str = domain["key"]
            domain_members = domain["members"]
            if domain_id not in seen_domains:
                seen_domains.add(domain_id)
                locations = []
                for loc in domain["locations"]:
                    entry_acc = loc["pfam"]
                    locations.append({
                        "accession": entry_acc,
                        "name": entries[entry_acc].short_name,
                        "coordinates": [{
                            "fragments": [{
                                "start": loc["start"],
                                "end": loc["end"]
                            }]
                        }]
                    })

                    if entry_acc := loc["interpro"]:
                        locations.append({
                            "accession": entry_acc,
                            "name": entries[entry_acc].short_name,
                            "coordinates": [{
                                "fragments": [{
                                    "start": loc["start"],
                                    "end": loc["end"]
                                }]
                            }]
                        })

                documents.append((
                    config.IDA_INDEX + version,
                    domain["id"],
                    {
                        "ida_id": domain["id"],
                        "ida": domain["key"],
                        "representative": {
                            "accession": domain["protein"],
                            "length": domain["length"],
                            "domains": locations
                        },
                        "counts": domain["count"]
                    }
                ))

        af_models = alphafold_store.get(protein_acc, [])
        if af_models:
            # list of tuples (AFDB ID, score)
            best_model = sorted(af_models, key=lambda x: x[1])[-1]
            af_score = best_model[1]
        else:
            af_score = -1

        # Creates an empty document (all properties set to None)
        doc = init_rel_doc()
        doc.update({
            "protein_acc": protein_acc.lower(),
            "protein_length": protein["length"],
            "protein_is_fragment": protein["fragment"],
            "protein_af_score": af_score,
            "protein_db": "reviewed" if protein["reviewed"] else "unreviewed",
            "text_protein": join(protein_acc,
                                 protein["identifier"],
                                 taxon["sci_name"]),

            # Taxonomy
            "tax_id": taxon_id,
            "tax_name": taxon["sci_name"],
            "tax_lineage": taxon["lineage"],
            "tax_rank": taxon["rank"],
            "text_taxonomy": join(taxon_id, taxon["full_name"], taxon["rank"])
        })

        proteome_id = proteomes_store.get(protein_acc)
        if proteome_id:
            # Adds proteome
            proteome = proteomes[proteome_id]
            doc.update({
                "proteome_acc": proteome_id.lower(),
                "proteome_name": proteome["name"],
                "proteome_is_reference": proteome["is_reference"],
                "text_proteome": join(proteome_id,
                                      proteome["name"],
                                      proteome["assembly"],
                                      proteome["taxon_id"],
                                      proteome["strain"]),
            })

        # Adds PDBe structures and chains
        pdb_chains = {}  # mapping PDBe-chain ID -> chain segments
        pdb_documents = {}  # mapping PDBe-chain ID -> ES document
        protein_structures = protein2structures.get(protein_acc, {})
        for pdb_id, chains in protein_structures.items():
            structure = structures[pdb_id]

            pdb_doc = doc.copy()
            pdb_doc.update({
                "structure_acc": pdb_id.lower(),
                "structure_resolution": structure["resolution"],
                "structure_date": structure["date"],
                "structure_evidence": structure["evidence"],
                "protein_structure": chains,
                "text_structure": join(pdb_id,
                                       structure["evidence"],
                                       structure["name"])
            })

            for chain_id, segments in chains.items():
                pdb_chain_id = f"{pdb_id}-{chain_id}"

                locations = []
                for segment in segments:
                    locations.append({
                        "fragments": [{
                            "start": segment["protein_start"],
                            "end": segment["protein_end"],
                        }]
                    })

                chain_doc = pdb_doc.copy()
                chain_doc.update({
                    "structure_chain_acc": chain_id,
                    "structure_protein_locations": locations,
                    "structure_chain": pdb_chain_id
                })

                pdb_chains[pdb_chain_id] = segments
                pdb_documents[pdb_chain_id] = chain_doc

        # Adds InterPro entries and member database signatures
        s_matches, e_matches = matches_store.get(protein_acc, ({}, {}))
        overlapping_chains = set()  # chains associated to at least one entry
        num_rel_docs = 0  # number of relationship documents

        for obj in [s_matches, e_matches]:
            for entry_acc, match in obj.items():
                seen_entries.add(entry_acc)
                locations = match["locations"]

                entry = entries[entry_acc]
                entry_database = entry.database.lower()

                if entry_database == "panther":
                    """
                    PANTHER: remove the node ID 
                    (other databases do not have a subfamily property)
                    """
                    for loc in locations:
                        try:
                            del loc["subfamily"]["node"]
                        except KeyError:
                            continue  # No subfamily annotation

                if entry.integrated_in:
                    integrated_in = entry.integrated_in.lower()
                else:
                    integrated_in = None

                entry_obj = {
                    "entry_acc": entry_acc.lower(),
                    "entry_db": entry_database,
                    "entry_type": entry.type.lower(),
                    "entry_date": entry.creation_date.strftime("%Y-%m-%d"),
                    "entry_protein_locations": locations,
                    "entry_go_terms": [t["identifier"] for t in
                                       entry.go_terms],
                    "entry_integrated": integrated_in,
                    "text_entry": join(entry_acc, entry.short_name, entry.name,
                                       entry.type.lower(), integrated_in),
                }

                if entry_acc in member2clan:
                    clan_acc, clan_name = member2clan[entry_acc]
                    entry_obj.update({
                        "set_acc": clan_acc.lower(),
                        "set_db": entry_database,
                        "text_set": join(clan_acc, clan_name),
                    })

                if entry_acc in domain_members:
                    entry_obj.update({
                        "ida_id": domain_id,
                        "ida": domain_str,
                    })

                # Tests if the entry overlaps PDBe chains
                entry_chains = set()
                for pdb_chain_id, segments in pdb_chains.items():
                    if overlaps_pdb_chain(locations, segments):
                        # Entry overlaps chain: associate entry to struct/chain
                        chain_doc = pdb_documents[pdb_chain_id]
                        entry_doc = chain_doc.copy()
                        entry_doc.update(entry_obj)

                        documents.append((
                            entry_database + version,
                            get_rel_doc_id(entry_doc),
                            entry_doc
                        ))

                        entry_chains.add(pdb_chain_id)
                        num_rel_docs += 1

                if entry_chains:
                    # Stores chains that overlap at least one entry
                    overlapping_chains |= entry_chains
                else:
                    # Associates the entry to the protein, directly
                    entry_doc = doc.copy()
                    entry_doc.update(entry_obj)
                    documents.append((
                        entry_database + version,
                        get_rel_doc_id(entry_doc),
                        entry_doc
                    ))
                    num_rel_docs += 1

        # Adds chains that do not overlap with any entry
        for chain_id, chain_doc in pdb_documents.items():
            if chain_id in overlapping_chains:
                continue

            chain_doc.update({
                "ida_id": domain_id,
                "ida": domain_str
            })

            documents.append((
                # Not overlapping any entry -> not associated to a member DB
                config.REL_DEFAULT_INDEX + version,
                get_rel_doc_id(chain_doc),
                chain_doc
            ))
            num_rel_docs += 1

        if num_rel_docs == 0:
            # No relationships for this protein: fallback to protein doc
            documents.append((
                config.REL_DEFAULT_INDEX + version,
                get_rel_doc_id(doc),
                doc
            ))

        while len(documents) >= cachesize:
            for dir_obj in directories:
                filepath = dir_obj.mktemp()

                with open(filepath, "wb") as fh:
                    pickle.dump(documents[:cachesize], fh)

                os.rename(filepath, f"{filepath}{config.EXTENSION}")

            del documents[:cachesize]
            num_documents += cachesize

        if (i + 1) % 1e7 == 0:
            logger.info(f"{i + 1:>15,}")

    proteins_store.close()
    matches_store.close()
    proteomes_store.close()
    alphafold_store.close()
    domorgs_store.close()
    logger.info(f"{i + 1:>15,}")

    logger.info("writing remaining documents")

    # Adds unseen entries
    for entry in entries.values():
        if entry.accession in seen_entries or not entry.public:
            continue

        if entry.integrated_in:
            integrated_in = entry.integrated_in.lower()
        else:
            integrated_in = None

        database = entry.database.lower()
        doc = init_rel_doc()
        doc.update({
            "entry_acc": entry.accession.lower(),
            "entry_db": database,
            "entry_type": entry.type.lower(),
            "entry_date": entry.creation_date.strftime("%Y-%m-%d"),
            "entry_protein_locations": [],
            "entry_go_terms": [t["identifier"] for t in entry.go_terms],
            "entry_integrated": integrated_in,
            "text_entry": join(entry.accession, entry.short_name,
                               entry.name,
                               entry.type.lower(), integrated_in),
        })

        if entry.accession in member2clan:
            clan_acc, clan_name = member2clan[entry.accession]
            doc.update({
                "set_acc": clan_acc.lower(),
                "set_db": database,
                "text_set": join(clan_acc, clan_name),
            })

        documents.append((
            database + version,
            get_rel_doc_id(doc),
            doc
        ))

    # Adds unseen taxa
    for taxon in taxa.values():
        if taxon["id"] in seen_taxa:
            continue

        doc = init_rel_doc()
        doc.update({
            "tax_id": taxon["id"],
            "tax_name": taxon["full_name"],
            "tax_lineage": taxon["lineage"],
            "tax_rank": taxon["rank"],
            "text_taxonomy": join(taxon["id"], taxon["full_name"],
                                  taxon["rank"])
        })

        documents.append((
            config.REL_DEFAULT_INDEX + version,
            get_rel_doc_id(doc),
            doc
        ))

    num_documents += len(documents)
    while documents:
        for dir_obj in directories:
            filepath = dir_obj.mktemp()

            with open(filepath, "wb") as fh:
                pickle.dump(documents[:cachesize], fh)

            os.rename(filepath, f"{filepath}{config.EXTENSION}")

        del documents[:cachesize]

    # Adds done sentinel file
    for path in outdirs:
        open(os.path.join(path, f"{version}{config.DONE_SUFFIX}"), "w").close()

    logger.info(f"complete: {num_documents:,} documents")
