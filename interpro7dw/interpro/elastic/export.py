import os
import pickle
import shelve
import shutil

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
                     protein2proteome_file: str, uniprot2pdb_file: str,
                     pdbmatches_file: str, alphafold_file: str,
                     proteomes_file: str, structures_file: str, clans_file: str,
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
        structures = pickle.load(fh)

    with open(uniprot2pdb_file, "rb") as fh:
        uniprot2pdb = pickle.load(fh)

    entry2pdb = {}
    with shelve.open(pdbmatches_file, writeback=False) as d:
        for pdb_chain, pdb_entry in d.items():
            for entry_acc, match in pdb_entry["matches"].items():
                try:
                    entry2pdb[entry_acc][pdb_chain] = match["locations"]
                except KeyError:
                    entry2pdb[entry_acc] = {pdb_chain: match["locations"]}

    logger.info("loading proteomes")
    with open(proteomes_file, "rb") as fh:
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
    proteomes_store = KVStore(protein2proteome_file)
    alphafold_store = KVStore(alphafold_file)
    domorgs_store = KVStore(domorgs_file)

    i = 0
    documents = []
    num_documents = 0
    seen_domains = set()
    seen_entries = set()
    seen_taxa = set()
    seen_structure2entry = set()
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
        pdb_documents = {}
        protein_structures = {}
        for pdb_chain, segments in uniprot2pdb.get(protein_acc, {}).items():
            pdb_id, chain_id = pdb_chain.split("_")

            try:
                protein_structures[pdb_id][chain_id] = segments
            except KeyError:
                protein_structures[pdb_id] = {chain_id: segments}

        for pdb_id, chains in protein_structures.items():
            structure = structures[pdb_id]
            pdb_base_doc = doc.copy()
            pdb_base_doc.update({
                "structure_acc": pdb_id.lower(),
                "structure_resolution": structure["resolution"],
                "structure_date": structure["date"],
                "structure_evidence": structure["evidence"],
                # "protein_structure": chains,
                "text_structure": join(pdb_id,
                                       structure["evidence"],
                                       structure["name"])
            })

            for chain_id, segments in chains.items():
                locations = []
                for segment in segments:
                    locations.append({
                        "fragments": [{
                            # Coordinates of UniProt entry on the PDB sequence
                            "start": segment["structure_start"],
                            "end": segment["structure_end"],
                            # Coordinates of PDB entry on the UniProt sequence
                            "protein_start": segment["protein_start"],
                            "protein_end": segment["protein_end"],
                        }]
                    })

                pdb_doc = pdb_base_doc.copy()
                pdb_doc.update({
                    "structure_chain_acc": chain_id,
                    "structure_protein_locations": locations,
                    "structure_chain": f"{pdb_id}-{chain_id}"
                })

                pdb_documents[f"{pdb_id}_{chain_id}"] = pdb_doc

        # Adds InterPro entries and member database signatures
        s_matches, e_matches = matches_store.get(protein_acc, ({}, {}))
        # structures/chains associated to at least one entry
        structures_with_entries = set()
        num_rel_docs = 0  # number of relationship documents

        for obj in [s_matches, e_matches]:
            for entry_acc, match in obj.items():
                locations = match["locations"]

                # TODO: remove after tests
                try:
                    entry = entries[entry_acc]
                except KeyError:
                    continue

                seen_entries.add(entry_acc)
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

                num_structures = 0
                entry_structures = entry2pdb.get(entry_acc, {})
                for pdb_chain, locations in entry_structures.items():
                    try:
                        pdb_doc = pdb_documents[pdb_chain]
                    except KeyError:
                        continue

                    # Associate entry to structure/chain
                    entry_doc = pdb_doc.copy()
                    entry_doc.update({
                        **entry_obj,
                        "entry_structure_locations": locations
                    })

                    documents.append((
                        entry_database + version,
                        get_rel_doc_id(entry_doc),
                        entry_doc
                    ))

                    seen_structure2entry.add((pdb_chain, entry_acc))
                    structures_with_entries.add(pdb_chain)
                    num_rel_docs += 1
                    num_structures += 1

                if not num_structures:
                    # Associates the entry to the protein, directly
                    entry_doc = doc.copy()
                    entry_doc.update(entry_obj)
                    documents.append((
                        entry_database + version,
                        get_rel_doc_id(entry_doc),
                        entry_doc
                    ))
                    num_rel_docs += 1

        for pdb_chain, pdb_doc in pdb_documents.items():
            if pdb_chain in structures_with_entries:
                continue

            # Structure without any entry
            pdb_doc.update({
                "ida_id": domain_id,
                "ida": domain_str
            })

            documents.append((
                # Not overlapping any entry -> not associated to a member DB
                config.REL_DEFAULT_INDEX + version,
                get_rel_doc_id(pdb_doc),
                pdb_doc
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

    """
    Add unseen entry-structure pairs 
    i.e. entries with PDB matches where the PDB structure is not associated
    to a protein in UniProtKB
    """
    for entry_acc, entry_structures in entry2pdb.items():
        for pdb_chain, locations in entry_structures.items():
            pdb_id, chain_id = pdb_chain.split("_")
            try:
                structure = structures[pdb_id]
            except KeyError:
                continue

            item = (pdb_chain, entry_acc)
            if item in seen_structure2entry:
                continue

            entry = entries[entry_acc]
            if not entry.public:
                continue

            database = entry.database.lower()
            if entry.integrated_in:
                integrated_in = entry.integrated_in.lower()
            else:
                integrated_in = None

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
                "structure_acc": pdb_id.lower(),
                "structure_resolution": structure["resolution"],
                "structure_date": structure["date"],
                "structure_evidence": structure["evidence"],
                "text_structure": join(pdb_id,
                                       structure["evidence"],
                                       structure["name"]),
                "structure_chain_acc": chain_id,
                "structure_chain": f"{pdb_id}-{chain_id}",
                "entry_structure_locations": locations
            })

            structure_taxa = structure["taxonomy"]
            if len(structure_taxa) == 1:
                taxon_id = structure_taxa[0]
                taxon = taxa[taxon_id]
                seen_taxa.add(taxon_id)

                doc.update({
                    "tax_id": taxon_id,
                    "tax_name": taxon["sci_name"],
                    "tax_lineage": taxon["lineage"],
                    "tax_rank": taxon["rank"],
                    "text_taxonomy": join(taxon_id, taxon["full_name"],
                                          taxon["rank"])
                })

            if entry_acc in member2clan:
                clan_acc, clan_name = member2clan[entry.accession]
                doc.update({
                    "set_acc": clan_acc.lower(),
                    "set_db": database,
                    "text_set": join(clan_acc, clan_name),
                })

            documents.append(doc)
            seen_entries.add(entry_acc)
            seen_structure2entry.add(item)

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
