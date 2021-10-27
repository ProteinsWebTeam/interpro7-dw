import os
import shutil
from typing import Sequence

from interpro7dw import alphafold
from interpro7dw.interpro.utils import overlaps_pdb_chain
from interpro7dw.utils import logger
from interpro7dw.utils.store import Directory, Store, dumpobj, loadobj
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
                     proteomes_file: str, entries_file: str,
                     proteomeinfo_file: str, structures_file: str,
                     taxa_file: str, alphafold_file: str,
                     outdirs: Sequence[str], version: str,
                     cache_size: int = 100000):
    directories = []
    for path in outdirs:
        try:
            shutil.rmtree(path)
        except FileNotFoundError:
            pass

        os.makedirs(path)
        directories.append(Directory(root=path))
        open(os.path.join(path, f"{version}{config.LOAD_SUFFIX}"), "w").close()

    logger.info("loading domain organisations")
    domorgs = {}
    with Store(domorgs_file, "r") as store:
        for domain in store.values():
            if domain["id"] not in domorgs:
                domorgs["id"] = domain

    logger.info("writing domain organisation documents")
    num_documents = 0
    domorgs = list(domorgs.values())
    for i in range(0, len(domorgs), cache_size):
        documents = []
        for dom in domorgs[i:i + cache_size]:
            locations = []
            for loc in dom["locations"]:
                locations.append({
                    "entry": loc["pfam"],
                    "coordinates": [{
                        "fragments": [{
                            "start": loc["start"],
                            "end": loc["end"]
                        }]
                    }]
                })

                if loc["interpro"]:
                    locations.append({
                        "entry": loc["interpro"],
                        "coordinates": [{
                            "fragments": [{
                                "start": loc["start"],
                                "end": loc["end"]
                            }]
                        }]
                    })

            documents.append((
                config.IDA_INDEX_PREFIX + version,
                dom["id"],
                {
                    "ida_id": dom["id"],
                    "ida": dom["key"],
                    "representative": {
                        "accession": dom["protein"],
                        "length": dom["length"],
                        "domains": locations
                    },
                    "counts": dom["count"]
                }
            ))

        num_documents += len(documents)
        for dir_obj in directories:
            filepath = dir_obj.mktemp()
            dumpobj(documents, filepath)
            os.rename(filepath, f"{filepath}{config.EXTENSION}")

    logger.info("loading proteins with AlphaFold predictions")
    af_proteins = alphafold.get_proteins(alphafold_file, keep_fragments=False)

    logger.info("loading other data")
    entries = loadobj(entries_file)
    proteomesinfo = loadobj(proteomeinfo_file)
    structures = loadobj(structures_file)
    taxa = loadobj(taxa_file)

    protein2structures = {}
    for pdbe_id, entry in loadobj(structures_file).items():
        for protein_acc in entry["proteins"]:
            try:
                protein2structures[protein_acc].add(pdbe_id)
            except KeyError:
                protein2structures[protein_acc] = {pdbe_id}

    logger.info("writing relationship documents")

    proteins = Store(proteins_file, "r")
    matches = Store(matches_file, "r")
    proteomes = Store(proteomes_file, "r")
    domorgs = Store(domorgs_file, "r")
    documents = []
    used_entries = set()
    used_taxa = set()
    i = 0
    for i, (protein_acc, protein) in enumerate(proteins.items()):
        taxon_id = protein["taxid"]
        taxon = taxa[taxon_id]
        used_taxa.add(taxon_id)

        try:
            domain = domorgs[protein_acc]
        except KeyError:
            dom_id = dom_key = None
            dom_members = []
        else:
            dom_id = domain["id"]
            dom_key = domain["key"]
            dom_members = domain["members"]

        # Creates an empty document (all properties set to None)
        doc = init_rel_doc()
        doc.update({
            "protein_acc": protein_acc.lower(),
            "protein_length": protein["length"],
            "protein_is_fragment": protein["fragment"],
            "protein_has_model": protein_acc in af_proteins,
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

        proteome_id = proteomes.get(protein_acc)
        if proteome_id:
            # Adds proteome
            proteome = proteomesinfo[proteome_id]
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
        pdbe_chains = {}     # mapping PDBe-chain ID -> chain segments
        pdbe_documents = {}  # mapping PDBe-chain ID -> ES document
        for pdbe_id in protein2structures.get(protein_acc, []):
            pdbe_entry = structures[pdbe_id]
            chains = pdbe_entry["proteins"][protein_acc]

            pdbe_doc = doc.copy()
            pdbe_doc.update({
                "structure_acc": pdbe_id.lower(),
                "structure_resolution": pdbe_entry["resolution"],
                "structure_date": pdbe_entry["date"],
                "structure_evidence": pdbe_entry["evidence"],
                "protein_structure": chains,
                "text_structure": join(pdbe_id,
                                       pdbe_entry["evidence"],
                                       pdbe_entry["name"])
            })

            for chain_id, segments in chains.items():
                pdbe_chain_id = f"{pdbe_id}-{chain_id}"

                locations = []
                for segment in segments:
                    locations.append({
                        "fragments": [{
                            "start": segment["protein_start"],
                            "end": segment["protein_end"],
                        }]
                    })

                chain_doc = pdbe_doc.copy()
                chain_doc.update({
                    "structure_chain_acc": chain_id,
                    "structure_protein_locations": locations,
                    "structure_chain": pdbe_chain_id
                })

                pdbe_chains[pdbe_chain_id] = segments
                pdbe_documents[pdbe_chain_id] = chain_doc

        # Adds InterPro entries and member database signatures
        protein_matches = matches.get(protein_acc, {})
        overlapping_chains = set()  # chains associated to at least one entry
        num_rel_docs = 0  # number of relationship documents

        for entry_acc, locations in protein_matches.items():
            used_entries.add(entry_acc)
            entry = entries[entry_acc]
            if entry.integrated_in:
                interpro_acc = entry.integrated_in.lower()
            else:
                interpro_acc = None

            entry_obj = {
                "entry_acc": entry_acc.lower(),
                "entry_db": entry.database,
                "entry_type": entry.type.lower(),
                "entry_date": entry.creation_date.strftime("%Y-%m-%d"),
                "entry_protein_locations": locations,
                "entry_go_terms": [t["identifier"] for t in entry.go_terms],
                "entry_integrated": interpro_acc,
                "text_entry": join(entry_acc, entry.short_name, entry.name,
                                   entry.type.lower(), interpro_acc),
            }

            if entry.clan:
                entry_obj.update({
                    "set_acc": entry.clan["accession"].lower(),
                    "set_db": entry.database,
                    "text_set": join(entry.clan["accession"],
                                     entry.clan["name"]),
                })

            if entry_acc in dom_members:
                entry_obj.update({
                    "ida_id": dom_id,
                    "ida": dom_key,
                })

            # Tests if the entry overlaps PDBe chains
            entry_chains = set()
            for pdbe_chain_id, segments in pdbe_chains.items():
                if overlaps_pdb_chain(locations, segments):
                    # Entry overlaps chain: associate entry to struct/chain
                    chain_doc = pdbe_documents[pdbe_chain_id]
                    entry_doc = chain_doc.copy()
                    entry_doc.update(entry_obj)

                    documents.append((
                        entry.database + version,
                        get_rel_doc_id(entry_doc),
                        entry_doc
                    ))

                    entry_chains.add(pdbe_chain_id)
                    num_rel_docs += 1

            if entry_chains:
                # Stores chains that overlap at least one entry
                overlapping_chains |= entry_chains
            else:
                # Associates the entry to the protein, directly
                entry_doc = doc.copy()
                entry_doc.update(entry_obj)
                documents.append((
                    entry.database + version,
                    get_rel_doc_id(entry_doc),
                    entry_doc
                ))
                num_rel_docs += 1

        # Adds chains that do not overlap with any entry
        for chain_id, chain_doc in pdbe_documents.items():
            if chain_id in overlapping_chains:
                continue

            chain_doc.update({
                "ida_id": dom_id,
                "ida": dom_key,
            })

            documents.append((
                # Not overlapping any entry -> not associated to a member DB
                config.REL_DEFAULT_INDEX_PREFIX + version,
                get_rel_doc_id(chain_doc),
                chain_doc
            ))
            num_rel_docs += 1

        if not num_rel_docs:
            # No relationships for this protein: fallback to protein doc
            documents.append((
                config.REL_DEFAULT_INDEX_PREFIX + version,
                get_rel_doc_id(doc),
                doc
            ))

        while len(documents) >= cache_size:
            for dir_obj in directories:
                filepath = dir_obj.mktemp()
                dumpobj(documents[:cache_size], filepath)
                os.rename(filepath, f"{filepath}{config.EXTENSION}")

            del documents[:cache_size]
            num_documents += cache_size

        if (i + 1) % 10e6 == 0:
            logger.info(f"{i + 1:>15,} proteins, {num_documents:,} documents")

    logger.info(f"{i + 1:>15,} proteins, {num_documents:,} documents")

    proteins.close()
    matches.close()
    proteomes.close()
    domorgs.close()

    logger.info("writing remaining documents")
    # Adds unused entries
    for entry in entries.values():
        if entry.is_public and entry.accession not in used_entries:
            if entry.integrated_in:
                interpro_acc = entry.integrated_in.lower()
            else:
                interpro_acc = None

            doc = init_rel_doc()
            doc.update({
                "entry_acc": entry.accession.lower(),
                "entry_db": entry.database,
                "entry_type": entry.type.lower(),
                "entry_date": entry.creation_date.strftime("%Y-%m-%d"),
                "entry_protein_locations": [],
                "entry_go_terms": [t["identifier"] for t in entry.go_terms],
                "entry_integrated": interpro_acc,
                "text_entry": join(entry.accession, entry.short_name,
                                   entry.name,
                                   entry.type.lower(), interpro_acc),
            })

            if entry.clan:
                doc.update({
                    "set_acc": entry.clan["accession"].lower(),
                    "set_db": entry.database,
                    "text_set": join(entry.clan["accession"],
                                     entry.clan["name"]),
                })

            documents.append((
                entry.database + version,
                get_rel_doc_id(doc),
                doc
            ))

    # Adds unused taxa
    for taxon in taxa.values():
        if taxon["id"] in used_taxa:
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
            config.REL_DEFAULT_INDEX_PREFIX + version,
            get_rel_doc_id(doc),
            doc
        ))

    num_documents += len(documents)
    while documents:
        for dir_obj in directories:
            filepath = dir_obj.mktemp()
            dumpobj(documents[:cache_size], filepath)
            os.rename(filepath, f"{filepath}{config.EXTENSION}")

        del documents[:cache_size]

    # Adds done sentinel file
    for path in outdirs:
        open(os.path.join(path, f"{version}{config.DONE_SUFFIX}"), "w").close()

    logger.info(f"complete: {num_documents:,} documents")
