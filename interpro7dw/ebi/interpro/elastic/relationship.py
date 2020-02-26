# -*- coding: utf-8 -*-

import os
import shutil
from multiprocessing import Process, Queue
from typing import Mapping, Optional, Sequence

from interpro7dw import logger
from interpro7dw.utils import DirectoryTree, Store, datadump, dataload
from interpro7dw.ebi.interpro.utils import overlaps_pdb_chain
from .utils import LOADING, join


def _init_document() -> dict:
    return {
        # Protein
        "protein_acc": None,
        "protein_length": None,
        "protein_is_fragment": None,
        "protein_db": None,
        "text_protein": None,

        # Domain architecture
        "ida_id": None,
        "ida": None,

        # Taxonomy
        "tax_id": None,
        "tax_name": None,
        "tax_lineage": None,
        "tax_rank": None,
        "text_taxonomy": None,

        # Proteome
        "proteome_acc": None,
        "proteome_name": None,
        "proteome_is_reference": None,      # todo: remove?
        "text_proteome": None,

        # Structure
        "structure_acc": None,
        "structure_resolution": None,
        "structure_date": None,             # todo: remove?
        "structure_evidence": None,
        "text_structure": None,

        # Chain
        "structure_chain_acc": None,
        "structure_protein_locations": None,
        "structure_chain": None,

        # Entry
        "entry_acc": None,
        "entry_db": None,
        "entry_type": None,
        "entry_date": None,                 # todo: remove?
        "entry_protein_locations": None,
        "entry_go_terms": None,
        "entry_integrated": None,
        "text_entry": None,

        # Set
        "set_acc": None,
        "set_db": None,
        "set_integrated": None,             # todo: remove?
        "text_set": None,
    }


def gen_prot_rels(uniprot_acc: str, info: Mapping, taxon: Mapping, proteome: Optional[Mapping], structures: Sequence[Mapping], matches: Mapping):
    doc = _init_document()
    doc.update({
        "protein_acc": uniprot_acc.lower(),
        "protein_length": info["length"],
        "protein_is_fragment": info["fragment"],
        "protein_db": "reviewed" if info["reviewed"] else "unreviewed",
        "text_protein": join(uniprot_acc, info["identifier"]),

        # Taxonomy
        "tax_id": taxon["id"],
        "tax_name": taxon["full_name"],
        "tax_lineage": taxon["lineage"],
        "tax_rank": taxon["rank"],
        "text_taxonomy": join(taxon["id"], taxon["full_name"], taxon["rank"])
    })

    if proteome:
        doc.update({
            "proteome_acc": proteome["id"].lower(),
            "proteome_name": proteome["name"],
            "proteome_is_reference": proteome["is_reference"],
            "text_proteome": join(proteome["id"],
                                  proteome["name"],
                                  proteome["assembly"],
                                  proteome["taxon_id"],
                                  proteome["strain"]),
        })

    # Adding PDBe structures/chains
    pdb_chains = {}  # mapping PDB-chain ID -> chain segments
    pdb_documents = {}  # mapping PDB-chain ID -> ES document
    for pdb_entry in structures:
        pdb_doc = doc.copy()
        pdb_doc.update({
            "structure_acc": pdb_entry["id"].lower(),
            "structure_resolution": pdb_entry["resolution"],
            "structure_date": pdb_entry["date"],
            "structure_evidence": pdb_entry["evidence"],
            "text_structure": join(pdb_entry["id"],
                                   pdb_entry["evidence"],
                                   pdb_entry["name"])
        })

        chains = pdb_entry["proteins"][uniprot_acc]
        for chain_id, segments in chains.items():
            pdb_chain_id = f"{pdb_entry['id']}-{chain_id}"

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

            pdb_documents[pdb_chain_id] = chain_doc
            pdb_chains[pdb_chain_id] = segments

    # Adding entries
    documents = []
    overlapping_chains = set()  # chains associated to at least one entry
    for entry_acc, locations in matches.items():
        entry = entries[entry_acc]

        if entry.integrated_in:
            interpro_acc = entry.integrated_in.lower()
        else:
            interpro_acc = None

        go_terms = [t["identifier"] for t in entry.go_terms]

        entry_obj = {
            "entry_acc": entry_acc.lower(),
            "entry_db": entry.database,
            "entry_type": entry.type,
            "entry_date": entry.creation_date.strftime("%Y-%m-%d"),
            "entry_protein_locations": locations,
            "entry_go_terms": go_terms,
            "entry_integrated": interpro_acc,
            "text_entry": join(entry_acc, entry.short_name, entry.name,
                               entry.type, interpro_acc),
        }

        # Test if the entry overlaps PDB chains
        entry_chains = set()
        for pdb_chain_id, segments in pdb_chains.items():
            if overlaps_pdb_chain(locations, segments):
                # Entry overlaps chain: associate entry to struct/chain
                chain_doc = pdb_documents[pdb_chain_id]
                entry_doc = chain_doc.copy()
                entry_doc.update(entry_obj)
                documents.append(entry_doc)
                entry_chains.add(pdb_chain_id)

        if entry_chains:
            # Entry overlaps at least one chain
            overlapping_chains |= entry_chains
        else:
            # Associate entry to protein directly
            entry_doc = doc.copy()
            entry_doc.update(entry_obj)
            documents.append(entry_doc)

    # Add non-overlapping chains
    for chain_id, chain_doc in pdb_documents.items():
        if chain_id in overlapping_chains:
            continue

        documents.append(chain_doc)


def _gen_documents(src_entries: str, src_proteomes: str, src_structures: str,
                   src_taxonomy: str, inqueue: Queue, dir: str):
    organizer = DirectoryTree(root=dir)

    entries = dataload(src_entries)
    proteomes = dataload(src_proteomes)
    structures = dataload(src_structures)
    taxonomy = dataload(src_taxonomy)

    for chunk_type, chunk in iter(inqueue.get, None):
        pass


def dump_entry_documents(src_proteins: str, src_entries: str,
                         src_proteomes: str, src_structures: str,
                         src_taxonomy: str, src_uniprot2ida: str,
                         src_uniprot2matches: str, src_uniprot2proteomes: str,
                         outdir: str, cache_size: int=10000):
    try:
        shutil.rmtree(outdir)
    except FileNotFoundError:
        pass
    finally:
        os.makedirs(outdir)
        organizer = DirectoryTree(outdir)
        open(os.path.join(outdir, LOADING), "w").close()

    proteins = Store(src_proteins)
    uniprot2ida = Store(src_uniprot2ida)
    uniprot2matches = Store(src_uniprot2matches)
    uniprot2proteomes = Store(src_uniprot2proteomes)

    entries = dataload(src_entries)
    proteomes = dataload(src_proteomes)
    structures = dataload(src_structures)
    taxonomy = dataload(src_taxonomy)

    uniprot2pdbe = {}
    for pdb_id, entry in structures.items():
        for uniprot_acc in entry["proteins"]:
            try:
                uniprot2pdbe[uniprot_acc].append(pdb_id)
            except KeyError:
                uniprot2pdbe[uniprot_acc] = [pdb_id]

    i = 0
    num_documents = 0
    cached_documents = []
    used_entries = set()
    used_taxa = set()
    for uniprot_acc, info in proteins.items():
        taxid = info["taxid"]

        # todo: raise Exception instead of skipping proteins
        try:
            taxon = taxonomy[taxid]
        except KeyError:
            continue
            # table.close()
            # con.close()
            # raise RuntimeError(f"{accession}: invalid taxon {taxid}")

        used_taxa.add(taxid)  # remember that this taxon has been used

        try:
            dom_arch, dom_arch_id = uniprot2ida[uniprot_acc]
        except KeyError:
            dom_arch = dom_arch_id = None

        doc = _init_document()
        doc.update({
            "protein_acc": uniprot_acc.lower(),
            "protein_length": info["length"],
            "protein_is_fragment": info["fragment"],
            "protein_db": "reviewed" if info["reviewed"] else "unreviewed",
            "text_protein": join(uniprot_acc, info["identifier"]),

            "ida_id": dom_arch_id,
            "ida": dom_arch,

            # Taxonomy
            "tax_id": taxid,
            "tax_name": taxon["full_name"],
            "tax_lineage": taxon["lineage"],
            "tax_rank": taxon["rank"],
            "text_taxonomy": join(taxid, taxon["full_name"], taxon["rank"])
        })

        proteome_id = uniprot2proteomes.get(uniprot_acc)
        if proteome_id:
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

        # Adding PDBe structures/chains
        pdb_chains = {}     # mapping PDB-chain ID -> chain segments
        pdb_documents = {}  # mapping PDB-chain ID -> ES document
        for pdb_id in uniprot2pdbe.get(uniprot_acc, []):
            pdb_entry = structures[pdb_id]

            pdb_doc = doc.copy()
            pdb_doc.update({
                "structure_acc": pdb_id.lower(),
                "structure_resolution": pdb_entry["resolution"],
                "structure_date": pdb_entry["date"],
                "structure_evidence": pdb_entry["evidence"],
                "text_structure": join(pdb_id,
                                       pdb_entry["evidence"],
                                       pdb_entry["name"])
            })

            chains = pdb_entry["proteins"][uniprot_acc]
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

                pdb_documents[pdb_chain_id] = chain_doc
                pdb_chains[pdb_chain_id] = segments

        # Adding entries
        documents = []
        overlapping_chains = set()  # chains associated to at least one entry
        matches = uniprot2matches.get(uniprot_acc, {})
        for entry_acc, locations in matches.items():
            used_entries.add(entry_acc)  # this entry has been used

            entry = entries[entry_acc]

            if entry.integrated_in:
                interpro_acc = entry.integrated_in.lower()
            else:
                interpro_acc = None

            entry_obj = {
                # We need entry_acc untouched for clans
                "entry_acc": entry_acc,
                "entry_db": entry.database,
                "entry_type": entry.type,
                "entry_date": entry.creation_date.strftime("%Y-%m-%d"),
                "entry_protein_locations": locations,
                "entry_go_terms": [t["identifier"] for t in entry.go_terms],
                "entry_integrated": interpro_acc,
                "text_entry": join(entry_acc, entry.short_name, entry.name,
                                   entry.type, interpro_acc),
            }

            # Test if the entry overlaps PDB chains
            entry_chains = set()
            for pdb_chain_id, segments in pdb_chains.items():
                if overlaps_pdb_chain(locations, segments):
                    # Entry overlaps chain: associate entry to struct/chain
                    chain_doc = pdb_documents[pdb_chain_id]
                    entry_doc = chain_doc.copy()
                    entry_doc.update(entry_obj)
                    documents.append(entry_doc)
                    entry_chains.add(pdb_chain_id)

            if entry_chains:
                # Entry overlaps at least one chain
                overlapping_chains |= entry_chains
            else:
                # Associate entry to protein directly
                entry_doc = doc.copy()
                entry_doc.update(entry_obj)
                documents.append(entry_doc)

        # Add non-overlapping chains
        for chain_id, chain_doc in pdb_documents.items():
            if chain_id in overlapping_chains:
                continue

            documents.append(chain_doc)

        if documents:
            # Add clans in documents with an entry
            for entry_doc in documents:
                entry_acc = entry_doc["entry_acc"]

                if entry_acc:
                    entry = entries[entry_acc]
                    if entry.clan:
                        entry_doc.update({
                            "entry_acc": entry_acc.lower(),
                            "set_acc": entry.clan["accession"].lower(),
                            "set_db": entry.database,
                            "text_set": join(entry.clan["accession"],
                                             entry.clan["name"]),
                        })
                    else:
                        entry_doc["entry_acc"] = entry_acc.lower()

                cached_documents.append(entry_doc)
        else:
            # No relationships for this protein: fallback to protein doc
            cached_documents.append(doc)

        while len(cached_documents) >= cache_size:
            filepath = organizer.mktemp()
            datadump(filepath, cached_documents[:cache_size])
            os.rename(filepath, f"{filepath}.dat")
            del cached_documents[:cache_size]
            num_documents += cache_size

        i += 1
        if not i % 1000:
            logger.info(f"{i:>12,}")

    logger.info(f"{i:>12,}")

    # Add unused entries
    for entry in entries.values():
        if entry.accession in used_entries or entry.is_deleted:
            continue

        if entry.integrated_in:
            interpro_acc = entry.integrated_in.lower()
        else:
            interpro_acc = None

        doc = _init_document()
        doc.update({
            "entry_acc": entry.accession.lower(),
            "entry_db": entry.database,
            "entry_type": entry.type,
            "entry_date": entry.creation_date.strftime("%Y-%m-%d"),
            "entry_protein_locations": [],
            "entry_go_terms": [t["identifier"] for t in entry.go_terms],
            "entry_integrated": interpro_acc,
            "text_entry": join(entry.accession, entry.short_name, entry.name,
                               entry.type, interpro_acc),
        })

        if entry.clan:
            doc.update({
                "set_acc": entry.clan["accession"].lower(),
                "set_db": entry.database,
                "text_set": join(entry.clan["accession"],
                                 entry.clan["name"]),
            })

        cached_documents.append(doc)

    # Add unused taxa
    for taxon in taxonomy.values():
        if taxon["id"] in used_taxa:
            continue

        doc = _init_document()
        doc.update({
            "tax_id": taxon["id"],
            "tax_name": taxon["full_name"],
            "tax_lineage": taxon["lineage"],
            "tax_rank": taxon["rank"],
            "text_taxonomy": join(taxon["id"], taxon["full_name"],
                                  taxon["rank"])
        })

        cached_documents.append(doc)

    num_documents += len(cached_documents)
    while cached_documents:
        filepath = organizer.mktemp()
        datadump(filepath, cached_documents[:cache_size])
        os.rename(filepath, f"{filepath}.dat")
        del cached_documents[:cache_size]

    proteins.close()
    uniprot2ida.close()
    uniprot2matches.close()
    uniprot2proteomes.close()

    # Delete flag file to notify loaders that all files are ready
    os.remove(os.path.join(outdir, LOADING))

    logger.info(f"complete ({num_documents:,} documents)")
