# -*- coding: utf-8 -*-

import os
import shutil
from typing import Optional, Sequence

from interpro7dw import logger
from interpro7dw.utils import DirectoryTree, Store, datadump, dataload
from interpro7dw.ebi.interpro.staging.database import get_entry_databases
from interpro7dw.ebi.interpro.utils import overlaps_pdb_chain
from . import utils


BODY = {
    "settings": {
        "analysis": {
            "analyzer": {
                "autocomplete": {
                    "tokenizer": "autocomplete",
                    "filter": [
                        "lowercase"
                    ]
                }
            },
            "tokenizer": {
                "autocomplete": {
                    "type": "edge_ngram",
                    "min_gram": 2,
                    "max_gram": 20,
                    "token_chars": [
                        "letter",
                        "digit"
                    ]
                }
            }
        }
    },
    "mappings": {
        "properties": {
            # Protein
            "protein_acc": {"type": "keyword"},
            "protein_length": {"type": "long"},
            "protein_is_fragment": {"type": "keyword"},
            "protein_db": {"type": "keyword"},
            "text_protein": {"type": "text", "analyzer": "autocomplete"},

            # Domain architecture
            "ida_id": {"type": "keyword"},
            "ida": {"type": "keyword"},

            # Taxonomy
            "tax_id": {"type": "long"},
            "tax_name": {"type": "keyword"},
            "tax_lineage": {"type": "keyword"},
            "tax_rank": {"type": "keyword"},
            "text_taxonomy": {"type": "text", "analyzer": "autocomplete"},

            # Proteome
            "proteome_acc": {"type": "keyword"},
            "proteome_name": {"type": "keyword"},
            "proteome_is_reference": {"type": "keyword"},   # todo: remove?
            "text_proteome": {"type": "text", "analyzer": "autocomplete"},

            # Structure
            "structure_acc": {"type": "keyword"},
            "structure_resolution": {"type": "float"},
            "structure_date": {"type": "date"},             # todo: remove?
            "structure_evidence": {"type": "keyword"},
            "protein_structure": {"type": "object", "enabled": False},
            "text_structure": {"type": "text", "analyzer": "autocomplete"},

            # Chain
            "structure_chain_acc": {"type": "text", "analyzer": "keyword"},
            "structure_protein_locations": {"type": "object", "enabled": False},
            "structure_chain": {"type": "text", "analyzer": "keyword", "fielddata": True},

            # Entry
            "entry_acc": {"type": "keyword"},
            "entry_db": {"type": "keyword"},
            "entry_type": {"type": "keyword"},
            "entry_date": {"type": "date"},                 # todo: remove?
            "entry_protein_locations": {"type": "object", "enabled" : False},
            "entry_go_terms": {"type": "keyword"},
            "entry_integrated": {"type": "keyword"},
            "text_entry": {"type": "text", "analyzer": "autocomplete"},

            # Clan/set
            "set_acc": {"type": "keyword"},
            "set_db": {"type": "keyword"},
            "set_integrated": {"type": "keyword"},         # todo: remove?
            "text_set": {"type": "text", "analyzer": "autocomplete"},
        }
    }
}
SHARDS = {
    "cathgene3d": 10,
    "interpro": 10,
    "panther": 10,
    "pfam": 10
}
DEFAULT_INDEX = "others"

# Aliases
STAGING = "rel_staging"
LIVE = "rel_current"
PREVIOUS = "rel_previous"


def join(*args, separator: str=' ') -> str:
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


def init_doc() -> dict:
    return {key: None for key in BODY["mappings"]["properties"].keys()}


def dump_documents(src_proteins: str, src_entries: str,
                   src_proteomes: str, src_structures: str,
                   src_taxonomy: str, src_uniprot2ida: str,
                   src_uniprot2matches: str, src_uniprot2proteomes: str,
                   outdir: str, version: str, cache_size: int=100000):
    logger.info("preparing data")
    try:
        shutil.rmtree(outdir)
    except FileNotFoundError:
        pass
    finally:
        os.makedirs(outdir)
        organizer = DirectoryTree(outdir)
        open(os.path.join(outdir, f"{version}{utils.LOAD_SUFFIX}"), "w").close()

    proteins = Store(src_proteins)
    uniprot2ida = Store(src_uniprot2ida)
    uniprot2matches = Store(src_uniprot2matches)
    uniprot2proteomes = Store(src_uniprot2proteomes)

    entries = dataload(src_entries)  # mem: ~1.5 GB
    proteomes = dataload(src_proteomes)  # mem: <1 GB
    structures = dataload(src_structures)  # mem: ~ 4GB
    taxonomy = dataload(src_taxonomy)  # mem: ~ 2.5GB

    uniprot2pdbe = {}  # mem: <1 GB
    for pdb_id, entry in structures.items():
        for uniprot_acc in entry["proteins"]:
            try:
                uniprot2pdbe[uniprot_acc].append(pdb_id)
            except KeyError:
                uniprot2pdbe[uniprot_acc] = [pdb_id]

    logger.info("starting")
    i = 0
    num_documents = 0
    documents = []
    used_entries = set()
    used_taxa = set()
    for uniprot_acc, info in proteins.items():
        taxid = info["taxid"]

        taxon = taxonomy[taxid]
        used_taxa.add(taxid)  # remember that this taxon has been used

        try:
            dom_members, dom_arch, dom_arch_id = uniprot2ida[uniprot_acc]
        except KeyError:
            dom_members = []
            dom_arch = dom_arch_id = None

        doc = init_doc()
        doc.update({
            "protein_acc": uniprot_acc.lower(),
            "protein_length": info["length"],
            "protein_is_fragment": info["fragment"],
            "protein_db": "reviewed" if info["reviewed"] else "unreviewed",
            "text_protein": join(uniprot_acc, info["identifier"]),

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
            chains = pdb_entry["proteins"][uniprot_acc]

            pdb_doc = doc.copy()
            pdb_doc.update({
                "structure_acc": pdb_id.lower(),
                "structure_resolution": pdb_entry["resolution"],
                "structure_date": pdb_entry["date"],
                "structure_evidence": pdb_entry["evidence"],
                "protein_structure": chains,
                "text_structure": join(pdb_id,
                                       pdb_entry["evidence"],
                                       pdb_entry["name"])
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

                pdb_documents[pdb_chain_id] = chain_doc
                pdb_chains[pdb_chain_id] = segments

        # Adding entries
        overlapping_chains = set()  # chains associated to at least one entry
        matches = uniprot2matches.get(uniprot_acc, {})
        num_protein_docs = 0
        for entry_acc, locations in matches.items():
            used_entries.add(entry_acc)  # this entry has been used

            entry = entries[entry_acc]

            if entry.integrated_in:
                interpro_acc = entry.integrated_in.lower()
            else:
                interpro_acc = None

            entry_obj = {
                "entry_acc": entry_acc.lower(),
                "entry_db": entry.database,
                "entry_type": entry.type,
                "entry_date": entry.creation_date.strftime("%Y-%m-%d"),
                "entry_protein_locations": locations,
                "entry_go_terms": [t["identifier"] for t in entry.go_terms],
                "entry_integrated": interpro_acc,
                "text_entry": join(entry_acc, entry.short_name, entry.name,
                                   entry.type, interpro_acc),
            }

            if entry.clan:
                entry_obj.update({
                    "set_acc": entry.clan["accession"].lower(),
                    "set_db": entry.database,
                    "text_set": join(entry.clan["accession"], entry.clan["name"]),
                })

            if entry_acc in dom_members:
                entry_obj.update({
                    "ida_id": dom_arch_id,
                    "ida": dom_arch,
                })

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
                    num_protein_docs += 1

            if entry_chains:
                # Entry overlaps at least one chain
                overlapping_chains |= entry_chains
            else:
                # Associate entry to protein directly
                entry_doc = doc.copy()
                entry_doc.update(entry_obj)
                documents.append(entry_doc)
                num_protein_docs += 1

        # Add non-overlapping chains
        for chain_id, chain_doc in pdb_documents.items():
            if chain_id in overlapping_chains:
                continue

            chain_doc.update({
                "ida_id": dom_arch_id,
                "ida": dom_arch,
            })

            documents.append(chain_doc)
            num_protein_docs += 1

        if not num_protein_docs:
            # No relationships for this protein: fallback to protein doc
            documents.append(doc)

        while len(documents) >= cache_size:
            filepath = organizer.mktemp()
            datadump(filepath, documents[:cache_size])
            os.rename(filepath, f"{filepath}{utils.EXTENSION}")
            del documents[:cache_size]
            num_documents += cache_size

        i += 1
        if not i % 10000000:
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

        doc = init_doc()
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

        documents.append(doc)

    # Add unused taxa
    for taxon in taxonomy.values():
        if taxon["id"] in used_taxa:
            continue

        doc = init_doc()
        doc.update({
            "tax_id": taxon["id"],
            "tax_name": taxon["full_name"],
            "tax_lineage": taxon["lineage"],
            "tax_rank": taxon["rank"],
            "text_taxonomy": join(taxon["id"], taxon["full_name"],
                                  taxon["rank"])
        })

        documents.append(doc)

    num_documents += len(documents)
    while documents:
        filepath = organizer.mktemp()
        datadump(filepath, documents[:cache_size])
        os.rename(filepath, f"{filepath}.dat")
        del documents[:cache_size]

    proteins.close()
    uniprot2ida.close()
    uniprot2matches.close()
    uniprot2proteomes.close()

    open(os.path.join(outdir, f"{version}{utils.DONE_SUFFIX}"), "w").close()
    os.remove(os.path.join(outdir, f"{version}{utils.LOAD_SUFFIX}"))

    logger.info(f"complete ({num_documents:,} documents)")


def index_documents(url: str, hosts: Sequence[str], indir: str,
                    version: str, outdir: Optional[str]=None,
                    writeback: bool=False, create_indices: bool=True):
    indices = [DEFAULT_INDEX]
    for name in get_entry_databases(url):
        indices.append(name)

    def wrap(doc: dict) -> dict:
        if doc["entry_db"]:
            index = doc["entry_db"] + version
        else:
            index = DEFAULT_INDEX + version

        if doc["protein_acc"]:
            doc_id = join(doc["protein_acc"],
                          doc["proteome_acc"],
                          doc["entry_acc"],
                          doc["set_acc"],
                          doc["structure_acc"],
                          doc["structure_chain_acc"],
                          separator='-')
        elif doc["entry_acc"]:
            doc_id = join(doc["entry_acc"], doc["set_acc"], separator='-')
        else:
            doc_id = doc["tax_id"]

        return {
            "_op_type": "index",
            "_index": index,
            "_id": doc_id,
            "_source": doc
        }

    es = utils.connect(hosts, verbose=False)
    if create_indices:
        logger.info("creating indices")
        for name in indices:
            body = BODY.copy()
            body["settings"].update({
                "index": {
                    # Static settings
                    "number_of_shards": SHARDS.get(name, utils.DEFAULT_SHARDS),
                    "codec": "best_compression",

                    # Dynamic settings
                    "number_of_replicas": 0,  # defaults to 1
                    "refresh_interval": -1  # defaults to 1s
                }
            })

            utils.create_index(es, name + version, body)

    if es.indices.exists_alias(name=PREVIOUS):
        for prev_index in es.indices.get_alias(name=PREVIOUS):
            utils.delete_index(es, prev_index)

    utils.index_documents(es, indir, version, callback=wrap, outdir=outdir,
                          threads=8, writeback=writeback)

    utils.add_alias(es, [idx+version for idx in indices], STAGING)


def publish(hosts: Sequence[str]):
    utils.publish(hosts, STAGING, LIVE, PREVIOUS)
