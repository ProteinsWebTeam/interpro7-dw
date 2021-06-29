# -*- coding: utf-8 -*-

import glob
import logging
import os
import shutil
import time
from typing import Sequence

from elasticsearch import Elasticsearch, exceptions
from elasticsearch.helpers import parallel_bulk as pbulk

from interpro7dw import logger
from interpro7dw.ebi.interpro.staging.database import get_entry_databases
from interpro7dw.ebi.interpro.utils import overlaps_pdb_chain
from interpro7dw.ebi.interpro.utils import parse_uniprot_struct_models
from interpro7dw.utils import DirectoryTree, Store, dumpobj, loadobj


IDA_BODY = {
    "mappings": {
        "properties": {
            "ida_id": {"type": "keyword"},
            "ida": {"type": "keyword"},
            "counts": {"type": "integer"}
        }
    }
}
IDA_INDEX = "ida"
REL_BODY = {
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
            "protein_has_model": {"type": "keyword"},
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
            "proteome_is_reference": {"type": "keyword"},
            "text_proteome": {"type": "text", "analyzer": "autocomplete"},

            # Structure
            "structure_acc": {"type": "keyword"},
            "structure_resolution": {"type": "float"},
            "structure_date": {"type": "date"},
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
            "entry_date": {"type": "date"},
            "entry_protein_locations": {"type": "object", "enabled": False},
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
REL_SHARDS = {
    "cathgene3d": 10,
    "interpro": 10,
    "panther": 10,
    "pfam": 10
}
REL_INDEX = "others"
DEFAULT_SHARDS = 5

IDA_BASE_ALIAS = "ida"
REL_BASE_ALIAS = "rel"
STAGING_ALIAS_SUFFIX = "_staging"
LIVE_ALIAS_SUFFIX = "_current"
PREVIOUS_ALIAS_SUFFIX = "_previous"

EXTENSION = ".dat"
LOAD_SUFFIX = ".load"
DONE_SUFFIX = ".done"


def add_alias(es: Elasticsearch, indices: Sequence[str], alias: str):
    if es.indices.exists_alias(name=alias):
        # Alias already exists: update it

        # Indices currently pointed by the alias
        old_indices = set(es.indices.get_alias(name=alias))

        actions = []
        for index in set(indices):
            if index in old_indices:
                # Index is already pointed by alias
                old_indices.remove(index)
            else:
                # Add alias to index
                actions.append({"add": {"index": index, "alias": alias}})

        # Remove alias from old indices
        for index in old_indices:
            actions.append({"remove": {"index": index, "alias": alias}})

        if actions:
            """
            Atomic operation:
            Alias removed from the old indices
            at the same time it's added to the new ones
            """
            es.indices.update_aliases(body={"actions": actions})
    else:
        # Creat new alias, then point to indices
        es.indices.put_alias(index=','.join(indices), name=alias)


def connect(hosts: Sequence[str], timeout: int = 10,
            verbose: bool = True) -> Elasticsearch:
    es = Elasticsearch(hosts=hosts, timeout=timeout)

    if not verbose:
        # Disable Elastic logger
        tracer = logging.getLogger("elasticsearch")
        tracer.setLevel(logging.CRITICAL + 1)

    return es


def create_indices(url: str, hosts: Sequence[str], version: str):
    es = connect(hosts, verbose=False)

    """
    Assuming we are creating indices for version 100.0,
    version 99.0 indices have the 'live' alias now,
    and will have the 'previous' alias when we make v100.0 live.
    However, right now, 'previous' should point to v98.0 (to be deleted)
    """
    for base_alias in (IDA_BASE_ALIAS, REL_BASE_ALIAS):
        alias = base_alias + PREVIOUS_ALIAS_SUFFIX
        if es.indices.exists_alias(name=alias):
            for idx in es.indices.get_alias(name=alias):
                delete_index(es, idx)

    # Create a list of all new indices to create
    indices = [
        (IDA_INDEX, IDA_BODY, IDA_BASE_ALIAS),
        (REL_INDEX, REL_BODY, REL_BASE_ALIAS)
    ]
    for name in get_entry_databases(url):
        indices.append((name, REL_BODY, REL_BASE_ALIAS))

    # Create new indices
    alias2indices = {}
    for index, body, base_alias in indices:
        try:
            settings = body["settings"]
        except KeyError:
            settings = body["settings"] = {}

        settings.update({
            "index": {
                # Static settings
                "number_of_shards": REL_SHARDS.get(index, DEFAULT_SHARDS),
                "codec": "best_compression",

                # Dynamic settings
                "number_of_replicas": 0,  # defaults to 1
                "refresh_interval": -1  # defaults to 1s
            }
        })

        index += version  # Use InterPro version as suffix

        delete_index(es, index)
        while True:
            try:
                es.indices.create(index=index, body=body)
            except exceptions.ConnectionTimeout:
                time.sleep(30)
            except exceptions.RequestError as exc:
                raise exc
            except Exception as exc:
                logger.warning(f"{index}: {exc}")
                time.sleep(30)
            else:
                break

        try:
            alias2indices[base_alias].append(index)
        except KeyError:
            alias2indices[base_alias] = [index]

    # Add an 'staging' alias to all newly created indices
    for base_alias, new_indices in alias2indices.items():
        add_alias(es, new_indices, base_alias + STAGING_ALIAS_SUFFIX)


def delete_index(es: Elasticsearch, name: str):
    while True:
        try:
            es.indices.delete(index=name)
        except exceptions.ConnectionTimeout:
            time.sleep(30)
        except exceptions.NotFoundError:
            break
        except Exception as exc:
            logger.error(f"{name}: {exc}")
            time.sleep(30)
        else:
            break


def join(*args, separator: str = ' ') -> str:
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


def init_rel_doc() -> dict:
    return {key: None for key in REL_BODY["mappings"]["properties"].keys()}


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
    else:
        return doc["tax_id"]


def export_documents(src_proteins: str, src_entries: str, src_proteomes: str,
                     src_structures: str, src_taxonomy: str,
                     src_uniprot2ida: str, src_uniprot2matches: str,
                     src_uniprot2proteomes: str, src_uniprot_models: str,
                     outdirs: Sequence[str], version: str,
                     cache_size: int = 100000):
    logger.info("preparing data")
    os.umask(0o002)
    organizers = []
    for path in outdirs:
        try:
            shutil.rmtree(path)
        except FileNotFoundError:
            pass

        os.makedirs(path, mode=0o775)
        organizers.append(DirectoryTree(path))
        open(os.path.join(path, f"{version}{LOAD_SUFFIX}"), "w").close()

    uniprot_models = set()
    if src_uniprot_models:
        logger.info("loading UniProt entries with structural model")
        uniprot_models = parse_uniprot_struct_models(src_uniprot_models)

    logger.info("loading domain architectures")
    domains = {}
    with Store(src_uniprot2ida) as u2ida:
        for dom_members, dom_arch, dom_arch_id in u2ida.values():
            try:
                dom = domains[dom_arch_id]
            except KeyError:
                domains[dom_arch_id] = {
                    "ida_id": dom_arch_id,
                    "ida": dom_arch,
                    "counts": 1
                }
            else:
                dom["counts"] += 1

    logger.info("writing IDA documents")
    num_documents = 0
    domains = list(domains.values())
    for i in range(0, len(domains), cache_size):
        documents = []
        for dom in domains[i:i + cache_size]:
            documents.append((
                IDA_INDEX + version,
                dom["ida_id"],
                dom
            ))

        num_documents += len(documents)
        for org in organizers:
            filepath = org.mktemp()
            dumpobj(filepath, documents)
            os.rename(filepath, f"{filepath}{EXTENSION}")

    domains = None

    proteins = Store(src_proteins)
    uniprot2ida = Store(src_uniprot2ida)
    uniprot2matches = Store(src_uniprot2matches)
    uniprot2proteomes = Store(src_uniprot2proteomes)

    entries = loadobj(src_entries)  # mem: ~1.5 GB
    proteomes = loadobj(src_proteomes)  # mem: <1 GB
    structures = loadobj(src_structures)  # mem: ~ 4GB
    taxonomy = loadobj(src_taxonomy)  # mem: ~ 2.5GB

    uniprot2pdbe = {}  # mem: <1 GB
    for pdb_id, entry in structures.items():
        for uniprot_acc in entry["proteins"]:
            try:
                uniprot2pdbe[uniprot_acc].append(pdb_id)
            except KeyError:
                uniprot2pdbe[uniprot_acc] = [pdb_id]

    logger.info("writing relationship documents")
    i = 0
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

        # Create an empty document (all properties set to None)
        doc = init_rel_doc()
        doc.update({
            "protein_acc": uniprot_acc.lower(),
            "protein_length": info["length"],
            "protein_is_fragment": info["fragment"],
            "protein_has_model": uniprot_acc in uniprot_models,
            "protein_db": "reviewed" if info["reviewed"] else "unreviewed",
            "text_protein": join(uniprot_acc, info["identifier"]),

            # Taxonomy
            "tax_id": taxid,
            "tax_name": taxon["sci_name"],
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
        pdb_chains = {}  # mapping PDB-chain ID -> chain segments
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

                    documents.append((
                        entry.database + version,
                        get_rel_doc_id(entry_doc),
                        entry_doc
                    ))

                    entry_chains.add(pdb_chain_id)
                    num_protein_docs += 1

            if entry_chains:
                # Entry overlaps at least one chain
                overlapping_chains |= entry_chains
            else:
                # Associate entry to protein directly
                entry_doc = doc.copy()
                entry_doc.update(entry_obj)
                documents.append((
                    entry.database + version,
                    get_rel_doc_id(entry_doc),
                    entry_doc
                ))
                num_protein_docs += 1

        # Add chains not overlapping any entry
        for chain_id, chain_doc in pdb_documents.items():
            if chain_id in overlapping_chains:
                continue

            chain_doc.update({
                "ida_id": dom_arch_id,
                "ida": dom_arch,
            })

            documents.append((
                # Not overlapping any entry -> not associated to a member DB
                REL_INDEX + version,
                get_rel_doc_id(chain_doc),
                chain_doc
            ))
            num_protein_docs += 1

        if not num_protein_docs:
            # No relationships for this protein: fallback to protein doc
            documents.append((
                REL_INDEX + version,
                get_rel_doc_id(doc),
                doc
            ))

        while len(documents) >= cache_size:
            for org in organizers:
                filepath = org.mktemp()
                dumpobj(filepath, documents[:cache_size])
                os.rename(filepath, f"{filepath}{EXTENSION}")

            del documents[:cache_size]
            num_documents += cache_size

        i += 1
        if not i % 10000000:
            logger.info(f"{i:>12,}")

    logger.info(f"{i:>12,}")

    logger.info("writing remaining documents")
    # Add unused entries
    for entry in entries.values():
        if entry.accession in used_entries or entry.is_deleted:
            continue

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
            "text_entry": join(entry.accession, entry.short_name, entry.name,
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

    # Add unused taxa
    for taxon in taxonomy.values():
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
            REL_INDEX + version,
            get_rel_doc_id(doc),
            doc
        ))

    num_documents += len(documents)
    while documents:
        for org in organizers:
            filepath = org.mktemp()
            dumpobj(filepath, documents[:cache_size])
            os.rename(filepath, f"{filepath}{EXTENSION}")

        del documents[:cache_size]

    proteins.close()
    uniprot2ida.close()
    uniprot2matches.close()
    uniprot2proteomes.close()

    for path in outdirs:
        open(os.path.join(path, f"{version}{DONE_SUFFIX}"), "w").close()

    logger.info(f"complete ({num_documents:,} documents)")


def iter_files(root: str, version: str):
    load_sentinel = os.path.join(root, f"{version}{LOAD_SUFFIX}")
    done_sentinel = os.path.join(root, f"{version}{DONE_SUFFIX}")

    if not os.path.isfile(done_sentinel):
        while not os.path.isfile(load_sentinel):
            # Wait until files start being generated
            time.sleep(60)

    pathname = os.path.join(root, "**", f"*{EXTENSION}")
    files = set()
    active = True
    while True:
        for path in glob.iglob(pathname, recursive=True):
            if path in files:
                continue

            files.add(path)
            yield path

        if not active:
            break
        elif os.path.isfile(done_sentinel):
            # All files ready: they will all be found in the next iteration
            active = False
        else:
            # Files are still being written
            time.sleep(60)


def index_documents(hosts: Sequence[str], indir: str, version: str,
                    threads: int = 4, step: int = 100e6):
    kwargs = {
        "thread_count": threads,
        "queue_size": threads,
        "raise_on_exception": False,
        "raise_on_error": False
    }

    es = connect(hosts, timeout=30, verbose=False)
    num_documents = 0
    num_indexed = 0
    first_pass = True
    while True:
        for filepath in iter_files(indir, version):
            docs = loadobj(filepath)

            if first_pass:
                # Count only once the number of documents to index
                num_documents += len(docs)

            actions = []
            for idx, doc_id, doc in docs:
                actions.append({
                    "_op_type": "index",
                    "_index": idx,
                    "_id": doc_id,
                    "_source": doc
                })

            failed = []
            for i, (ok, info) in enumerate(pbulk(es, actions, **kwargs)):
                if ok:
                    num_indexed += 1
                    if not num_indexed % 100e6:
                        logger.info(f"{num_indexed:>14,} / {num_documents:,}")
                else:
                    failed.append(docs[i])

                    # try:
                    #     is_429 = info["index"]["status"] == 429
                    # except (KeyError, IndexError):
                    #     is_429 = False
                    #
                    # try:
                    #     exc = info["index"]["exception"]
                    # except (KeyError, TypeError):
                    #     exc = None
                    #
                    # if is_429 or isinstance(exc, exceptions.ConnectionTimeout):
                    #     pause = True
                    # else:
                    #     logger.debug(info)

            if failed:
                # Overwrite file with failed documents
                dumpobj(filepath, failed)
            else:
                # Remove file as all documents have been successfully indexed
                os.remove(filepath)

        logger.info(f"{num_indexed:>14,} / {num_documents:,}")
        first_pass = False

        if num_indexed == num_documents:
            break

    # Update index settings
    for base_alias in (IDA_BASE_ALIAS, REL_BASE_ALIAS):
        alias = base_alias + STAGING_ALIAS_SUFFIX

        # This assumes there are indices with the 'staging' alias
        for index in es.indices.get_alias(name=alias):
            es.indices.put_settings(
                body={
                    "number_of_replicas": 1,
                    "refresh_interval": None  # default (1s)
                },
                index=index
            )


def publish(hosts: Sequence[str]):
    es = connect(hosts, verbose=False)

    for base in (IDA_BASE_ALIAS, REL_BASE_ALIAS):
        live_alias = base + LIVE_ALIAS_SUFFIX

        # Add the 'previous' alias to current 'live' indices
        if es.indices.exists_alias(name=live_alias):
            indices = es.indices.get_alias(name=live_alias)

            prev_alias = base + PREVIOUS_ALIAS_SUFFIX
            add_alias(es, indices, prev_alias)

        # Add the 'live' alias to current 'staging' indices
        staging_alias = base + STAGING_ALIAS_SUFFIX
        indices = es.indices.get_alias(name=staging_alias)
        add_alias(es, indices, live_alias)
