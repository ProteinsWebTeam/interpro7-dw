# -*- coding: utf-8 -*-

import glob
import json
import logging
import os
import shutil
import time
from multiprocessing import Process, Queue
from typing import Dict, List, Optional

from elasticsearch import Elasticsearch, exceptions, helpers

from i7dw import io, logger
from i7dw.interpro import DomainArchitecture, mysql


LOADING_FILE = "loading"
NODB_INDEX = "others"
MOBIDBLITE = "mobidblt"
BODY_PATH = os.path.join(os.path.dirname(__file__), "body.json")
SPECIAL_SHARDS = {
    'cathgene3d': 10,
    'interpro': 10,
    'panther': 10,
    'pfam': 10
}
DEFAULT_SHARDS = 5


class DocumentProducer(Process):
    def __init__(self, url: str, task_queue: Queue,
                 done_queue: Queue, outdir: str, min_overlap: int=20):
        super().__init__()
        self.url = url
        self.task_queue = task_queue
        self.done_queue = done_queue
        self.min_overlap = min_overlap
        self.organizer = io.JsonFileOrganizer(outdir)

        self.entries = {}
        self.integrated = {}
        self.entry_set = {}
        self.proteomes = {}
        self.structures = {}
        self.protein_structures = {}
        self.dom_arch = None

    def run(self):
        # Get PDBe structures, entries, sets, and proteomes
        for s in mysql.structures.iter_structures(self.url):
            pdbe_id = s["accession"]
            self.structures[pdbe_id] = s

            for protein_acc in s["proteins"]:
                try:
                    self.protein_structures[protein_acc].add(pdbe_id)
                except KeyError:
                    self.protein_structures[protein_acc] = {pdbe_id}

        self.entries = mysql.entries.get_entries(self.url)
        self.dom_arch = DomainArchitecture(self.entries)

        for s in mysql.entries.iter_sets(self.url):
            members = s.pop("members")
            for entry_acc in members:
                """
                Sets are checked after the initial document has be created.
                As we use lowercase for accessions, we need to use
                lowercase here as well
                """
                self.entry_set[entry_acc.lower()] = s

        for p in mysql.proteomes.iter_proteomes(self.url):
            self.proteomes[p["accession"]] = p

        n_docs = 0
        for doc_type, chunk in iter(self.task_queue.get, None):
            if doc_type == "protein":
                fn = self.process_protein
            elif doc_type == "entry":
                fn = self.process_entry
            elif doc_type == "taxon":
                fn = self.process_taxon
            else:
                raise ValueError(f"'{doc_type}' is not a valid doc type")

            for args in chunk:
                for doc in fn(*args):
                    if doc["entry_db"] != MOBIDBLITE:
                        self.organizer.add(doc)
                        n_docs += 1

        self.organizer.flush()
        self.done_queue.put(n_docs)

    def process_protein(self, accession: str, identifier: str, name: str,
                        database: str, is_fragment: bool, length: int,
                        comments: List, matches: Dict[str, List[Dict]],
                        upid: Optional[str], taxon: Dict) -> List[Dict]:

        self.dom_arch.update(matches)
        protein_doc = self.init_document()

        # Add protein info
        protein_doc.update({
            "protein_acc": accession.lower(),
            "protein_length": length,
            "protein_is_fragment": is_fragment,
            "protein_db": database,
            "text_protein": joinitems(accession, identifier, name, database,
                                      comments, taxon["full_name"]),

            "tax_id": taxon["id"],
            "tax_name": taxon["scientific_name"],
            "tax_lineage": taxon["lineage"],
            "tax_rank": taxon["rank"],
            "text_taxonomy": joinitems(taxon["id"], taxon["full_name"],
                                       taxon["rank"]),

            "ida": self.dom_arch.accession,
            "ida_id": self.dom_arch.identifier
        })

        # Add proteome, in any
        if upid:
            p = self.proteomes[upid]
            protein_doc.update({
                "proteome_acc": upid.lower(),
                "proteome_name": p["name"],
                "proteome_is_reference": p["is_reference"],
                "text_proteome": joinitems(*list(p.values()))
            })

        # Adding PDBe structures and chains
        chain_documents = {}
        struct_chains = {}
        for pdbe_id in self.protein_structures.get(accession, []):
            s = self.structures[pdbe_id]
            text = joinitems(pdbe_id, s["evidence"], s["name"],
                             *(pub["title"]
                               for pub in s["citations"].values()
                               if pub.get("title")))

            chains = s["proteins"][accession]

            struct_doc = protein_doc.copy()
            struct_doc.update({
                "structure_acc": pdbe_id.lower(),
                "structure_resolution": s["resolution"],
                "structure_date": s["date"].strftime("%Y-%m-%d"),
                "structure_evidence": s["evidence"],
                "protein_structure": chains
            })

            for chain_id, chain_fragments in chains.items():
                chain_doc = struct_doc.copy()
                struct_chain_id = pdbe_id + '-' + chain_id
                chain_doc.update({
                    "structure_chain_acc": chain_id,
                    "structure_protein_locations": [
                        {
                            "fragments": [
                                {
                                    "start": fragment["protein_start"],
                                    "end": fragment["protein_end"]
                                }
                            ]
                        } for fragment in chain_fragments
                    ],
                    "structure_chain": struct_chain_id,
                    "text_structure": f"{chain_id} {text}"
                })
                chain_documents[struct_chain_id] = chain_doc
                struct_chains[struct_chain_id] = chains

        # Adding entries
        documents = []
        overlapping_chains = set()
        for entry_acc, locations in matches.items():
            entry = self.entries[entry_acc]

            if entry["integrated"]:
                integrated_acc = entry["integrated"].lower()
            else:
                integrated_acc = None

            go_terms = [t["identifier"] for t in entry["go_terms"]]

            entry_obj = {
                "entry_acc": entry["accession"].lower(),
                "entry_db": entry["database"],
                "entry_type": entry["type"],
                "entry_date": entry["date"].strftime("%Y-%m-%d"),
                "entry_integrated": integrated_acc,
                "text_entry": joinitems(entry["accession"], entry["name"],
                                        entry["type"], entry["descriptions"],
                                        *go_terms),
                "entry_protein_locations": locations,
                "entry_go_terms": go_terms
            }

            # Associate entry to structure/chain if they overlap
            overlaps = False
            for struct_chain_id, chain_doc in chain_documents.items():
                chains = struct_chains[struct_chain_id]
                if mysql.entries.overlaps_with_structure(locations, chains):
                    entry_doc = chain_doc.copy()
                    entry_doc.update(entry_obj)
                    documents.append(entry_doc)
                    overlaps = True
                    overlapping_chains.add(struct_chain_id)

            if not overlaps:
                # Entry does not overlap with ANY chain: associate with protein
                entry_doc = protein_doc.copy()
                entry_doc.update(entry_obj)
                documents.append(entry_doc)

        # Associated not used chains with protein
        for struct_chain_id in set(chain_documents) - overlapping_chains:
            documents.append(chain_documents[struct_chain_id])

        # Finally add sets
        for entry_doc in documents:
            try:
                entry_set = self.entry_set[entry_doc["entry_acc"]]
            except KeyError:
                pass
            else:
                entry_doc.update({
                    "set_acc": entry_set["accession"].lower(),
                    "set_db": entry_set["database"],
                    "text_set": joinitems(*list(entry_set.values())),
                    # todo: implement set integration (e.g. pathways)
                    "set_integrated": []
                })

        return documents if documents else [protein_doc]

    def process_entry(self, accession: str) -> List[dict]:
        entry = self.entries[accession]
        entry_doc = self.init_document()

        if entry["integrated"]:
            integrated_acc = entry["integrated"].lower()
        else:
            integrated_acc = None

        go_terms = [t["identifier"] for t in entry["go_terms"]]

        entry_doc.update({
            "entry_acc": entry["accession"].lower(),
            "entry_db": entry["database"],
            "entry_type": entry["type"],
            "entry_date": entry["date"].strftime("%Y-%m-%d"),
            "entry_integrated": integrated_acc,
            "text_entry": joinitems(entry["accession"], entry["name"],
                                    entry["type"], entry["descriptions"],
                                    *go_terms),
            "entry_protein_locations": [],
            "entry_go_terms": go_terms
        })

        try:
            entry_set = self.entry_set[entry_doc["entry_acc"]]
        except KeyError:
            pass
        else:
            entry_doc.update({
                "set_acc": entry_set["accession"].lower(),
                "set_db": entry_set["database"],
                "text_set": joinitems(*list(entry_set.values())),
                # todo: implement set integration (e.g. pathways)
                "set_integrated": []
            })

        return [entry_doc]

    def process_taxon(self, taxon: dict) -> List[dict]:
        tax_doc = self.init_document()
        tax_doc.update({
            "tax_id": taxon["id"],
            "tax_name": taxon["scientific_name"],
            "tax_lineage": taxon["lineage"],
            "tax_rank": taxon["rank"],
            "text_taxonomy": joinitems(taxon["id"], taxon["full_name"],
                                       taxon["rank"]),
        })
        return [tax_doc]

    @staticmethod
    def init_document() -> dict:
        return {
            # Protein
            "protein_acc": None,
            "protein_length": None,
            "protein_is_fragment": None,
            "protein_db": None,
            "text_protein": None,

            # Taxonomy
            "tax_id": None,
            "tax_name": None,
            "tax_lineage": None,
            "tax_rank": None,
            "text_taxonomy": None,

            # Proteome
            "proteome_acc": None,
            "proteome_name": None,
            "proteome_is_reference": None,
            "text_proteome": None,

            # Entry
            "entry_acc": None,
            "entry_db": None,
            "entry_type": None,
            "entry_date": None,
            "entry_protein_locations": None,
            "entry_go_terms": None,
            "entry_integrated": None,
            "text_entry": None,

            # Set
            "set_acc": None,
            "set_db": None,
            "set_integrated": None,
            "text_set": None,

            # Structure
            "structure_acc": None,
            "structure_resolution": None,
            "structure_date": None,
            "structure_evidence": None,

            # Chain
            "structure_chain_acc": None,
            "structure_protein_locations": None,
            "structure_chain": None,
            "text_structure": None,

            # Domain architecture
            "ida_id": None,
            "ida": None
        }


def joinitems(*args, separator: str=' ') -> str:
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


def create_indices(hosts: List[str], indices: List[str], suffix: str=''):
    # Load default settings and property mapping
    with open(BODY_PATH, "rt") as fh:
        body = json.load(fh)

    # Establish connection
    es = Elasticsearch(hosts, timeout=30)

    # Disable Elastic logger
    tracer = logging.getLogger("elasticsearch")
    tracer.setLevel(logging.CRITICAL + 1)

    # Create indices
    for index in indices:
        """
        Change settings for large bulk imports:

        https://www.elastic.co/guide/en/elasticsearch/reference/current/indices-update-settings.html
        https://www.elastic.co/guide/en/elasticsearch/guide/current/indexing-performance.html
        """
        body["settings"].update({
            "index": {
                # Static settings
                "number_of_shards": SPECIAL_SHARDS.get(index, DEFAULT_SHARDS),
                "codec": "best_compression",

                # Dynamic settings
                "number_of_replicas": 0,  # defaults to 1
                "refresh_interval": -1    # defaults to 1s
            }
        })

        index += suffix

        # Make sure the index is deleted
        while True:
            try:
                es.indices.delete(index)
            except exceptions.NotFoundError:
                break
            except Exception as exc:
                logger.error(f"{type(exc)}: {exc}")
                time.sleep(10)
            else:
                break

        # And make sure it's created
        while True:
            try:
                es.indices.create(index, body=body)
            except exceptions.RequestError as exc:
                raise exc
            except Exception as exc:
                logger.error(f"{type(exc)}: {exc}")
                time.sleep(10)
            else:
                break


def iter_json_files(root: str, seconds: int=60):
    flag_file = os.path.join(root, LOADING_FILE)
    pathname = os.path.join(root, "**", "*.json")
    files = set()

    active = True
    while True:
        for path in glob.iglob(pathname, recursive=True):
            if path not in files:
                files.add(path)
                yield path

        if os.path.isfile(flag_file):
            time.sleep(seconds)
        elif active:
            # All files ready, but loop one last time
            active = False
        else:
            break


def postprocess(task_queue: Queue, done_queue: Queue,
                outdir: Optional[str]=None, write_back: bool=False):
    if outdir:
        """
        Ensure the directory does not exist
        as we don't want files from a previous run to be considered
        """
        try:
            shutil.rmtree(outdir)
        except FileNotFoundError:
            pass
        finally:
            os.makedirs(outdir)
            organizer = io.JsonFileOrganizer(outdir)
    else:
        organizer = None

    milestone_step = 100e6
    while True:
        num_indexed = 0
        num_failed = 0
        num_docs = 0
        next_milestone = milestone_step

        for path, num_successful, errors in iter(task_queue.get, None):
            num_indexed += num_successful
            num_failed += len(errors)

            if errors:
                if organizer:
                    for doc in errors:
                        organizer.add(doc)
                elif write_back:
                    with open(path, "wt") as fh:
                        json.dump(errors, fh)
            elif write_back:
                os.remove(path)

            num_docs = num_indexed + num_failed
            if num_docs >= next_milestone:
                next_milestone += milestone_step
                logger.info(f"{num_docs:>15,} "
                            f"(success: {num_indexed/num_docs*100:>5.1f}%)")

        if organizer:
            organizer.flush()
        logger.info(f"{num_docs:>15,} "
                    f"(success: {num_indexed/num_docs*100:>5.1f}%)")

        done_queue.put(num_failed)

        """
        Wait for instruction from parent:
            - True: we continue
            - False: max retries reached: quit
        """
        if not task_queue.get():
            break


class DocumentLoader(Process):
    def __init__(self, hosts: List, task_queue: Queue, done_queue: Queue,
                 chunk_size: int=500, max_bytes: int=100*1024*1024,
                 suffix: str='', threads: int=4, timeout: int=4):
        super().__init__()
        self.hosts = hosts
        self.task_queue = task_queue
        self.done_queue = done_queue
        self.chunk_size = chunk_size
        self.max_bytes = max_bytes
        self.suffix = suffix
        self.threads = threads
        self.timeout = timeout

    def run(self):
        es = Elasticsearch(self.hosts, timeout=self.timeout)

        # Disable Elastic logger
        tracer = logging.getLogger("elasticsearch")
        tracer.setLevel(logging.CRITICAL + 1)

        for filepath in iter(self.task_queue.get, None):
            with open(filepath, "rt") as fh:
                documents = json.load(fh)

            bulk = helpers.parallel_bulk(
                es, map(self.wrap, documents),
                thread_count=self.threads,
                queue_size=self.threads,
                chunk_size=self.chunk_size,
                max_chunk_bytes=self.max_bytes,
                raise_on_exception=False,
                raise_on_error=False
            )

            num_successful = 0
            errors = []
            for i, (status, item) in enumerate(bulk):
                if status:
                    num_successful += 1
                else:
                    errors.append(documents[i])
                    try:
                        del item["index"]["data"]
                    except KeyError:
                        pass
                    finally:
                        logger.debug(item)

            self.done_queue.put((filepath, num_successful, errors))

    def wrap(self, document: dict) -> dict:
        if document["entry_db"]:
            idx = document["entry_db"] + self.suffix
        else:
            idx = NODB_INDEX + self.suffix

        if document["protein_acc"]:
            doc_id = joinitems(document["protein_acc"],
                               document["proteome_acc"],
                               document["entry_acc"],
                               document["set_acc"],
                               document["structure_acc"],
                               document["structure_chain_acc"],
                               separator='-')
        elif document["entry_acc"]:
            doc_id = joinitems(document["entry_acc"],
                               document["set_acc"],
                               separator='-')
        else:
            doc_id = document["tax_id"]

        return {
            "_op_type": "index",
            "_index": idx,
            "_id": doc_id,
            "_source": document
        }


def index_documents(hosts: List[str], src: str, suffix: str,
                    outdir: Optional[str], max_retries: int, processes: int,
                    write_back: bool) -> int:
    inqueue = Queue()
    outqueue = Queue(maxsize=processes)
    cntqueue = Queue()
    postprocessor = Process(target=postprocess,
                            args=(outqueue, cntqueue, outdir, write_back))
    postprocessor.start()

    processes = max(1, processes-2)  # parent process + postprocessor
    num_attempts = 0
    while True:
        num_attempts += 1
        logger.info(f"indexing documents (try #{num_attempts})")
        workers = []

        for i in range(processes):
            p = DocumentLoader(hosts, inqueue, outqueue, suffix=suffix)
            p.start()
            workers.append(p)

        for path in iter_json_files(src):
            inqueue.put(path)

        # All files enqueued
        for _ in workers:
            inqueue.put(None)

        # Wait for workers to complete
        for p in workers:
            p.join()

        # Inform postprocessor that we want the count of failed documents
        outqueue.put(None)
        num_errors = cntqueue.get()

        if num_errors and num_attempts <= max_retries:
            outqueue.put(True)   # continue
        else:
            outqueue.put(False)  # stop
            postprocessor.join()
            break

    return num_errors


def update_alias(hosts: List[str], indices: List[str], alias: str,
                 suffix: str='', keep_prev_indices: bool=True):
    new_indices = {index + suffix for index in indices}
    es = Elasticsearch(hosts, timeout=30)

    # Disable Elastic logger
    tracer = logging.getLogger("elasticsearch")
    tracer.setLevel(logging.CRITICAL + 1)

    exists = es.indices.exists_alias(name=alias)
    if exists:
        # Alias already exists: update it

        # Indices currently using the alias
        current_indices = set(es.indices.get_alias(name=alias))

        actions = []
        for index in new_indices:
            try:
                # If passes: new index is already using the alias
                current_indices.remove(index)
            except KeyError:
                # Otherwise, add the alias to the new index
                actions.append({
                    "add": {
                        "index": index,
                        "alias": alias
                    }
                })

        # Remove the alias from the current indices
        for index in current_indices:
            actions.append({
                "remove": {
                    "index": index,
                    "alias": alias
                }
            })

        if actions:
            """
            Atomic operation:
            Alias removed from the old indices
            at the same time it's added to the new ones
            """
            es.indices.update_aliases(body={"actions": actions})

        if not keep_prev_indices:
            # Delete old indices that used the alias
            for index in current_indices:
                while True:
                    try:
                        es.indices.delete(index)
                    except exceptions.NotFoundError:
                        break
                    except Exception as exc:
                        logger.error(f"{type(exc)}: {exc}")
                    else:
                        break
    else:
        # Create alias, and add it to new indices
        es.indices.put_alias(index=','.join(new_indices), name=alias)

    # Update index settings
    for index in new_indices:
        es.indices.put_settings({
            # "number_of_replicas": 1,
            "refresh_interval": None  # default (1s)
        }, index)
