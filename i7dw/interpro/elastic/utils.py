# -*- coding: utf-8 -*-

import glob
import hashlib
import json
import logging
import os
import shutil
import time
from multiprocessing import Process, Queue
from typing import List, Optional

from elasticsearch import Elasticsearch, exceptions, helpers

from i7dw import io, logger, pdbe
from i7dw.interpro import condense, mysql


LOADING_FILE = "loading"
NODB_INDEX = "others"


class DocumentProducer(Process):
    def __init__(self, ora_ipr: str, my_ipr: str, task_queue: Queue,
                 done_queue: Queue, outdir: str, min_overlap: int=20):
        super().__init__()
        self.ora_ipr = ora_ipr
        self.my_ipr = my_ipr
        self.task_queue = task_queue
        self.done_queue = done_queue
        self.min_overlap = min_overlap
        self.organizer = io.JsonFileOrganizer(outdir)

        self.entries = {}
        self.integrated = {}
        self.entry2set = {}
        self.proteomes = {}
        self.pfam = set()
        self.structures = {}
        self.protein2pdb = {}

    def run(self):
        # Get PDBe structures, entries, sets, and proteomes
        self.structures = pdbe.get_structures(self.ora_ipr)

        for pdb_id, s in self.structures.items():
            for acc in s["proteins"]:
                if acc in self.protein2pdb:
                    self.protein2pdb[acc].add(pdb_id)
                else:
                    self.protein2pdb[acc] = {pdb_id}

        self.entries = mysql.entry.get_entries(self.my_ipr)
        self.integrated = {
            acc: e["integrated"]
            for acc, e in self.entries.items()
            if e["integrated"]
        }
        self.entry2set = {
            entry_ac: (set_ac, s["database"], s["name"], s["description"])
            for set_ac, s in mysql.entry.get_sets(self.my_ipr).items()
            for entry_ac in s["members"]
        }
        self.proteomes = mysql.proteome.get_proteomes(self.my_ipr)

        # List Pfam entries (for IDA)
        self.pfam = {
            e["accession"]
            for e in self.entries.values()
            if e["database"] == "pfam"
        }

        cnt = 0
        types = {
            "protein": self.process_protein,
            "entry": self.process_entry,
            "taxonomy": self.process_taxonomy
        }

        for _type, chunk in iter(self.task_queue.get, None):
            fn = types[_type]

            for args in chunk:
                for doc in fn(*args):
                    self.organizer.add(doc)
                    cnt += 1

        self.organizer.flush()
        self.done_queue.put(cnt)

    def process_protein(self, accession: str, identifier: str, name: str,
                        database: str, is_fragment: bool, length: int,
                        comments: list, matches: list, proteome_id: str,
                        taxon: dict) -> List[dict]:
        entry_matches = {}
        to_condense = {}
        dom_arch = []
        dom_entries = set()
        for m in matches:
            method_ac = m["method_ac"]
            if method_ac in entry_matches:
                e = entry_matches[method_ac]
            else:
                e = entry_matches[method_ac] = []

            e.append({
                "fragments": m["fragments"],
                "model_acc": m["model_ac"]
            })

            entry_ac = self.integrated.get(method_ac)
            if method_ac in self.pfam:
                dom_entries.add(method_ac)
                if entry_ac:
                    dom_entries.add(entry_ac)
                    dom_arch.append("{}:{}".format(method_ac, entry_ac))
                else:
                    dom_arch.append("{}".format(method_ac))

            if entry_ac:
                if entry_ac in to_condense:
                    to_condense[entry_ac].append(m["fragments"])
                else:
                    to_condense[entry_ac] = [m["fragments"]]

        for entry_ac, locations in condense(to_condense).items():
            entry_matches[entry_ac] = locations

        if dom_arch:
            dom_arch = '-'.join(dom_arch)
            dom_arch_id = hashlib.sha1(dom_arch.encode("utf-8")).hexdigest()
        else:
            dom_arch = dom_arch_id = None

        if length <= 100:
            size = "small"
        elif length <= 1000:
            size = "medium"
        else:
            size = "large"

        doc = self.init_document()

        # Add protein info
        doc.update({
            "protein_acc": accession.lower(),
            "protein_length": length,
            "protein_is_fragment": is_fragment,
            "protein_size": size,
            "protein_db": database,
            "text_protein": joinitems(
                accession, identifier, name, database, comments
            ),

            "tax_id": taxon["taxId"],
            "tax_name": taxon["scientificName"],
            "tax_lineage": taxon["lineage"],
            "tax_rank": taxon["rank"],
            "text_taxonomy": joinitems(
                taxon["taxId"], taxon["fullName"], taxon["rank"]
            ),
        })

        if proteome_id:
            p = self.proteomes[proteome_id]
            doc.update({
                "proteome_acc": proteome_id.lower(),
                "proteome_name": p["name"],
                "proteome_is_reference": p["is_reference"],
                "text_proteome": joinitems(proteome_id, *list(p.values()))
            })

        documents = []
        if entry_matches:
            # Add entries
            for entry_ac in entry_matches:
                entry = self.entries[entry_ac]
                if entry["integrated"]:
                    entry["integrated"] = entry["integrated"].lower()

                go_terms = [t["identifier"] for t in entry["go_terms"]]
                _doc = doc.copy()
                _doc.update({
                    "entry_acc": entry["accession"].lower(),
                    "entry_db": entry["database"],
                    "entry_type": entry["type"],
                    "entry_date": entry["date"].strftime("%Y-%m-%d"),
                    "entry_integrated": entry["integrated"],
                    "text_entry": joinitems(
                        entry["accession"], entry["name"],
                        entry["type"], entry["descriptions"], *go_terms
                    ),
                    "entry_protein_locations": entry_matches[entry_ac],
                    "entry_go_terms": go_terms
                })

                if entry["accession"] in dom_entries:
                    _doc.update({
                        "ida_id": dom_arch_id,
                        "ida": dom_arch
                    })

                _set = self.entry2set.get(entry_ac)
                if _set:
                    set_ac, set_db, set_name, set_descr = _set
                    _doc.update({
                        "set_acc": set_ac.lower(),
                        "set_db": set_db,
                        # todo: implement set integration (e.g. pathways)
                        "set_integrated": [],
                        "text_set": joinitems(set_ac, set_db, set_name,
                                              set_descr)
                    })

                documents.append(_doc)
        else:
            documents.append(doc)

        # Add PDBe structures (and chains)
        _documents = []
        for pdbe_id in self.protein2pdb.get(accession, []):
            structure = self.structures[pdbe_id]
            text = joinitems(
                pdbe_id,
                structure["evidence"],
                structure["name"],
                ' '.join([
                    pub["title"]
                    for pub in structure["citations"].values()
                    if pub.get("title")
                ])
            )

            for doc in documents:
                chains = structure["proteins"][accession]

                _doc = doc.copy()
                _doc.update({
                    "structure_acc": pdbe_id.lower(),
                    "structure_resolution": structure["resolution"],
                    "structure_date": structure["date"].strftime("%Y-%m-%d"),
                    "structure_evidence": structure["evidence"],
                    "protein_structure": chains
                })

                for chain_id, chain in chains.items():
                    _doc_chain = _doc.copy()
                    _doc_chain.update({
                        "structure_chain_acc": chain_id,
                        "structure_protein_locations": [
                            {
                                "fragments": [
                                    {
                                        "start": fragment["protein_start"],
                                        "end": fragment["protein_end"]
                                    }
                                ]
                            } for fragment in chain
                        ],
                        "structure_chain": "{} - {}".format(
                            pdbe_id, chain_id
                        ),
                        "text_structure": "{} {}".format(chain_id, text)
                    })
                    _documents.append(_doc_chain)

        return _documents if _documents else documents

    def process_entry(self, accession: str) -> List[dict]:
        entry = self.entries[accession]
        if entry["integrated"]:
            entry["integrated"] = entry["integrated"].lower()

        go_terms = [t["identifier"] for t in entry["go_terms"]]
        doc = self.init_document()
        doc.update({
            "entry_acc": entry["accession"].lower(),
            "entry_db": entry["database"],
            "entry_type": entry["type"],
            "entry_date": entry["date"].strftime("%Y-%m-%d"),
            "entry_integrated": entry["integrated"],
            "text_entry": joinitems(
                entry["accession"], entry["name"],
                entry["type"], entry["descriptions"], *go_terms
            ),
            "entry_protein_locations": [],
            "entry_go_terms": go_terms
        })

        _set = self.entry2set.get(accession)
        if _set:
            set_ac, set_db, set_name, set_descr = _set
            doc.update({
                "set_acc": set_ac.lower(),
                "set_db": set_db,
                # todo: implement set integration (e.g. pathways)
                "set_integrated": [],
                "text_set": joinitems(set_ac, set_db, set_name, set_descr)
            })

        return [doc]

    def process_taxonomy(self, taxon: dict) -> List[dict]:
        doc = self.init_document()
        doc.update({
            "tax_id": taxon["taxId"],
            "tax_name": taxon["scientificName"],
            "tax_lineage": taxon["lineage"],
            "tax_rank": taxon["rank"],
            "text_taxonomy": joinitems(
                taxon["taxId"], taxon["fullName"], taxon["rank"]
            )
        })
        return [doc]

    @staticmethod
    def init_document() -> dict:
        return {
            # Protein
            "protein_acc": None,
            "protein_length": None,
            "protein_is_fragment": None,
            "protein_size": None,
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
            "protein_structure_chain": None,
            "structure_chain": None,
            "text_structure": None,

            # Domain architecture
            "ida_id": None,
            "ida": None
        }


def is_ready(path: str) -> bool:
    return not os.path.isfile(os.path.join(path, LOADING_FILE))


def set_ready(path: str):
    os.remove(os.path.join(path, LOADING_FILE))


def joinitems(*args, separator: str = ' ') -> str:
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


def create_indices(hosts: List[str], indices: List[str], body_path: str, **kwargs):
    delete_all = kwargs.get("delete_all", False)
    num_shards = kwargs.get("num_shards", 5)
    shards_path = kwargs.get("shards_path")
    suffix = kwargs.get("suffix", "")

    # Load default settings and property mapping
    with open(body_path, "rt") as fh:
        body = json.load(fh)

    # Custom number of shards
    if shards_path:
        with open(shards_path, "rt") as fh:
            shards = json.load(fh)

        # Force keys to be in lower case
        shards = {k.lower(): v for k, v in shards.items()}
    else:
        shards = {}

    # Establish connection
    es = Elasticsearch(hosts, timeout=30)

    # Disable Elastic logger
    tracer = logging.getLogger("elasticsearch")
    tracer.setLevel(logging.CRITICAL + 1)

    if delete_all:
        # Delete all existing indices
        es.indices.delete('*', allow_no_indices=True)

    # Create indices
    for index in indices:
        try:
            n_shards = shards[index]
        except KeyError:
            # No custom number of shards for this index
            n_shards = num_shards

        """
        Change settings for large bulk imports:

        https://www.elastic.co/guide/en/elasticsearch/reference/current/indices-update-settings.html
        https://www.elastic.co/guide/en/elasticsearch/guide/current/indexing-performance.html
        """
        body["settings"].update({
            "index": {
                # Static settings
                "number_of_shards": n_shards,
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
                res = es.indices.delete(index)
            except exceptions.NotFoundError:
                break
            except Exception as exc:
                logger.error("{}: {}".format(type(exc), exc))
                time.sleep(10)
            else:
                break

        # And make sure it's created
        while True:
            try:
                res = es.indices.create(index, body=body)
            except exceptions.RequestError as exc:
                raise exc
            except Exception as exc:
                logger.error("{}: {}".format(type(exc), exc))
                time.sleep(10)
            else:
                break


def iter_json_files(src: str, seconds: int=60):
    pathname = os.path.join(src, "**", "*.json")
    files = set()
    active = True

    while True:
        for path in glob.iglob(pathname, recursive=True):
            if path not in files:
                files.add(path)
                yield path

        if not active:
            break
        if is_ready(src):
            # All files ready, but loop one last time
            active = False
        else:
            time.sleep(seconds)


def save_failed_docs(task_queue: Queue, done_queue: Queue,
                     dst: Optional[str]=None, write_back: bool=False):
    if dst:
        """
        Ensure the directory does not exist
        as we don't want files from a previous run to be considered
        """
        try:
            shutil.rmtree(dst)
        except FileNotFoundError:
            pass
        finally:
            os.makedirs(dst)
            organizer = io.JsonFileOrganizer(dst)
    else:
        organizer = None

    while True:
        total_successful = 0
        total_failed = 0
        num_files = 0
        it = iter(task_queue.get, None)

        for path, num_successful, failed, errors in it:
            total_successful += num_successful
            total_failed += len(failed)

            if failed:
                if organizer:
                    for doc in failed:
                        organizer.add(doc)
                    organizer.flush()
                elif write_back:
                    with open(path, "wt") as fh:
                        json.dump(failed, fh)

                for item in errors:
                    logger.debug(item)
            elif write_back:
                os.remove(path)

            num_files += 1
            if not num_files % 1000:
                logger.info("documents indexed: {:>15,} "
                            "({:,} failed)".format(total_successful,
                                                   total_failed))

        logger.info("documents indexed: {:>15,} "
                    "({:,} failed)".format(total_successful, total_failed))

        done_queue.put(total_failed)

        """
        Wait for instruction from parent:
            - True: we continue
            - False: max retries reached: quit
        """
        if not task_queue.get():
            break


class DocumentController(object):
    def __init__(self, suffix: str):
        self.suffix = suffix

    def wrap(self, doc: dict) -> dict:
        if doc["entry_db"]:
            idx = doc["entry_db"] + self.suffix
        else:
            idx = NODB_INDEX + self.suffix

        if doc["protein_acc"]:
            args = (doc["protein_acc"], doc["proteome_acc"], doc["entry_acc"],
                    doc["set_acc"], doc["structure_acc"],
                    doc["structure_chain_acc"])
        elif doc["entry_acc"]:
            args = (doc["entry_acc"], doc["set_acc"])
        else:
            args = (doc["tax_id"],)

        return {
            "_op_type": "index",
            "_index": idx,
            "_id": joinitems(*args, separator='-'),
            "_source": doc
        }


class DocumentLoader(Process):
    def __init__(self, hosts: List, suffix: str, task_queue: Queue,
                 done_queue: Queue, **kwargs):
        super().__init__()
        self.hosts = hosts
        self.controller = DocumentController(suffix)
        self.task_queue = task_queue
        self.done_queue = done_queue

        # elasticsearch-py defaults
        self.chunk_size = kwargs.get("chunk_size", 500)
        self.max_bytes = kwargs.get("max_bytes", 100 * 1024 * 1024)
        self.threads = kwargs.get("threads", 4)
        self.timeout = kwargs.get("timeout", 10)

    def run(self):
        es = Elasticsearch(self.hosts, timeout=self.timeout)

        # Disable Elastic logger
        tracer = logging.getLogger("elasticsearch")
        tracer.setLevel(logging.CRITICAL + 1)

        for filepath in iter(self.task_queue.get, None):
            with open(filepath, "rt") as fh:
                documents = json.load(fh)

            bulk = helpers.parallel_bulk(
                es, map(self.controller.wrap, documents),
                thread_count=self.threads,
                queue_size=self.threads,
                chunk_size=self.chunk_size,
                max_chunk_bytes=self.max_bytes,
                raise_on_exception=False,
                raise_on_error=False
            )

            num_successful = 0
            failed = []
            errors = []
            for i, (status, item) in enumerate(bulk):
                if status:
                    num_successful += 1
                else:
                    failed.append(documents[i])
                    try:
                        del item["index"]["data"]
                    except KeyError:
                        pass
                    finally:
                        errors.append(item)

            self.done_queue.put((filepath, num_successful, failed, errors))


def index_documents(hosts: List[str], src: str, **kwargs) -> bool:
    suffix = kwargs.get("suffix", "")
    dst = kwargs.get("dst")
    max_retries = kwargs.get("max_retries", 0)
    processes = kwargs.get("processes", 1)
    raise_on_error = kwargs.get("raise_on_error", True)
    write_back = kwargs.get("write_back", False)

    file_queue = Queue()
    fail_queue = Queue(maxsize=processes)
    count_queue = Queue()
    organizer = Process(target=save_failed_docs,
                        args=(fail_queue, count_queue, dst, write_back))
    organizer.start()

    processes = max(1, processes-2)  # parent process + organizer
    num_retries = 0
    while True:
        logger.info("indexing documents (try #{})".format(num_retries + 1))
        workers = []

        for i in range(processes):
            w = DocumentLoader(hosts, suffix, file_queue, fail_queue, **kwargs)
            w.start()
            workers.append(w)

        for path in iter_json_files(src):
            file_queue.put(path)

        # All files enqueued
        for _ in workers:
            file_queue.put(None)

        # Wait for workers to complete
        for w in workers:
            w.join()

        # Inform organizer that we want the count of failed documents
        fail_queue.put(None)
        num_failed = count_queue.get()

        if num_failed and num_retries < max_retries:
            num_retries += 1
            fail_queue.put(True)   # continue
        else:
            fail_queue.put(False)  # stop
            organizer.join()
            break

    if num_failed:
        if raise_on_error:
            raise RuntimeError("{:,} documents not indexed".format(num_failed))
        else:
            logger.error("{:,} documents not indexed".format(num_failed))
            return False
    else:
        logger.info("complete")
        return True


def update_alias(hosts: List[str], indices: List[str], alias: str, **kwargs):
    delete_removed = kwargs.get("delete_removed", False)
    suffix = kwargs.get("suffix", "")

    new_indices = {index + suffix for index in indices}
    es = Elasticsearch(hosts, timeout=30)

    # Disable Elastic logger
    tracer = logging.getLogger("elasticsearch")
    tracer.setLevel(logging.CRITICAL + 1)

    exists = es.indices.exists_alias(name=alias)
    if exists:
        # Alias already exists: update it

        # Indices currently using the alias
        cur_indices = set(es.indices.get_alias(name=alias))

        actions = []
        for index in new_indices:
            try:
                # If passes: new index is already using the alias
                cur_indices.remove(index)
            except KeyError:
                # Otherwise, add the alias to the new index
                actions.append({
                    "add": {
                        "index": index,
                        "alias": alias
                    }
                })

        # Remove the alias from the current indices
        for index in cur_indices:
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

        if delete_removed:
            # Delete old indices that used the alias
            for index in cur_indices:
                while True:
                    try:
                        res = es.indices.delete(index)
                    except exceptions.NotFoundError:
                        break
                    except Exception as exc:
                        logger.error("{}: {}".format(type(exc), exc))
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
