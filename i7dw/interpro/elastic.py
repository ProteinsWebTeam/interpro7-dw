#!/usr/bin/env python
# -*- coding: utf-8 -*-

import glob
import hashlib
import json
import logging
import os
import shutil
import time
from multiprocessing import Process, Queue
from tempfile import mkdtemp, mkstemp
from typing import Generator, Optional

from elasticsearch import Elasticsearch, helpers, exceptions

from . import mysql, supermatch
from .. import io, pdbe

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s: %(levelname)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)

LOADING_FILE = "loading"
EXTRA_INDEX = "others"


def init_dir(path):
    try:
        shutil.rmtree(path)
    except FileNotFoundError:
        pass
    finally:
        os.makedirs(path)
        open(os.path.join(path, LOADING_FILE), "w").close()


def parse_host(host: str) -> dict:
    host = host.split(':')
    if len(host) == 2:
        return {
            "host": host[0],
            "port": int(host[1])
        }
    else:
        return {
            "host": host,
            "port": 9200
        }


class DocumentProducer(Process):
    def __init__(self, ora_ippro: str, my_ippro: str, task_queue: Queue,
                 done_queue: Queue, outdir: str, min_overlap: int=20,
                 chunk_size: int=10000):
        super().__init__()
        self.ora_ippro = ora_ippro
        self.my_ippro = my_ippro
        self.task_queue = task_queue
        self.done_queue = done_queue
        self.outdir = mkdtemp(dir=outdir)
        self.min_overlap = min_overlap
        self.chunk_size = chunk_size

        self.dir_limit = 1000
        self.dir_count = 0

        self.entries = {}
        self.integrated = {}
        self.entry2set = {}
        self.proteomes = {}
        self.pfam = set()
        self.structures = {}
        self.protein2pdb = {}

    def run(self):
        # Get PDBe structures, entries, sets, and proteomes
        self.structures = pdbe.get_structures(self.ora_ippro)

        for pdb_id, s in self.structures.items():
            for acc in s["proteins"]:
                if acc in self.protein2pdb:
                    self.protein2pdb[acc].add(pdb_id)
                else:
                    self.protein2pdb[acc] = {pdb_id}

        self.entries = mysql.get_entries(self.my_ippro)
        self.integrated = {
            acc: e["integrated"]
            for acc, e in self.entries.items()
            if e["integrated"]
        }
        self.entry2set = {
            entry_ac: (set_ac, s["database"])
            for set_ac, s in mysql.get_sets(self.my_ippro).items()
            for entry_ac in s["members"]
        }
        self.proteomes = mysql.get_proteomes(self.my_ippro)

        # List Pfam entries (for IDA)
        self.pfam = {
            e["accession"]
            for e in self.entries.values()
            if e["database"] == "pfam"
        }

        documents = []
        cnt = 0
        types = {
            "protein": self.process_protein,
            "entry": self.process_entry,
            "taxonomy": self.process_taxonomy
        }

        for _type, chunk in iter(self.task_queue.get, None):
            fn = types[_type]

            for args in chunk:
                documents += fn(*args)

                if len(documents) >= self.chunk_size:
                    cnt += len(documents)
                    documents = self.dump(documents)
                    cnt -= len(documents)

        if documents:
            cnt += len(documents)
            self.dump(documents, strict_size=False)

        self.done_queue.put(cnt)

    def process_protein(self, accession: str, identifier: str, name: str,
                        database: str, length: int, comments: list,
                        matches: list, proteome_id: str, taxon: dict) -> list:
        # Prepare matches/supermatches
        entry_matches = {}
        supermatches = []
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
                "model_acc": m["model_ac"],
                "seq_feature": m["seq_feature"]
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
                entry = self.entries[entry_ac]
                pos_start = None
                pos_end = None
                for f in m["fragments"]:
                    if pos_start is None or f["start"] < pos_start:
                        pos_start = f["start"]
                    if pos_end is None or f["end"] > pos_end:
                        pos_end = f["end"]

                supermatches.append(
                    supermatch.Supermatch(
                        entry_ac,
                        entry["root"],
                        pos_start,
                        pos_end
                    )
                )

        if dom_arch:
            dom_arch = '-'.join(dom_arch)
            dom_arch_id = hashlib.sha1(dom_arch.encode("utf-8")).hexdigest()
        else:
            dom_arch = dom_arch_id = None

        # Merge overlapping supermatches
        sm_sets = supermatch.merge_supermatches(supermatches,
                                                self.min_overlap)
        for s in sm_sets:
            for sm in s.supermatches:
                for entry_ac in sm.get_entries():
                    # Add supermatches to Elastic matches
                    if entry_ac in entry_matches:
                        e = entry_matches[entry_ac]
                    else:
                        e = entry_matches[entry_ac] = []

                    e.append({
                        "fragments": [{"start": sm.start, "end": sm.end}],
                        "model_acc": None,
                        "seq_feature": None
                    })

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
            "protein_size": size,
            "protein_db": database,
            "text_protein": self._join(
                accession, identifier, name, database, comments
            ),

            "tax_id": taxon["taxId"],
            "tax_name": taxon["scientificName"],
            "tax_lineage": taxon["lineage"].strip().split(),
            "tax_rank": taxon["rank"],
            "text_taxonomy": self._join(
                taxon["taxId"], taxon["fullName"], taxon["rank"]
            ),
        })

        if proteome_id:
            p = self.proteomes[proteome_id]
            doc.update({
                "proteome_acc": proteome_id.lower(),
                "proteome_name": p["name"],
                "proteome_is_reference": p["is_reference"],
                "text_proteome": self._join(proteome_id, *list(p.values()))
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
                    "text_entry": self._join(
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
                    set_ac, set_db = _set
                    _doc.update({
                        "set_acc": set_ac.lower(),
                        "set_db": set_db,
                        # todo: implement set integration (e.g. pathways)
                        "set_integrated": [],
                        "text_set": self._join(set_ac, set_db)
                    })

                documents.append(_doc)
        else:
            documents.append(doc)

        # Add PDBe structures (and chains)
        _documents = []
        for pdbe_id in self.protein2pdb.get(accession, []):
            structure = self.structures[pdbe_id]
            text = self._join(
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

        if _documents:
            documents = _documents

        for doc in documents:
            doc["id"] = self._join(
                doc["protein_acc"], doc["proteome_acc"], doc["entry_acc"],
                doc["set_acc"], doc["structure_acc"],
                doc["structure_chain_acc"],
                separator='-'
            )

        return documents

    def process_entry(self, accession: str) -> list:
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
            "text_entry": self._join(
                entry["accession"], entry["name"],
                entry["type"], entry["descriptions"], *go_terms
            ),
            "entry_protein_locations": [],
            "entry_go_terms": go_terms
        })

        _set = self.entry2set.get(accession)
        if _set:
            set_ac, set_db = _set
            doc.update({
                "set_acc": set_ac.lower(),
                "set_db": set_db,
                # todo: implement set integration (e.g. pathways)
                "set_integrated": [],
                "text_set": self._join(set_ac, set_db)
            })
            doc["id"] = doc["entry_acc"] + '-' + doc["set_acc"]
        else:
            doc["id"] = doc["entry_acc"]

        return [doc]

    def process_taxonomy(self, taxon: dict) -> list:
        doc = self.init_document()
        doc.update({
            "tax_id": taxon["taxId"],
            "tax_name": taxon["scientificName"],
            "tax_lineage": taxon["lineage"].strip().split(),
            "tax_rank": taxon["rank"],
            "text_taxonomy": self._join(
                taxon["taxId"], taxon["fullName"], taxon["rank"]
            ),
            "id": taxon["taxId"]
        })
        return [doc]

    def dump(self, documents: list, strict_size: bool=True) -> list:
        for i in range(0, len(documents), self.chunk_size):
            chunk = documents[i:i+self.chunk_size]
            if len(chunk) == self.chunk_size or not strict_size:
                if self.dir_count + 1 == self.dir_limit:
                    # Too many files in directory: create a subdirectory
                    self.outdir = mkdtemp(dir=self.outdir)
                    self.dir_count = 0

                fd, path = mkstemp(dir=self.outdir)
                os.close(fd)

                with open(path, "wt") as fh:
                    json.dump(chunk, fh)

                os.rename(path, path + ".json")
                self.dir_count += 1
            else:
                return chunk

        return []

    @staticmethod
    def init_document() -> dict:
        return {
            "id": None,

            # Protein
            "protein_acc": None,
            "protein_length": None,
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

    @staticmethod
    def _join(*args, separator: str=' ') -> str:
        # Use underscore to NOT override Process.join()
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

        return separator.join(items).lower()


def _iter_proteins(store: io.Store, keys: list=list()) -> Generator:
    if keys:
        for k in sorted(keys):
            v = store.get(k)
            if v:
                yield k, v
    else:
        for k, v in store:
            yield k, v


def create_documents(ora_ippro: str, my_ippro: str, src_proteins: str,
                     src_names: str, src_comments: str, src_proteomes: str,
                     src_matches: str, outdir: str, **kwargs):
    processes = kwargs.get("processes", 1)
    chunk_size = kwargs.get("chunk_size", 100000)
    limit = kwargs.get("limit", 0)
    keys = kwargs.get("keys", [])

    processes = max(1, processes-1)  # minus one for parent process
    task_queue = Queue(processes)
    done_queue = Queue()
    workers = []
    for _ in range(processes):
        p = DocumentProducer(ora_ippro, my_ippro, task_queue, done_queue,
                             outdir)
        p.start()
        workers.append(p)

    # MySQL data
    logging.info("loading data from MySQL")
    taxa = mysql.get_taxa(my_ippro, lineage=True)
    integrated = {}
    entry_accessions = set()
    for entry_ac, e in mysql.get_entries(my_ippro).items():
        entry_accessions.add(entry_ac)
        if e["integrated"]:
            integrated[entry_ac] = e["integrated"]

    # Open stores
    proteins = io.Store(src_proteins)
    protein2names = io.Store(src_names)
    protein2comments = io.Store(src_comments)
    protein2proteome = io.Store(src_proteomes)
    protein2matches = io.Store(src_matches)

    tax_ids = set(taxa.keys())

    logging.info("starting")
    n_proteins = 0
    chunk = []
    entries_with_matches = set()
    enqueue_time = 0
    ts = time.time()

    for acc, protein in _iter_proteins(proteins, keys):
        tax_id = protein["taxon"]
        taxon = taxa[tax_id]

        name, other_names = protein2names.get(acc, (None, None))
        matches = protein2matches.get(acc, [])

        # Enqueue protein
        chunk.append((
            acc,
            protein["identifier"],
            name,
            "reviewed" if protein["isReviewed"] else "unreviewed",
            protein["length"],
            protein2comments.get(acc, []),
            matches,
            protein2proteome.get(acc),
            taxon
        ))

        if len(chunk) == chunk_size:
            t = time.time()
            task_queue.put(("protein", chunk))
            enqueue_time += time.time() - t
            chunk = []

        # Keep track of taxa associated to at least one protein
        try:
            tax_ids.remove(tax_id)
        except KeyError:
            pass

        # Keep track of entries with protein matches
        for m in matches:
            method_ac = m["method_ac"]
            entries_with_matches.add(method_ac)

            if method_ac in integrated:
                entries_with_matches.add(integrated[method_ac])

        n_proteins += 1
        if n_proteins == limit:
            break
        elif not n_proteins % 10000000:
            logging.info("{:>12,} ({:.0f} proteins/sec)".format(
                n_proteins, n_proteins / (time.time() - ts)
            ))

    if chunk:
        t = time.time()
        task_queue.put(("protein", chunk))
        enqueue_time += time.time() - t

    logging.info("{:>12,} ({:.0f} proteins/sec)".format(
        n_proteins, n_proteins / (time.time() - ts)
    ))

    if not keys:
        # Add entries without matches
        chunk = [
            (entry_ac,)
            for entry_ac in entry_accessions - entries_with_matches
        ]
        for i in range(0, len(chunk), chunk_size):
            t = time.time()
            task_queue.put(("entry", chunk[i:i+chunk_size]))
            enqueue_time += time.time() - t

        # Add taxa without proteins
        chunk = [(taxa[tax_id],) for tax_id in tax_ids]
        for i in range(0, len(chunk), chunk_size):
            t = time.time()
            task_queue.put(("taxonomy", chunk[i:i+chunk_size]))
            enqueue_time += time.time() - t

    logging.info("enqueue time: {:>10.0f} seconds".format(enqueue_time))

    # Poison pill
    for _ in workers:
        task_queue.put(None)

    # Closing stores
    proteins.close()
    protein2names.close()
    protein2comments.close()
    protein2proteome.close()
    protein2matches.close()

    n_docs = sum([done_queue.get() for _ in workers])

    # Wait for workers to finish
    for p in workers:
        p.join()

    # Delete loading file so Loaders know that all files are generated
    os.remove(os.path.join(outdir, LOADING_FILE))

    logging.info("complete: {:,} documents".format(n_docs))


class DocumentLoader(Process):
    def __init__(self, host: dict, doc_type: str, queue_in: Queue,
                 queue_out: Queue, **kwargs):
        super().__init__()
        self.host = host
        self.type = doc_type
        self.queue_in = queue_in
        self.queue_out = queue_out
        self.max_bytes = kwargs.get("max_bytes", 100 * 1024 * 1024)
        self.suffix = kwargs.get("suffix", "")
        self.threads = kwargs.get("threads", 4)
        self.err_log = kwargs.get("err_log")

    def run(self):
        es = Elasticsearch([self.host])

        # Disable Elastic logger
        tracer = logging.getLogger("elasticsearch")
        tracer.setLevel(logging.CRITICAL + 1)

        if self.err_log:
            err_dst = self.err_log
        else:
            err_dst = os.devnull

        with open(err_dst, "wt") as err:
            for filepath in iter(self.queue_in.get, None):
                with open(filepath, "rt") as fh:
                    documents = json.load(fh)

                actions = []
                for doc in documents:
                    # Define which index to use
                    if doc["entry_db"]:
                        index = doc["entry_db"] + self.suffix
                    else:
                        index = EXTRA_INDEX + self.suffix

                    actions.append({
                        "_op_type": "index",
                        "_index": index,
                        "_type": self.type,
                        "_id": doc["id"],
                        "_source": doc
                    })

                gen = helpers.parallel_bulk(
                    es, actions,
                    thread_count=self.threads,
                    queue_size=self.threads,
                    # disable chunk_size (num of docs)
                    # to only rely on max_chunk_bytes (bytes)
                    chunk_size=-1,
                    max_chunk_bytes=self.max_bytes,
                    raise_on_exception=False,
                    raise_on_error=False
                )

                success = True
                for status, item in gen:
                    if not status:
                        success = False
                        try:
                            err.write("{}\n".format(item))
                        except:
                            pass

                self.queue_out.put((filepath, success))


def listen_files(src: str, seconds: int=60):
    pathname = os.path.join(src, "**", "*.json")
    files = set()
    stop = False
    while True:
        for filepath in glob.iglob(pathname, recursive=True):
            if filepath not in files:
                files.add(filepath)
                yield filepath

        if stop:
            break
        elif not os.path.isfile(os.path.join(src, LOADING_FILE)):
            # All files ready, but loop one last time
            stop = True
        else:
            time.sleep(seconds)


def create_indices(body_f: str, shards_f: str, my_ippro: str, host: dict,
                   shards: int=5, suffix: str=""):
    # Load default settings and property mapping
    with open(body_f, "rt") as fh:
        body = json.load(fh)

    if shards_f:
        # Custom number of shards
        with open(shards_f, "rt") as fh:
            custom_shards = json.load(fh)

        # Force keys to be in lower case
        custom_shards = {k.lower(): v for k, v in custom_shards.items()}
    else:
        custom_shards = {}

    # Load databases (base of indices)
    databases = list(mysql.get_entry_databases(my_ippro).keys())

    # Establish connection
    es = Elasticsearch([host])

    # Disable Elastic logger
    tracer = logging.getLogger("elasticsearch")
    tracer.setLevel(logging.CRITICAL + 1)

    # Create indices
    for index in databases + [EXTRA_INDEX]:
        try:
            n_shards = custom_shards[index]
        except KeyError:
            # No custom number of shards for this index
            n_shards = shards

        """
        Change settings for large bulk imports:

        https://www.elastic.co/guide/en/elasticsearch/reference/current/indices-update-settings.html
        https://www.elastic.co/guide/en/elasticsearch/guide/current/indexing-performance.html
        """
        try:
            body["settings"].update({
                "number_of_shards": n_shards,
                "number_of_replicas": 0,    # default: 1
                "refresh_interval": -1      # default: 1s
            })
        except KeyError:
            body["settings"] = {
                "number_of_shards": n_shards,
                "number_of_replicas": 0,    # default: 1
                "refresh_interval": -1      # default: 1s
            }

        index += suffix

        # Make sure the index is deleted
        while True:
            try:
                res = es.indices.delete(index)
            except exceptions.NotFoundError:
                break
            except Exception as e:
                logging.error("{}: {}".format(type(e), e))
                time.sleep(10)
            else:
                break

        # And make sure it's created
        while True:
            try:
                es.indices.create(index, body=body)
            except exceptions.RequestError as e:
                break  # raised if index exists
            except Exception as e:
                logging.error("{}: {}".format(type(e), e))
                time.sleep(10)
            else:
                break


def index_documents(my_ippro: str, host: str, doc_type: str, src: str,
                    **kwargs) -> list:
    alias = kwargs.get("alias")
    body_f = kwargs.get("body")
    _create_indices = kwargs.get("create_indices", False)
    shards_f = kwargs.get("custom_shards")
    default_shards = kwargs.get("default_shards", 5)
    err_prefix = kwargs.get("err_prefix")
    files = kwargs.get("files", [])
    max_retries = kwargs.get("max_retries", 0)
    processes = kwargs.get("processes", 1)
    raise_on_error = kwargs.get("raise_on_error", True)
    suffix = kwargs.get("suffix", "").lower()

    # Parse Elastic host (str -> dict)
    _host = parse_host(host)

    logging.info("indexing documents to {}".format(_host["host"]))

    if _create_indices and body_f:
        create_indices(body_f, shards_f, my_ippro, _host, default_shards,
                       suffix)

    processes = max(1, processes - 1)  # consider parent process

    n_retries = 0
    while True:
        queue_in = Queue()
        queue_out = Queue()
        workers = []
        for i in range(processes):
            if err_prefix:
                err_log = err_prefix + str(i+1) + ".err"
            else:
                err_log = None

            w = DocumentLoader(_host, doc_type, queue_in, queue_out,
                               suffix=suffix, err_log=err_log)
            w.start()
            workers.append(w)

        n_files = 0
        if files:
            for filepath in files:
                queue_in.put(filepath)
                n_files += 1
        else:
            for filepath in listen_files(src):
                queue_in.put(filepath)
                n_files += 1

        for _ in workers:
            queue_in.put(None)

        # Get files and their status from child processes
        files = []
        i = 0
        for i in range(n_files):
            filepath, status = queue_out.get()
            if not status:
                files.append(filepath)

            if not (i + 1) % 1000:
                logging.info("files: {:>10,} / {:,} "
                             "({:,} failed)".format(i+1, n_files, len(files)))
        logging.info("files: {:>10,} / {:,} "
                     "({:,} failed)".format(i+1, n_files, len(files)))

        # Wait for workers to terminate
        for w in workers:
            w.join()

        if not files or n_retries == max_retries:
            break

        n_retries += 1

    if files:
        if raise_on_error:
            for filepath in files:
                logging.critical("could not load {}".format(filepath))

            raise RuntimeError("{:,} files not loaded".format(len(files)))
        else:
            return files
    elif alias:
        """
        Do NOT delete old indices as they are used in production
        They will be deleted when we switch transparently between them
            and the new indices
        """
        update_alias(my_ippro, host, alias=alias, suffix=suffix, delete=False)

    logging.info("complete")
    return files


def update_alias(my_ippro: str, host: str, alias: str, **kwargs):
    suffix = kwargs.get("suffix", "").lower()
    delete = kwargs.get("delete", False)

    databases = list(mysql.get_entry_databases(my_ippro).keys())
    new_indices = set()
    for index in databases + [EXTRA_INDEX]:
        new_indices.add(index + suffix)

    es = Elasticsearch([parse_host(host)])

    # Disable Elastic logger
    tracer = logging.getLogger("elasticsearch")
    tracer.setLevel(logging.CRITICAL + 1)

    exists = es.indices.exists_alias(name=alias)
    if exists:
        # Alias already exists: update it

        # Indices currently using the alias
        indices = set(es.indices.get_alias(name=alias))

        actions = []
        for index in new_indices:
            try:
                # If passes: new index is already using the alias
                indices.remove(index)
            except KeyError:
                # Otherwise, add the alias to the new index
                actions.append({
                    "add": {
                        "index": index,
                        "alias": alias
                    }
                })

        # Remove the alias from the old indices
        for index in indices:
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

        if delete:
            # Delete old indices that used the alias
            for index in indices:
                while True:
                    try:
                        res = es.indices.delete(index)
                    except exceptions.NotFoundError:
                        break
                    except Exception as e:
                        logging.error("{}: {}".format(type(e), e))
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
