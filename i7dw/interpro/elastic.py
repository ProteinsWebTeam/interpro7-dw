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

from elasticsearch import Elasticsearch, helpers, exceptions

from . import mysql, supermatch
from .. import io, pdbe

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s: %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
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
    def __init__(self, ora_ippro: str, my_ippro: str, queue_in: Queue,
                 outdir: str, **kwargs):
        super().__init__()
        self.ora_ippro = ora_ippro
        self.my_ippro = my_ippro
        self.queue_in = queue_in
        self.outdir = mkdtemp(dir=outdir)
        self.min_overlap = kwargs.get("min_overlap", 20)
        self.chunk_size = kwargs.get("chunk_size", 10000)

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
        logging.info("{} ({}) started".format(self.name, os.getpid()))

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

        while True:
            task = self.queue_in.get()
            if task is None:
                break

            _type, chunk = task
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

        logging.info("{} ({}) terminated ({:,} documents)".format(
            self.name, os.getpid(), cnt)
        )

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
                        "fragments": [{'start': sm.start, 'end': sm.end}],
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

        # Accession in lower case
        accession_lc = accession.lower()

        # Add protein info
        doc.update({
            "protein_acc": accession_lc,
            "protein_length": length,
            "protein_size": size,
            "protein_db": database,
            "text_protein": self._join(
                accession_lc, identifier, name, database, comments
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
                "proteome_acc": proteome_id,
                "proteome_name": p["name"],
                "proteome_is_reference": p["is_reference"],
                "text_proteome": self._join(proteome_id, *list(p.values()))
            })

        documents = []
        if entry_matches:
            # Add entries
            for entry_ac in entry_matches:
                entry = self.entries[entry_ac]
                go_terms = [t["identifier"] for t in entry["go_terms"]]
                _doc = doc.copy()
                _doc.update({
                    "entry_acc": entry["accession"],
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
                        "set_acc": set_ac,
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
                _doc = doc.copy()
                _doc.update({
                    "structure_acc": pdbe_id,
                    "structure_resolution": structure["resolution"],
                    "structure_date": structure["date"].strftime("%Y-%m-%d"),
                    "structure_evidence": structure["evidence"]
                })

                for chain_id in structure["proteins"][accession]:
                    chain = structure["proteins"][accession][chain_id]
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
        go_terms = [t["identifier"] for t in entry["go_terms"]]
        doc = self.init_document()
        doc.update({
            "entry_acc": entry["accession"],
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
                "set_acc": set_ac,
                "set_db": set_db,
                # todo: implement set integration (e.g. pathways)
                "set_integrated": [],
                "text_set": self._join(set_ac, set_db)
            })
            doc["id"] = "{}-{}".format(entry["accession"], set_ac)
        else:
            doc["id"] = entry["accession"]

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


def create_documents(ora_ippro: str, my_ippro: str, src_proteins: str,
                     src_names: str, src_comments: str, src_proteomes: str,
                     src_matches: str, outdir: str, **kwargs):
    processes = kwargs.get("processes", 1)
    chunk_size = kwargs.get("chunk_size", 100000)
    limit = kwargs.get("limit", 0)

    doc_queue = Queue(processes)

    workers = [
        DocumentProducer(ora_ippro, my_ippro, doc_queue, outdir)
        for _ in range(max(1, processes-1))
    ]

    for p in workers:
        p.start()

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
    for acc, protein in proteins:
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
            doc_queue.put(("protein", chunk))
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
        elif not n_proteins % 1000000:
            logging.info("{:>12,} ({:.0f} proteins/sec)".format(
                n_proteins, n_proteins / (time.time() - ts)
            ))

    if chunk:
        t = time.time()
        doc_queue.put(("protein", chunk))
        enqueue_time += time.time() - t

    logging.info("{:>12,} ({:.0f} proteins/sec)".format(
        n_proteins, n_proteins / (time.time() - ts)
    ))

    # Add entries without matches
    chunk = [
        (entry_ac,)
        for entry_ac in entry_accessions - entries_with_matches
    ]
    for i in range(0, len(chunk), chunk_size):
        t = time.time()
        doc_queue.put(("entry", chunk[i:i+chunk_size]))
        enqueue_time += time.time() - t

    # Add taxa without proteins
    chunk = [(taxa[tax_id],) for tax_id in tax_ids]
    for i in range(0, len(chunk), chunk_size):
        t = time.time()
        doc_queue.put(("taxonomy", chunk[i:i+chunk_size]))
        enqueue_time += time.time() - t

    logging.info("enqueue time: {:>10.0f} seconds".format(enqueue_time))

    # Poison pill
    for _ in workers:
        doc_queue.put(None)

    # Closing stores
    proteins.close()
    protein2names.close()
    protein2comments.close()
    protein2proteome.close()
    protein2matches.close()

    # Wait for workers to finish
    for p in workers:
        p.join()

    # Delete loading file so Loaders know that all files are generated
    os.remove(os.path.join(outdir, LOADING_FILE))

    logging.info("complete")


class DocumentLoader(Process):
    def __init__(self, host, doc_type, queue_in, queue_out, **kwargs):
        super().__init__()
        self.host = host
        self.type = doc_type
        self.queue_in = queue_in
        self.queue_out = queue_out
        self.max_bytes = kwargs.get("max_bytes", 100 * 1024 * 1024)
        self.suffix = kwargs.get("suffix", "")
        self.threads = kwargs.get("threads", 4)

    def run(self):
        logging.info("{} ({}) started".format(self.name, os.getpid()))
        es = Elasticsearch([self.host])

        # Disable Elastic logger
        tracer = logging.getLogger("elasticsearch")
        tracer.setLevel(logging.CRITICAL + 1)

        failed_files = []
        total_documents = 0
        total_errors = 0
        while True:
            filepath = self.queue_in.get()
            if filepath is None:
                break

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
                    '_op_type': 'index',
                    '_index': index,
                    '_type': self.type,
                    '_id': doc["id"],
                    '_source': doc
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

            total_documents += len(actions)
            n_errors = len([item for status, item in gen if not status])
            total_errors += n_errors

            if n_errors:
                failed_files.append(filepath)

            # TODO: remove after debug
            logging.info(
                "{} ({}) loaded {}: {}".format(
                    self.name, os.getpid(), filepath,
                    "error" if n_errors else "success"
                )
            )

        self.queue_out.put(failed_files)
        logging.info(
            "{} ({}) terminated "
            "({} documents, {} errors)".format(
                self.name, os.getpid(), total_documents, total_errors
            )
        )


def index_documents(my_ippro: str, host: str, doc_type: str,
                    properties: str, src: str, **kwargs):
    indices = kwargs.get("indices")
    create_indices = kwargs.get("create_indices", True)
    processes = kwargs.get("processes", 1)
    shards = kwargs.get("shards", 5)
    suffix = kwargs.get("suffix", "").lower()
    limit = kwargs.get("limit", 0)
    files = kwargs.get("files", [])

    # Parse Elastic host (str -> dict)
    _host = parse_host(host)

    logging.info("indexing documents to {}".format(_host["host"]))

    if create_indices:
        # Load property mapping
        with open(properties, "rt") as fh:
            properties = json.load(fh)

        if indices:
            # Custom number of shards
            with open(indices, "rt") as fh:
                indices = json.load(fh)

            # Force keys to be in lower case
            indices = {k.lower(): v for k, v in indices.items()}
        else:
            indices = {}

        # Load databases (base of indices)
        databases = list(mysql.get_entry_databases(my_ippro).keys())

        # Establish connection
        es = Elasticsearch([_host])

        # Disable Elastic logger
        tracer = logging.getLogger("elasticsearch")
        tracer.setLevel(logging.CRITICAL + 1)

        # Create indices
        for index in databases + [EXTRA_INDEX]:
            try:
                n_shards = indices[index]
            except KeyError:
                n_shards = shards

            index += suffix

            # https://www.elastic.co/guide/en/elasticsearch/reference/current/indices-update-settings.html
            # https://www.elastic.co/guide/en/elasticsearch/guide/current/indexing-performance.html
            while True:
                try:
                    res = es.indices.delete(index)
                except exceptions.ConnectionTimeout:
                    pass
                except exceptions.NotFoundError:
                    break
                else:
                    break

            body = {
                "mappings": {
                    doc_type: {
                        "properties": properties
                    }
                },
                "settings": {
                    "number_of_shards": n_shards,
                    "number_of_replicas": 0,    # default: 1
                    "refresh_interval": -1      # default: 1s
                }
            }

            while True:
                try:
                    es.indices.create(index, body=body)
                except exceptions.ConnectionTimeout:
                    pass
                except exceptions.RequestError:
                    break
                else:
                    break

    queue_in = Queue()
    queue_out = Queue()
    workers = [
        DocumentLoader(_host, doc_type, queue_in, queue_out, suffix=suffix)
        for _ in range(max(1, processes-1))
    ]
    for l in workers:
        l.start()

    if files:
        for filepath in files:
            queue_in.put(filepath)
    else:
        pathname = os.path.join(src, "**", "*.json")
        files = set()
        stop = False
        while True:
            for filepath in glob.iglob(pathname, recursive=True):
                if filepath not in files:
                    files.add(filepath)
                    queue_in.put(filepath)

                    if len(files) == limit:
                        stop = True
                        break

            if stop:
                logging.info("{:,} files to load".format(len(files)))
                break
            elif not os.path.isfile(os.path.join(src, LOADING_FILE)):
                # All files ready, but loop one last time
                stop = True
            else:
                time.sleep(60)

    """
    Get files that failed to load BEFORE joining child processesses
    to avoid deadlocks.

    From the `multiprocessing` docs:
        if a child process has put items on a queue [...],
        then that process will not terminate until all buffered items
        have been flushed to the pipe.
    """
    files = []
    for _ in workers:
        queue_in.put(None)
        files += queue_out.get()
        logging.info("now {:,} files".format(len(files)))

    # Join child-processes
    for l in workers:
        l.join()

    # Repeat until all files are loaded
    while files:
        logging.info("{:,} files to load".format(len(files)))
        queue_in = Queue()
        queue_out = Queue()
        workers = [
            DocumentLoader(_host, doc_type, queue_in, queue_out, suffix=suffix)
            for _ in range(min(processes, len(files)))
        ]
        for l in workers:
            l.start()

        for filepath in files:
            queue_in.put(filepath)

        files = []
        for _ in workers:
            queue_in.put(None)
            files += queue_out.get()
            logging.info("now {:,} files".format(len(files)))

        for l in workers:
            l.join()

    """
    Do NOT delete old indices:
    Old indices are those used in production, 
    we just want them not to be mapped to the `next` alias any more.
    They will be deleted when the new indices are mapped to the `current` alias
    """
    logging.info("creating temporary alias")
    update_alias(my_ippro, host, alias="next", suffix=suffix, delete=False)

    logging.info("complete")


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
                    except exceptions.ConnectionTimeout:
                        pass
                    except exceptions.NotFoundError:
                        break
                    else:
                        break
    else:
        # Create alias, and add it to new indices
        es.indices.put_alias(index=','.join(new_indices), name=alias)

    # Update index settings
    for index in new_indices:
        es.indices.put_settings({
            # 'number_of_replicas': 1,
            'refresh_interval': None  # default (1s)
        }, index)
