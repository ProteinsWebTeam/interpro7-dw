# -*- coding: utf-8 -*-

from multiprocessing import Process, Queue
from tempfile import mkdtemp
from typing import Generator, List, Optional

from . import index, set_ready
from .. import mysql
from ... import logger
from ...io import JsonFileOrganizer, Store


SRCH_INDEX = "iprsearch"


class DocumentProducer(Process):
    def __init__(self, my_url: str, task_queue: Queue, outdir: str,
                 max_references: int=1000000):
        super().__init__()
        self.url = my_url
        self.task_queue = task_queue
        self.max_references = max_references
        """
        Disable `items_per_file` because we want to flush manually
        (we do not care about the number of entry/item per file,
        but we do about the number of cross-references per file)
        """
        self.organizer = JsonFileOrganizer(outdir, items_per_file=0)

    def run(self):
        # Loading MySQL data
        entries = mysql.entry.get_entries(self.url)
        entry2set = {
            entry_ac: set_ac
            for set_ac, s in mysql.entry.get_sets(self.url).items()
            for entry_ac in s["members"]
        }

        num_references = 0
        for acc, xrefs in iter(self.task_queue.get, None):
            entry = entries.pop(acc)
            set_acc = entry2set.get(acc)
            gen = self.chunk_xrefs(entry, xrefs, set_acc, self.max_references)
            for chunk in gen:
                self.organizer.add(chunk)
                num_references += len(chunk["references"])

                if num_references >= self.max_references:
                    self.organizer.flush()
                    num_references = 0

        self.organizer.flush()

    @staticmethod
    def chunk_xrefs(entry: dict, xrefs: Optional[dict] = None,
                    set_acc: Optional[str] = None,
                    max_references: int = 0) -> Generator[dict, None, None]:

        refs = set()
        if entry["database"] == "interpro":
            for src_db, signatures in entry["member_databases"].items():
                refs.add(src_db)

                for acc, name in signatures.items():
                    refs.add(acc)
                    refs.add(name)

            for ref_db, ref_ids in entry["cross_references"].items():
                refs.add(ref_db)
                for ref_id in ref_ids:
                    refs.add(ref_id)

            for term in entry["go_terms"]:
                refs.add(term["identifier"])

            for entry_acc in entry["relations"]:
                refs.add(entry_acc)
        else:
            if entry["integrated"]:
                refs.add(entry["integrated"])

        for pub in entry["citations"].values():
            if pub.get("PMID"):
                refs.add(pub["PMID"])

        if xrefs:
            for protein_acc, protein_id in xrefs.get("proteins", []):
                refs.add(protein_acc)
                refs.add(protein_id)

            for tax_id in xrefs.get("taxa", []):
                refs.add(tax_id)

            for upid in xrefs.get("proteomes", []):
                refs.add(upid)

            for pdbe_id in xrefs.get("structures", []):
                refs.add(pdbe_id)

        if set_acc:
            refs.add(set_acc)

        if not refs:
            yield {
                "entry_acc": entry["accession"],
                "entry_db": entry["database"],
                "entry_type": entry["type"],
                "entry_name": entry["name"],
                "references": ""
            }
        elif max_references > 0:
            refs = list(refs)
            for i in range(0, len(refs), max_references):
                yield {
                    "entry_acc": entry["accession"],
                    "entry_db": entry["database"],
                    "entry_type": entry["type"],
                    "entry_name": entry["name"],
                    "references": ' '.join(
                        map(str, refs[i:i + max_references]))
                }
        else:
            yield {
                "entry_acc": entry["accession"],
                "entry_db": entry["database"],
                "entry_type": entry["type"],
                "entry_name": entry["name"],
                "references": ' '.join(map(str, refs))
            }


def create_documents(uri: str, src_entries: str, outdir: str,
                     processes: int=4, include_mobidblite: bool=False):
    logger.info("starting")
    processes = max(1, processes-1)  # minus one for parent process

    task_queue = Queue(processes)
    workers = []
    for _ in range(processes):
        p = DocumentProducer(uri, task_queue, mkdtemp(dir=outdir))
        p.start()
        workers.append(p)

    entries = set(mysql.entry.get_entries(uri))
    n_entries = len(entries)
    cnt = 0
    with Store(src_entries) as store:
        for acc, xrefs in store:
            entries.remove(acc)

            if acc != "mobidb-lite" or include_mobidblite:
                task_queue.put((acc, xrefs))

            cnt += 1
            if not cnt % 10000:
                logger.info("{:>8,} / {:>8,}".format(cnt, n_entries))

    # Remaining entries (without protein matches)
    for acc in entries:
        task_queue.put((acc, None))

        cnt += 1
        if not cnt % 10000:
            logger.info("{:>8,} / {:>8,}".format(cnt, n_entries))

    logger.info("{:>8,} / {:>8,}".format(cnt, n_entries))

    for _ in workers:
        task_queue.put(None)

    for p in workers:
        p.join()

    set_ready(outdir)
    logger.info("complete")


class DocumentController(index.DocumentController):
    def dump(self, doc: dict) -> dict:
        return {
            "_op_type": "update",
            "_index": SRCH_INDEX + self.suffix,
            "_type": "search",
            "_id": doc["entry_acc"],
            "_source": {
                "script": {
                    "source": "ctx._source.references += params.references",
                    "lang": "painless",
                    "params": {
                        # prefix by a whitespace to avoid string concatenation
                        "references": ' ' + doc["references"]
                    }
                },
                "upsert": doc
            }
        }

    def parse(self, item: dict):
        try:
            del item["update"]["data"]
        except KeyError:
            pass


def index_documents(hosts: List[str], src: str, **kwargs):
    indices = [SRCH_INDEX]

    if kwargs.get("body_path"):
        # Create indices
        index.create_indices(hosts=hosts,
                             indices=indices,
                             body_path=kwargs.pop("body_path"),
                             doc_type="search",
                             **kwargs
                             )

    controller = DocumentController(**kwargs)
    alias = kwargs.get("alias")
    if index.index_documents(hosts, controller, src, **kwargs) and alias:
        index.update_alias(hosts, indices, alias, **kwargs)


def update_alias(hosts: List[str], alias: str, **kwargs):
    index.update_alias(hosts, [SRCH_INDEX], alias, **kwargs)
