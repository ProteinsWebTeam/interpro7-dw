# -*- coding: utf-8 -*-

from multiprocessing import Process, Queue
from tempfile import mkdtemp
from typing import Optional

from .organize import JsonFileOrganizer, set_ready
from .. import mysql
from ... import logger
from ...io import Store


def _create_doc(entry: dict, xrefs: Optional[dict]=None,
                set_acc: Optional[str]=None) -> dict:

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

    """
    DocumentLoader:
        - expects an `id` field
        - indexes the document in the `SRCH_INDEX` if no `entry_db` is found
    """
    return {
        "id": entry["accession"],
        "database": entry["database"],
        "type": entry["type"],
        "name": entry["name"],
        "references": list(refs)
    }


def _create_docs(uri: str, task_queue: Queue, outdir: str,
                 max_references: int=1000000):

    # Disable `items_per_file`
    organizer = JsonFileOrganizer(mkdtemp(dir=outdir), items_per_file=0)

    # Loading MySQL data
    entries = mysql.entry.get_entries(uri)
    entry2set = {
        entry_ac: set_ac
        for set_ac, s in mysql.entry.get_sets(uri).items()
        for entry_ac in s["members"]
    }

    num_references = 0

    for acc, xrefs in iter(task_queue.get, None):
        doc = _create_doc(entries.pop(acc), xrefs, entry2set.get(acc))
        organizer.add(doc)
        num_references += len(doc["references"])

        if num_references >= max_references:
            organizer.flush()
            num_references = 0

    organizer.flush()


def create_documents(uri: str, src_entries: str, outdir: str,
                     processes: int=4, include_mobidblite: bool=False):
    logger.info("starting")
    processes = max(1, processes - 2)  # -2: parent process and organizer

    task_queue = Queue()
    workers = []
    for _ in range(processes):
        w = Process(target=_create_docs, args=(uri, task_queue, outdir))
        w.start()
        workers.append(w)

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

    for w in workers:
        w.join()

    set_ready(outdir)

    logger.info("complete")
