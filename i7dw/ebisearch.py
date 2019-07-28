# -*- coding: utf-8 -*-

import os
import shutil
from multiprocessing import Process, Queue
from tempfile import mkdtemp

from . import io, logger
from .interpro import mysql


def format_entry(entry: dict, databases: dict, xrefs: dict=None,
                 set_ac: str=None) -> dict:
    database = entry["database"]
    fields = [
        {
            "name": "id",
            "value": entry["accession"]},
        {
            "name": "short_name",
            "value": entry["short_name"]},
        {
            "name": "name",
            "value": entry["name"]},
        {
            "name": "type",
            "value": entry["type"]},
        {
            "name": "creation_date",
            "value": entry["date"].strftime("%Y-%m-%d")},
        {
            "name": "source_database",
            "value": databases[database]["name_long"]
        },
        {
            "name": "description",
            "value": " ".join(entry["descriptions"])
        }
    ]

    if set_ac:
        fields.append({
            "name": "set",
            "value": set_ac
        })

    cross_refs = []
    if database == "interpro":
        for dbname, dbkeys in entry["member_databases"].items():
            fields.append({
                "name": "contributing_database",
                "value": databases[dbname]["name_long"]
            })

            for dbkey in dbkeys:
                cross_refs.append({
                    "dbname": dbname,
                    "dbkey": dbkey
                })

        for dbname, dbkeys in entry["cross_references"].items():
            for dbkey in dbkeys:
                cross_refs.append({
                    "dbname": dbname,
                    "dbkey": dbkey
                })

        for pub in entry["citations"].values():
            if pub.get("PMID"):
                cross_refs.append({
                    "dbname": "PUBMED",
                    "dbkey": pub["PMID"]
                })

        for term in entry["go_terms"]:
            cross_refs.append({
                "dbname": "GO",
                "dbkey": term["identifier"]
            })

        for acc in entry["relations"]:
            cross_refs.append({
                "dbname": "INTERPRO",
                "dbkey": acc
            })
    else:
        # Member DB signature
        if entry["integrated"]:
            cross_refs.append({
                "dbname": "INTERPRO",
                "dbkey": entry["integrated"]
            })

        for pub in entry["citations"].values():
            if pub.get("PMID"):
                cross_refs.append({
                    "dbname": "PUBMED",
                    "dbkey": pub["PMID"]
                })

    if xrefs:
        for (protein_ac, protein_id) in xrefs.get("proteins", []):
            cross_refs.append({
                "dbname": "UNIPROT",
                "dbkey": protein_ac
            })

            cross_refs.append({
                "dbname": "UNIPROT",
                "dbkey": protein_id
            })

        for tax_id in xrefs.get("taxa", []):
            cross_refs.append({
                "dbname": "TAXONOMY",
                "dbkey": tax_id
            })

        for upid in xrefs.get("proteomes", []):
            cross_refs.append({
                "dbname": "PROTEOMES",
                "dbkey": upid
            })

        for pdbe_id in xrefs.get("structures", []):
            cross_refs.append({
                "dbname": "PDB",
                "dbkey": pdbe_id
            })

    return {
        "fields": fields,
        "cross_references": cross_refs
    }


class JsonWrapper(object):
    def __init__(self, project: str, version: str, release: str):
        self.project = project
        self.version = version
        self.release = release

    def wrap(self, entries: list) -> dict:
        return {
            "name": self.project,
            "release": self.version,
            "release_date": self.release,
            "entry_count": len(entries),
            "entries": entries
        }


def _write(uri: str, outdir: str, task_queue: Queue, wrapper: JsonWrapper,
           done_queue: Queue, by_type: bool=False, max_references: int=100000):

    # Loading MySQL data
    entries = mysql.entry.get_entries(uri)
    entry2set = {
        entry_ac: set_ac
        for set_ac, s in mysql.entry.get_sets(uri).items()
        for entry_ac in s["members"]
    }
    databases = mysql.database.get_databases(uri)

    organizers = {}
    counters = {}
    if by_type:
        organizer = None
    else:
        organizer = io.JsonFileOrganizer(mkdtemp(dir=outdir),
                                         items_per_file=0,
                                         func=wrapper.wrap,
                                         indent=4)

    num_references = tot_references = 0
    for acc, xrefs in iter(task_queue.get, None):
        entry = entries.pop(acc)
        item = format_entry(entry, databases, xrefs, entry2set.get(acc))
        tot_references += len(item["cross_references"])

        if by_type:
            _type = entry["type"]
            if _type in organizers:
                counters[_type] += len(item["cross_references"])
            else:
                workdir = os.path.join(outdir, _type)
                os.makedirs(workdir, exist_ok=True)
                organizers[_type] = io.JsonFileOrganizer(mkdtemp(dir=workdir),
                                                         items_per_file=0,
                                                         func=wrapper.wrap,
                                                         indent=4)
                counters[_type] = len(item["cross_references"])

            organizers[_type].add(item)
            if counters[_type] >= max_references:
                organizers[_type].flush()
                counters[_type] = 0
        else:
            organizer.add(item)
            num_references += len(item["cross_references"])

            if num_references >= max_references:
                organizer.flush()
                num_references = 0

    if by_type:
        for organizer in organizers.values():
            organizer.flush()
    else:
        organizer.flush()
    done_queue.put(tot_references)


def dump(uri: str, src_entries: str, project_name: str, version: str,
         release_date: str, outdir: str, processes: int=4,
         by_type: bool=False):
    logger.info("starting")

    # Create the directory (if needed), and remove its content
    os.makedirs(outdir, exist_ok=True)
    for item in os.listdir(outdir):
        path = os.path.join(outdir, item)
        try:
            os.remove(path)
        except IsADirectoryError:
            shutil.rmtree(path)

    wrapper = JsonWrapper(project_name, version, release_date)
    task_queue = Queue(maxsize=1)
    done_queue = Queue()
    writers = []
    for _ in range(max(1, processes-1)):
        p = Process(target=_write,
                    args=(uri, outdir, task_queue, wrapper, done_queue,
                          by_type))
        p.start()
        writers.append(p)

    entries = set(mysql.entry.get_entries(uri))
    n_entries = len(entries)
    cnt = 0
    with io.Store(src_entries) as store:
        for acc, xrefs in store:
            entries.remove(acc)

            if acc != "mobidb-lite":
                task_queue.put((acc, xrefs))

            cnt += 1
            if not cnt % 10000:
                logger.info(f"{cnt:>8,}/{n_entries:,}")

    # Remaining entries (without protein matches)
    for acc in entries:
        task_queue.put((acc, None))

        cnt += 1
        if not cnt % 10000:
            logger.info(f"{cnt:>8,}/{n_entries:,}")

    for _ in writers:
        task_queue.put(None)

    n_refs = sum([done_queue.get() for _ in writers])

    for p in writers:
        p.join()

    logger.info(f"{cnt:>8,}/{n_entries:,} "
                f"({n_refs:,} cross-references)")


def exchange(src: str, dst: str):
    try:
        shutil.rmtree(dst)
    except FileNotFoundError:
        pass
    finally:
        shutil.move(src, dst)
