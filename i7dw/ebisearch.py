#!/usr/bin/env python
# -*- coding: utf-8 -*-

import json
import logging
import os
import shutil
from multiprocessing import Process, Queue
from tempfile import mkstemp

from . import interpro, io

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s: %(levelname)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)


def prepare_entry(entry: dict, databases: dict, xrefs: dict=None,
                  set_ac=None) -> tuple:
    database = entry["database"]
    fields = {
        "id": entry["accession"].upper(),
        "short_name": entry["short_name"],
        "name": entry["name"],
        "type": entry["type"],
        "creation_date": entry["date"].strftime("%Y-%m-%d"),
        "source_database": databases[database]["name_long"],
        "description": " ".join(entry["descriptions"])
    }

    if set_ac:
        fields["set"] = set_ac

    cross_refs = {
        "UNIPROT": set(),
        "TAXONOMY": set(),
        "PROTEOMES": set(),
        "PDB": set()
    }
    if database == "interpro":
        dbs = fields["contributing_database"] = set()
        for dbname, dbkeys in entry["member_databases"].items():
            dbs.add(databases[dbname]["name_long"])

            s = cross_refs[dbname.upper()] = set()
            for dbkey in dbkeys:
                s.add(dbkey)

        for dbname, dbkeys in entry["cross_references"].items():
            k = dbname.upper()
            if k in cross_refs:
                s = cross_refs[k]
            else:
                s = cross_refs[k] = set()

            for dbkey in dbkeys:
                s.add(dbkey)

        s = cross_refs["PUBMED"] = set()
        for pub in entry["citations"].values():
            if pub.get("PMID"):
                s.add(pub["PMID"])

        s = cross_refs["GO"] = set()
        for term in entry["go_terms"]:
            s.add(term["identifier"])

        if "INTERPRO" in cross_refs:
            s = cross_refs["INTERPRO"]
        else:
            s = cross_refs["INTERPRO"] = set()
        for acc in entry["relations"]:
            s.add(acc)
    else:
        # Member DB signature
        if entry["integrated"]:
            cross_refs["INTERPRO"] = entry["integrated"].upper()

        s = cross_refs["PUBMED"] = set()
        for pub in entry["citations"].values():
            if pub.get("PMID"):
                s.add(pub["PMID"])

    if xrefs:
        for (protein_ac, protein_id) in xrefs.get("proteins", []):
            cross_refs["UNIPROT"].add(protein_ac)
            cross_refs["UNIPROT"].add(protein_id)

        for tax_id in xrefs.get("taxa", []):
            cross_refs["TAXONOMY"].add(tax_id)

        for upid in xrefs.get("proteomes", []):
            cross_refs["PROTEOMES"].add(upid)

        for pdbe_id in xrefs.get("structures", []):
            cross_refs["PDB"].add(pdbe_id)

    return fields, cross_refs


def format_entry(fields: dict, cross_refs: dict) -> dict:
    _fields = []
    for k, v in fields.items():
        if isinstance(v, set):
            for item in v:
                _fields.append({
                    "name": k,
                    "value": item
                })
        else:
            _fields.append({
                "name": k,
                "value": v
            })

    _cross_refs = []
    for dbname, dbkeys in cross_refs.items():
        for dbkey in dbkeys:
            _cross_refs.append({
                "dbname": dbname,
                "dbkey": dbkey
            })

    return {
        "fields": _fields,
        "cross_references": _cross_refs
    }


def write_json(uri: str, project_name: str, version: str, release_date: str,
               queue_in: Queue, chunk_size: int, queue_out: Queue):

    # Loading MySQL data
    entries = interpro.get_entries(uri)
    entry2set = {
        entry_ac: set_ac
        for set_ac, s in interpro.get_sets(uri).items()
        for entry_ac in s["members"]
    }
    databases = interpro.get_entry_databases(uri)
    chunk = []

    while True:
        args = queue_in.get()
        if args is None:
            break

        acc, xrefs = args
        chunk.append(prepare_entry(entries.pop(acc), databases, xrefs,
                                   entry2set.get(acc)))

        if len(chunk) == chunk_size:
            fd, filepath = mkstemp()
            os.close(fd)

            with open(filepath, "wt") as fh:
                json.dump({
                    "name": project_name,
                    "release": version,
                    "release_date": release_date,
                    "entry_count": len(chunk),
                    "entries": [format_entry(*e) for e in chunk]
                }, fh, indent=4)

            queue_out.put(filepath)

            chunk = []

    if chunk:
        fd, filepath = mkstemp()
        os.close(fd)

        with open(filepath, "wt") as fh:
            json.dump({
                "name": project_name,
                "release": version,
                "release_date": release_date,
                "entry_count": len(chunk),
                "entries": [format_entry(*e) for e in chunk]
            }, fh, indent=4)

        queue_out.put(filepath)


def move_files(outdir: str, queue: Queue, dir_limit: int):
    dir_count = 1
    n_chars = len(str(dir_limit))  # 1000 -> 4 chars -> 0001.json
    while True:
        src = queue.get()
        if src is None:
            break

        if dir_count == dir_limit:
            dirname = "{:0{}d}".format(dir_count, n_chars)
            outdir = os.path.join(outdir, dirname)
            os.mkdir(outdir)
            dir_count = 1

        filename = "{:0{}d}.json".format(dir_count, n_chars)
        dst = os.path.join(outdir, filename)
        shutil.move(src, dst)
        dir_count += 1


def dump(uri: str, src_entries: str, project_name: str, version: str,
         release_date: str, outdir: str, chunk_size: int=50,
         dir_limit: int=1000, n_readers: int=3, n_writers=3):
    logging.info("starting")

    # Create the directory (if needed), and remove its content
    os.makedirs(outdir, exist_ok=True)
    for item in os.listdir(outdir):
        path = os.path.join(outdir, item)
        try:
            os.remove(path)
        except IsADirectoryError:
            os.rmdir(path)

    queue_entries = Queue(maxsize=n_writers*chunk_size)
    queue_files = Queue()
    writers = [
        Process(target=write_json, args=(uri, project_name, version,
                                         release_date, queue_entries,
                                         chunk_size, queue_files))
        for _ in range(n_writers)
    ]
    for w in writers:
        w.start()

    organizer = Process(target=move_files,
                        args=(outdir, queue_files, dir_limit))
    organizer.start()

    entries = set(interpro.get_entries(uri))
    n_entries = len(entries)
    cnt = 0
    with io.Store(src_entries) as store:
        if n_readers > 1:
            fn = store.iter(n_readers)
        else:
            fn = store

        for acc, xrefs in fn:
            entries.remove(acc)
            queue_entries.put((acc, xrefs))

            cnt += 1
            if not cnt % 1000:
                logging.info("{:>8,} / {:>8,}".format(cnt, n_entries))

    # Remaining entries (without protein matches)
    for acc in entries:
        queue_entries.put((acc, None))

        cnt += 1
        if not cnt % 1000:
            logging.info("{:>8,} / {:>8,}".format(cnt, n_entries))

    logging.info("{:>8,} / {:>8,}".format(cnt, n_entries))

    for _ in writers:
        queue_entries.put(None)

    for w in writers:
        w.join()

    queue_files.put(None)
    organizer.join()

    logging.info("complete")
