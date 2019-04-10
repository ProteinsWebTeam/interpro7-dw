#!/usr/bin/env python
# -*- coding: utf-8 -*-

import json
import os
import shutil
from multiprocessing import Process, Queue
from tempfile import mkstemp

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

        # for tax_id in xrefs.get("taxa", []):
        #     cross_refs.append({
        #         "dbname": "TAXONOMY",
        #         "dbkey": tax_id
        #     })

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


def write_json(uri: str, project_name: str, version: str, release_date: str,
               queue_in: Queue, chunk_size: int, queue_out: Queue):

    # Loading MySQL data
    entries = mysql.entry.get_entries(uri)
    entry2set = {
        entry_ac: set_ac
        for set_ac, s in mysql.entry.get_sets(uri).items()
        for entry_ac in s["members"]
    }
    databases = mysql.database.get_databases(uri)
    chunk = []

    while True:
        args = queue_in.get()
        if args is None:
            break

        acc, xrefs = args
        chunk.append(format_entry(entries.pop(acc), databases, xrefs,
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
                    "entries": chunk
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
                "entries": chunk
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
            os.chmod(outdir, mode=0o775)
            dir_count = 1

        filename = "{:0{}d}.json".format(dir_count, n_chars)
        dst = os.path.join(outdir, filename)
        shutil.move(src, dst)
        os.chmod(dst, mode=0o777)
        dir_count += 1


def dump(uri: str, src_entries: str, project_name: str, version: str,
         release_date: str, outdir: str, chunk_size: int=10,
         dir_limit: int=1000, processes: int=4,
         include_mobidblite: bool=False):
    logger.info("starting")

    # Create the directory (if needed), and remove its content
    os.makedirs(outdir, exist_ok=True)
    for item in os.listdir(outdir):
        path = os.path.join(outdir, item)
        try:
            os.remove(path)
        except IsADirectoryError:
            shutil.rmtree(path)

    processes = max(1, processes - 2)  # -2: parent process and organizer

    queue_entries = Queue(maxsize=processes*chunk_size)
    queue_files = Queue()
    writers = []
    for _ in range(processes):
        w = Process(target=write_json,
                    args=(uri, project_name, version, release_date,
                          queue_entries, chunk_size, queue_files))
        w.start()
        writers.append(w)

    organizer = Process(target=move_files,
                        args=(outdir, queue_files, dir_limit))
    organizer.start()

    entries = set(mysql.entry.get_entries(uri))
    n_entries = len(entries)
    cnt = 0
    with io.Store(src_entries) as store:
        for acc, xrefs in store:
            entries.remove(acc)

            if acc != "mobidb-lite" or include_mobidblite:
                queue_entries.put((acc, xrefs))

            cnt += 1
            if not cnt % 10000:
                logger.info("{:>8,} / {:>8,}".format(cnt, n_entries))

    # Remaining entries (without protein matches)
    for acc in entries:
        queue_entries.put((acc, None))

        cnt += 1
        if not cnt % 10000:
            logger.info("{:>8,} / {:>8,}".format(cnt, n_entries))

    logger.info("{:>8,} / {:>8,}".format(cnt, n_entries))

    for _ in writers:
        queue_entries.put(None)

    for w in writers:
        w.join()

    queue_files.put(None)
    organizer.join()

    logger.info("complete")


def dump_per_type(uri: str, src_entries: str, project_name: str, version: str,
                  release_date: str, outdir: str, chunk_size: int=50,
                  dir_limit: int=1000, processes: int=4,
                  include_mobidblite: bool=True):
    logger.info("starting")

    # Create the directory (if needed), and remove its content
    os.makedirs(outdir, exist_ok=True)
    for item in os.listdir(outdir):
        path = os.path.join(outdir, item)
        try:
            os.remove(path)
        except IsADirectoryError:
            shutil.rmtree(path)

    processes = max(1, processes - 2)  # -2: parent process and organizer

    queue_entries = Queue(maxsize=processes*chunk_size)
    queue_files = Queue()
    writers = []
    for _ in range(processes):
        w = Process(target=write_json_per_type,
                    args=(uri, project_name, version, release_date,
                          queue_entries, chunk_size, queue_files))
        w.start()
        writers.append(w)

    organizer = Process(target=move_files_per_type,
                        args=(outdir, queue_files, dir_limit))
    organizer.start()

    entries = set(mysql.entry.get_entries(uri))
    n_entries = len(entries)
    cnt = 0
    with io.Store(src_entries) as store:
        for acc, xrefs in store:
            entries.remove(acc)

            if acc != "mobidb-lite" or include_mobidblite:
                queue_entries.put((acc, xrefs))

            cnt += 1
            if not cnt % 10000:
                logger.info("{:>8,} / {:>8,}".format(cnt, n_entries))

    # Remaining entries (without protein matches)
    for acc in entries:
        queue_entries.put((acc, None))

        cnt += 1
        if not cnt % 10000:
            logger.info("{:>8,} / {:>8,}".format(cnt, n_entries))

    logger.info("{:>8,} / {:>8,}".format(cnt, n_entries))

    for _ in writers:
        queue_entries.put(None)

    for w in writers:
        w.join()

    queue_files.put(None)
    organizer.join()

    logger.info("complete")


def move_files_per_type(outdir: str, queue: Queue, dir_limit: int):
    type_dirs = {}
    n_chars = len(str(dir_limit))  # 1000 -> 4 chars -> 0001.json
    while True:
        args = queue.get()
        if args is None:
            break

        _type, src = args

        if _type in type_dirs:
            type_dir = type_dirs[_type]
        else:
            type_dir = type_dirs[_type] = {
                "dir": os.path.join(outdir, _type),
                "count": 1
            }
            os.mkdir(os.path.join(outdir, _type))

        if type_dir["count"] == dir_limit:
            dirname = "{:0{}d}".format(type_dir["count"], n_chars)
            type_dir.update({
                "dir": os.path.join(type_dir["dir"], dirname),
                "count": 1
            })
            os.mkdir(type_dir["dir"])

        filename = "{:0{}d}.json".format(type_dir["count"], n_chars)
        dst = os.path.join(type_dir["dir"], filename)
        shutil.move(src, dst)
        type_dir["count"] += 1


def write_json_per_type(uri: str, project_name: str, version: str,
                        release_date: str, queue_in: Queue, chunk_size: int,
                        queue_out: Queue):

    # Loading MySQL data
    entries = mysql.entry.get_entries(uri)
    entry2set = {
        entry_ac: set_ac
        for set_ac, s in mysql.entry.get_sets(uri).items()
        for entry_ac in s["members"]
    }
    databases = mysql.database.get_databases(uri)
    type_entries = {}

    while True:
        args = queue_in.get()
        if args is None:
            break

        acc, xrefs = args
        entry = entries.pop(acc)

        _type = entry["type"]
        if _type not in type_entries:
            type_entries[_type] = []

        type_entries[_type].append(format_entry(entry, databases, xrefs,
                                                entry2set.get(acc)))

        if len(type_entries[_type]) == chunk_size:
            fd, filepath = mkstemp()
            os.close(fd)

            with open(filepath, "wt") as fh:
                json.dump({
                    "name": project_name,
                    "release": version,
                    "release_date": release_date,
                    "entry_count": len(type_entries[_type]),
                    "entries": type_entries[_type]
                }, fh, indent=4)

            queue_out.put((_type, filepath))
            type_entries[_type] = []

    for _type, chunk in type_entries.items():
        if chunk:
            fd, filepath = mkstemp()
            os.close(fd)

            with open(filepath, "wt") as fh:
                json.dump({
                    "name": project_name,
                    "release": version,
                    "release_date": release_date,
                    "entry_count": len(chunk),
                    "entries": chunk
                }, fh, indent=4)

            queue_out.put((_type, filepath))
