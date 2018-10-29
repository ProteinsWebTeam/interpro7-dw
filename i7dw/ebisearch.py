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


def format_entry(entry: dict, databases: dict, xrefs: dict=None,
                 set_ac=None) -> dict:
    database = entry["database"]
    fields = [
        {
            "name": "id",
            "value": entry["accession"].upper()},
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
                    "dbname": dbname.upper(),
                    "dbkey": dbkey
                })

        for dbname, dbkeys in entry["cross_references"].items():
            for dbkey in dbkeys:
                cross_refs.append({
                    "dbname": dbname.upper(),
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
                "dbkey": entry["integrated"].upper()
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


def write_json(project_name: str, version: str, release_date: str,
               queue_in: Queue, queue_out: Queue):
    while True:
        chunk = queue_in.get()
        if chunk is None:
            break

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
            dir_count = 1

        filename = "{:0{}d}.json".format(dir_count, n_chars)
        dst = os.path.join(outdir, filename)
        shutil.move(src, dst)
        dir_count += 1


def dump(uri: str, src_entries: str, project_name: str, version: str,
         release_date: str, outdir: str, chunk_size: int=50,
         dir_limit: int=1000, n_readers: int=3, n_writers=3):

    logging.info("preparing output directory")

    # Create the directory (if needed), and remove its content
    os.makedirs(outdir, exist_ok=True)
    for item in os.listdir(outdir):
        path = os.path.join(outdir, item)
        try:
            os.remove(path)
        except IsADirectoryError:
            os.rmdir(path)

    queue_chunks = Queue()
    queue_files = Queue()
    writers = [
        Process(target=write_json, args=(project_name, version, release_date,
                                         queue_chunks, queue_files))
        for _ in range(n_writers)
    ]
    for w in writers:
        w.start()

    organizer = Process(target=move_files,
                        args=(outdir, queue_files, dir_limit))
    organizer.start()

    logging.info("loading MySQL data")
    entries = interpro.get_entries(uri)
    entry2set = {
        entry_ac: set_ac
        for set_ac, s in interpro.get_sets(uri).items()
        for entry_ac in s["members"]
    }
    databases = interpro.get_entry_databases(uri)

    logging.info("starting")
    n_entries = len(entries)
    cnt = 0
    chunk = []
    with io.Store(src_entries) as store:
        if n_readers > 1:
            fn = store.iter(n_readers)
        else:
            fn = store

        for acc, xrefs in fn:
            entry = entries.pop(acc)

            chunk.append(format_entry(entry, databases, xrefs,
                                      entry2set.get(acc)))

            if len(chunk) == chunk_size:
                queue_chunks.put(chunk)
                chunk = []

            cnt += 1
            if not cnt % 1000:
                logging.info("{:>6} / {}".format(cnt, n_entries))

    # Remaining entries (without protein matches)
    for acc, entry in entries.items():
        chunk.append(format_entry(entry, databases, None, entry2set.get(acc)))

        if len(chunk) == chunk_size:
            queue_chunks.put(chunk)
            chunk = []

        cnt += 1
        if not cnt % 1000:
            logging.info("{:>6} / {}".format(cnt, n_entries))

    if chunk:
        queue_chunks.put(chunk)

    logging.info("{:>6} / {}".format(cnt, n_entries))

    for _ in writers:
        queue_chunks.put(None)

    for w in writers:
        w.join()

    queue_files.put(None)
    organizer.join()

    logging.info("complete")
