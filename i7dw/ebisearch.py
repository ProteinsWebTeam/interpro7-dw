#!/usr/bin/env python
# -*- coding: utf-8 -*-

import json
import logging
import os

from . import interpro, io

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s: %(levelname)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)


def dump(uri: str, src_entries: str, project_name: str, version: str,
         release_date: str, outdir: str, chunk_size: int=100,
         dir_limit: int=1000):

    logging.info("starting")

    # Create the directory (if needed), and remove its content
    os.makedirs(outdir, exist_ok=True)
    for item in os.listdir(outdir):
        path = os.path.join(outdir, item)
        try:
            os.remove(path)
        except IsADirectoryError:
            os.rmdir(path)

    entries = interpro.get_entries(uri)
    entry2set = {
        entry_ac: set_ac
        for set_ac, s in interpro.get_sets(uri)
        for entry_ac in s["members"]
    }

    databases = interpro.get_entry_databases(uri)

    i = 0
    dir_count = 1
    chunk = []
    n_chars = len(str(dir_limit))  # 1000 -> 4 chars -> 0001.json
    with io.Store(src_entries) as store:
        for i, acc in enumerate(sorted(entries)):
            e = entries[acc]
            database = e["database"]

            fields = [
                {
                    "name": "id",
                    "value": acc.upper()},
                {
                    "name": "short_name",
                    "value": e["short_name"]},
                {
                    "name": "name",
                    "value": e["name"]},
                {
                    "name": "type",
                    "value": e["type"]},
                {
                    "name": "creation_date",
                    "value": e["date"].strftime("%Y-%m-%d")},
                {
                    "name": "source_database",
                    "value": databases[database]["name_long"]},
                {
                    "name": "description",
                    "value": " ".join(e["descriptions"])
                }
            ]

            if acc in entry2set:
                fields.append({
                    "name": "set",
                    "value": entry2set[acc]
                })

            cross_refs = []
            if database == "interpro":
                for dbname, dbkeys in e["member_databases"].items():
                    fields.append({
                        "name": "contributing_database",
                        "value": databases[dbname]["name_long"]
                    })

                    for dbkey in dbkeys:
                        cross_refs.append({
                            "dbname": dbname.upper(),
                            "dbkey": dbkey
                        })

                for dbname, dbkeys in e["cross_references"].items():
                    for dbkey in dbkeys:
                        cross_refs.append({
                            "dbname": dbname.upper(),
                            "dbkey": dbkey
                        })

                for pub in e["citations"].values():
                    if pub.get("PMID"):
                        cross_refs.append({
                            "dbname": "PUBMED",
                            "dbkey": pub["PMID"]
                        })

                for term in e["go_terms"]:
                    cross_refs.append({
                        "dbname": "GO",
                        "dbkey": term["identifier"]
                    })

                for acc2 in e["relations"]:
                    cross_refs.append({
                        "dbname": "INTERPRO",
                        "dbkey": acc2
                    })
            else:
                # Member DB signature
                if e["integrated"]:
                    cross_refs.append({
                        "dbname": "INTERPRO",
                        "dbkey": e["integrated"].upper()
                    })

                for pub in e["citations"].values():
                    if pub.get("PMID"):
                        cross_refs.append({
                            "dbname": "PUBMED",
                            "dbkey": pub["PMID"]
                        })

            data = store.get(acc, {})
            for (protein_ac, protein_id) in data.get("proteins", []):
                cross_refs.append({
                    "dbname": "UNIPROT",
                    "dbkey": protein_ac
                })

                cross_refs.append({
                    "dbname": "UNIPROT",
                    "dbkey": protein_id
                })

            for tax_id in data.get("taxa", []):
                cross_refs.append({
                    "dbname": "TAXONOMY",
                    "dbkey": tax_id
                })

            for upid in data.get("proteomes", []):
                cross_refs.append({
                    "dbname": "PROTEOMES",
                    "dbkey": upid
                })

            for pdbe_id in data.get("structures", []):
                cross_refs.append({
                    "dbname": "PDB",
                    "dbkey": pdbe_id
                })

            chunk.append({
                "fields": fields,
                "cross_references": cross_refs
            })

            if len(chunk) == chunk_size:
                filename = "{:0{}d}.json".format(dir_count, n_chars)
                with open(os.path.join(outdir, filename), "wt") as fh:
                    json.dump({
                        "name": project_name,
                        "release": version,
                        "release_date": release_date,
                        "entry_count": len(chunk),
                        "entries": chunk
                    }, fh, indent=4)

                chunk = []
                dir_count += 1

                if dir_count == dir_limit:
                    dirname = "{:0{}d}".format(dir_count, n_chars)
                    outdir = os.path.join(outdir, dirname)
                    os.mkdir(outdir)
                    dir_count = 1

                logging.info("{:>6} / {}".format(i+1, len(entries)))

    if chunk:
        filename = "{:0{}d}.json".format(dir_count, n_chars)
        with open(os.path.join(outdir, filename), "wt") as fh:
            json.dump({
                "name": project_name,
                "release": version,
                "release_date": release_date,
                "entry_count": len(chunk),
                "entries": chunk
            }, fh, indent=4)

    logging.info("{:>6} / {}".format(i+1, len(entries)))
