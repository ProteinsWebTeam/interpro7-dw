# -*- coding: utf-8 -*-

import os
import shutil
from multiprocessing import Process, Queue
from tempfile import mkdtemp
from typing import Optional

from i7dw import io, logger
from i7dw.interpro import mysql


class XrefsWriter(Process):
    def __init__(self, url: str, version: str, release_date: str,
                 outdir: str, inqueue: Queue, outqueue: Queue,
                 buffer_size: int = 100000):
        super().__init__()
        self.url = url
        self.version = version
        self.release_date = release_date
        self.outdir = outdir
        self.inqueue = inqueue
        self.outqueue = outqueue
        self.buffer_size = buffer_size
        self.entries = None
        self.entry_set = None
        self.databases = None
        self.taxa = None

    def run(self):
        # Loading MySQL data
        self.entries = mysql.entries.get_entries(self.url)

        self.entry_set = {}
        for s in mysql.entries.iter_sets(self.url):
            set_acc = s["accession"]

            for entry_acc in s["members"]:
                self.entry_set[entry_acc] = set_acc

        self.databases = {}
        for name, info in mysql.databases.get_databases(self.url).items():
            self.databases[name] = info["name_long"]

        self.taxa = {}
        for taxon in mysql.taxonomy.iter_taxa(self.url, lineage=False):
            self.taxa[taxon["id"]] = taxon["scientific_name"]

        organizers = {}
        counters = {}
        total_references = 0
        for entry_acc, xrefs in iter(self.inqueue.get, None):
            entry = self.entries.pop(entry_acc)
            entry_type = entry["type"]
            entry = self.format_entry(entry, xrefs=xrefs,
                                      set_acc=self.entry_set.get(entry_acc))
            cnt_references = len(entry["cross_references"])
            total_references += cnt_references

            try:
                counters[entry_type] += cnt_references
            except KeyError:
                counters[entry_type] = cnt_references
                type_dir = os.path.join(self.outdir, entry_type)
                os.makedirs(type_dir, exist_ok=True)

                """
                Create a subdirectory 
                (so each process writes to its own subdirectory)
                """
                outdir = mkdtemp(dir=type_dir)

                """
                `items_per_file=0` disables the auto-flush:
                we do not want to have a fixed number of entries per file,
                we want to avoid having too many cross-references in one file
                """
                jfo = io.JsonFileOrganizer(outdir, items_per_file=0,
                                           func=self.wrap, indent=4)
                organizers[entry_type] = jfo

            organizers[entry_type].add(entry)
            if counters[entry_type] >= self.buffer_size:
                organizers[entry_type].flush()
                counters[entry_type] = 0

        for jfo in organizers.values():
            jfo.flush()

        self.outqueue.put(total_references)

    def wrap(self, entries: list) -> dict:
        return {
            "name": "InterPro",
            "release": self.version,
            "release_date": self.release_date,
            "entry_count": len(entries),
            "entries": entries
        }

    def format_entry(self, entry: dict, xrefs: dict,
                     set_acc: Optional[str]=None) -> dict:
        fields = [
            {
                "name": "id",
                "value": entry["accession"]},
            {
                "name": "short_name",
                "value": entry["short_name"] or entry["accession"]},
            {
                "name": "name",
                "value": entry["name"] or entry["accession"]},
            {
                "name": "type",
                "value": entry["type"]},
            {
                "name": "creation_date",
                "value": entry["date"].strftime("%Y-%m-%d")},
            {
                "name": "source_database",
                "value": self.databases[entry["database"]]
            },
            {
                "name": "description",
                "value": " ".join(entry["descriptions"])
            }
        ]

        if set_acc:
            fields.append({"name": "set", "value": set_acc})

        cross_refs = []
        if entry["database"] == "interpro":
            for dbname, dbkeys in entry["member_databases"].items():
                fields.append({
                    "name": "contributing_database",
                    "value": self.databases[dbname]
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

        for protein_ac, protein_id in xrefs.get("proteins", []):
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

            cross_refs.append({
                "dbname": "TAXONOMY",
                "dbkey": self.taxa[tax_id]
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


def dump(url: str, src_entries: str, version: str,
         release_date: str, outdir: str, processes: int=4):
    logger.info("starting")

    # Create the directory (if needed), and remove its content
    os.makedirs(outdir, exist_ok=True)
    for item in os.listdir(outdir):
        shutil.rmtree(os.path.join(outdir, item))

    task_queue = Queue(maxsize=1)
    done_queue = Queue()
    writers = []
    for _ in range(max(1, processes-1)):
        p = XrefsWriter(url, version, release_date, outdir,
                        task_queue, done_queue)
        p.start()
        writers.append(p)

    cnt_entries = 0
    with io.Store(src_entries) as store:
        for acc, xrefs in store:

            # TODO: create global variable somewhere for mobidb-lite
            if acc != "mobidb-lite":
                task_queue.put((acc, xrefs))

            cnt_entries += 1
            if not cnt_entries % 10000:
                logger.info(f"{cnt_entries:>8,}")

    for _ in writers:
        task_queue.put(None)

    cnt_xrefs = sum([done_queue.get() for _ in writers])

    for p in writers:
        p.join()

    logger.info(f"{cnt_entries:>8,} ({cnt_xrefs:,} cross-references)")


def exchange(src: str, dst: str):
    try:
        shutil.rmtree(dst)
    except FileNotFoundError:
        pass
    finally:
        shutil.move(src, dst)
