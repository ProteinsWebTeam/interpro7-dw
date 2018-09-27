#!/usr/bin/env python
# -*- coding: utf-8 -*-

import json
import logging
import math
import os
import shutil
import time
from multiprocessing import Process, Queue, Value

from .. import disk, mysql

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s: %(levelname)s: %(message)s",
    datefmt="%y-%m-%d %H:%M:%S"
)


class XrefWriter(Process):
    def __init__(self, uri: str, name: str, version: str, release_date: str,
                 root: str, queue: Queue, outdir: str, n_files: Value,
                 **kwargs):
        super().__init__()
        self.uri = uri
        self.name = name
        self.version = version
        self.date = release_date
        self.xrefs = disk.XrefStore(root=root)
        self.queue = queue
        self.outdir = outdir
        self.n_files = n_files
        self.max_files = kwargs.get("max_files", 1000)
        self.max_xrefs = kwargs.get("max_xrefs", 1000000)

    def run(self):
        logging.info("{} ({}) started".format(self.name, os.getpid()))
        sets = mysql.get_sets(self.uri)
        databases = mysql.get_entry_databases(self.uri)

        # Keep one set per entry (can entries belong to more than one set?)
        sets = {acc: sets[acc].keys()[0] for acc in sets}

        entries = []
        n_entries = 0   # total number of processed entries
        n_xrefs = 0     # total number of cross-references
        mem_xrefs = 0   # number of cross-references in memory
        while True:
            entry = self.queue.get()
            if entry is None:
                break

            n_entries += 1
            accession = entry["accession"]
            if entry["name"]:
                name = entry["name"]
            else:
                name = entry["short_name"]

            fields = {
                "id": accession.upper(),
                "short_name": entry["short_name"],
                "name": name,
                "type": entry["type"],
                "creation_date": entry["date"].strftime("%Y-%m-%d"),
                "source_database": entry["database"],
                "description": " ".join(entry["descriptions"]),
                "set": sets.get(accession)
            }

            entry_xrefs = {
                "INTERPRO": set(),
                "PUBMED": set(),
                "GO": set()
            }

            """
            Add cross-references that are not related to proteins 
            (so could not be added before)
            """
            if entry["member_databases"]:
                # InterPro entry

                # Convert database names (e.g. cathgene3d -> CATH-Gene3D)
                fields["contributing_database"] = [
                    databases[dbname]["name_long"]
                    for dbname in entry["member_databases"]
                ]

                for dbname in entry["member_databases"]:
                    for acc in entry["member_databases"][dbname]:
                        if dbname in entry_xrefs:
                            entry_xrefs[dbname].add(acc.upper())
                        else:
                            entry_xrefs[dbname] = {acc.upper()}

                for dbname in entry["cross_references"]:
                    for acc in entry["cross_references"][dbname]:
                        if dbname in entry_xrefs:
                            entry_xrefs[dbname].add(acc.upper())
                        else:
                            entry_xrefs[dbname] = {acc.upper()}

                for pub in entry["citations"].values():
                    if pub.get("PMID"):
                        entry_xrefs["PUBMED"].add(pub["PMID"])

                for term in entry["go_terms"]:
                    entry_xrefs["GO"].add(term["identifier"])

                for acc in entry["relations"]:
                    entry_xrefs["INTERPRO"].add(acc.upper())
            else:
                # Member DB signature
                if entry["integrated"]:
                    entry_xrefs["INTERPRO"] = {entry["integrated"].upper()}

                for pub in entry["citations"].values():
                    if pub.get("PMID"):
                        entry_xrefs["PUBMED"].add(pub["PMID"])

            # Read protein-related cross-references
            for dbname, dbkeys in self.xrefs.get(accession.upper()):
                entry_xrefs[dbname] = dbkeys

            # Transform into EBI Search JSON data format
            l_fields = []
            for name in fields:
                if isinstance(fields[name], (list, tuple)):
                    for value in fields[name]:
                        l_fields.append({"name": name, "value": value})
                else:
                    l_fields.append({"name": name, "value": fields[name]})

            l_xrefs = [
                {"dbname": dbname, "dbkey": dbkey}
                for dbname in entry_xrefs
                for dbkey in entry_xrefs[dbname]
            ]

            if mem_xrefs + len(l_xrefs) >= self.max_xrefs:
                """
                Cannot add this entry without exceeding ratio
                of number of cross-ref / JSON file
                    -> dump entries
                """
                self.dump(entries)
                entries = []
                mem_xrefs = 0

            # Add entry
            entries.append({
                "fields": l_fields,
                "cross_references": l_xrefs
            })

            n_xrefs += len(l_xrefs)
            mem_xrefs += len(l_xrefs)

        self.dump(entries)
        logging.info(
            "{} ({}) terminated: "
            "{} cross-references in {} entries".format(
                self.name, os.getpid(), n_xrefs, n_entries
            )
        )

    def dump(self, entries: list):
        # number of characters in files/dirs name (e.g. 1000 -> 4)
        n_chars = math.ceil(math.log10(self.max_files)) + 1

        if self.n_files.value + 1 == self.max_files:
            # Too many files in current directory
            # +1 because one reserved for the new sub-dir

            with self.n_files.get_lock():
                self.outdir = os.path.join(
                    self.outdir,
                    "{:0{}d}".format(self.n_files.value+1, n_chars)
                )
                self.n_files.value = 0
                os.mkdir(self.outdir)

        with self.n_files.get_lock():
            self.n_files.value += 1
            filepath = os.path.join(
                self.outdir,
                "{:0{}d}.json".format(self.n_files.value, n_chars)
            )

        with open(filepath, "wt") as fh:
            json.dump({
                "name": self.name,
                "release": self.version,
                "release_date": self.date,
                "entry_count": len(entries),
                "entries": entries
            }, fh, indent=4)


def create_index(uri: str, proteins_f: str, prot_matches_f: str,
                 struct_matches_f: str, proteomes_f: str, name: str,
                 version: str, release_date: str, outdir: str, **kwargs):
    tmpdir = kwargs.get("tmpdir")

    # max number of cross-references in memory (soft limit)
    max_xrefs = kwargs.get("max_xrefs", 1000000)

    # max number of files/subdirs per directory
    max_files = kwargs.get("max_files", 1000)

    # max number of cross-references per JSON file
    json_size = kwargs.get("json_size", 1000000)

    limit = kwargs.get("limit", 0)

    n_writers = kwargs.get("writers", 1)

    if max_files < 2:
        raise ValueError("max_files cannot be lesser than 2")

    if tmpdir:
        try:
            os.makedirs(tmpdir)
        except FileExistsError:
            pass

    entries_info = mysql.get_entries(uri)
    xrefs = disk.XrefStore(
        accessions=sorted(entries_info.keys()),
        n_aisles=max_files, tmpdir=tmpdir, max_xrefs=max_xrefs
    )

    proteins_s = disk.Store(proteins_f)
    prot_matches_s = disk.Store(prot_matches_f)
    struct_matches_s = disk.Store(struct_matches_f)
    proteomes_s = disk.Store(proteomes_f)

    n_proteins = 0
    ts = time.time()
    logging.info("storing protein-related cross-references")
    for accession, protein in proteins_s.iter():
        tax_id = protein["taxon"]
        prot_matches = prot_matches_s.get(accession, [])

        """
        Flatten structures (dict to list)
        {feature/prediction -> dbname -> acc} to [(dbname ,acc), ...]
        """
        structures = set()
        for databases in struct_matches_s.get(accession, {}).values():
            for dbname, domains in databases.items():
                for _id, dom in domains.items():
                    structures.add((dbname.upper(), dom["domain_id"]))

        proteomes = set(map(str.upper, proteomes_s.get(accession, [])))

        entries = set()
        for m in prot_matches:
            method_ac = m["method_ac"]
            entry_ac = m["entry_ac"]

            entries.add(method_ac)
            if entry_ac:
                entries.add(entry_ac)

        entry_xrefs = {}
        for entry_ac in entries:
            e = entry_xrefs[entry_ac] = {
                "TAXONOMY": {tax_id},
                "UNIPROT": {accession, protein["identifier"]},
                "PROTEOME": proteomes
            }

            for dbname, acc in structures:
                if dbname in e:
                    e[dbname].add(acc)
                else:
                    e[dbname] = {acc}

        xrefs.add(entry_xrefs)
        n_proteins += 1

        if n_proteins == limit:
            break
        elif not n_proteins % 1000000:
            logging.info(
                "{:>12} ({:.0f} proteins/sec)".format(
                    n_proteins, n_proteins // (time.time() - ts)
                )
            )

    xrefs.clean()

    proteins_s.close()
    prot_matches_s.close()
    struct_matches_s.close()
    proteomes_s.close()

    logging.info(
        "{:>12} ({:.0f} proteins/sec)".format(
            n_proteins, n_proteins // (time.time() - ts)
        )
    )

    logging.info("temporary files: {} bytes".format(xrefs.get_size()))
    logging.info("dumping entries to JSON files")

    try:
        shutil.rmtree(outdir)
    except FileNotFoundError:
        pass
    finally:
        os.makedirs(outdir)

    n_files = Value("I", 0)
    queue = Queue()

    writers = [
        XrefWriter(uri, name, version, release_date, xrefs.root,
                   queue, outdir, n_files,
                   max_files=max_files, max_xrefs=json_size)
        for _ in range(n_writers)
    ]
    
    for w in writers:
        w.start()

    for accession in sorted(entries_info):
        queue.put(entries_info[accession])

    for _ in writers:
        queue.put(None)

    for w in writers:
        w.join()

    xrefs.close()
    logging.info("complete")
