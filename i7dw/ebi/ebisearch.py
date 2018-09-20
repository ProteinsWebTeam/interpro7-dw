#!/usr/bin/env python
# -*- coding: utf-8 -*-

import json
import logging
import math
import os
import time

from .. import disk, mysql


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s: %(levelname)s: %(message)s",
    datefmt="%y-%m-%d %H:%M:%S"
)


class JsonWriter(object):
    def __init__(self, name, version, release_date, dst, **kwargs):
        self.name = name
        self.version = version
        self.date = release_date
        self.dir = dst
        self.max_dir_size = kwargs.get("dir_size", 1000)
        self.cur_dir_size = 0

    def write(self, entries):
        # number of characters in files/dirs name
        n_chars = math.ceil(math.log10(self.max_dir_size)) + 1

        if self.cur_dir_size == self.max_dir_size:
            # Too many files in current directory

            self.dir = os.path.join(
                self.dir,
                "{:0{}d}".format(self.cur_dir_size + 1, n_chars)
            )
            self.cur_dir_size = 0
            os.mkdir(self.dir)

        self.cur_dir_size += 1
        filepath = os.path.join(
            self.dir,
            "{:0{}d}.json".format(self.cur_dir_size, n_chars)
        )
        with open(filepath, "wt") as fh:
            json.dump({
                "name": self.name,
                "release": self.version,
                "release_date": self.date,
                "entry_count": len(entries),
                "entries": entries
            }, fh, indent=4)


def export(uri, proteins_f, prot_matches_f, struct_matches_f, proteomes_f,
           name, version, release_date, dst, **kwargs):
    tmpdir = kwargs.get("tmpdir")

    # max number of cross-references in memory (soft limit)
    max_xref = kwargs.get("max_xref", 50000000)

    # max number of files/subdirs per directory
    dir_size = kwargs.get("dir_size", 1000)

    # max number of cross-references per JSON file
    json_size = kwargs.get("json_size", 4000000)

    if dir_size < 2:
        raise ValueError("dir_size cannot be lesser than 2")

    if tmpdir and not os.path.isdir(tmpdir):
        os.makedirs(tmpdir)

    writer = JsonWriter(name, version, release_date, dst, dir_size=dir_size)

    logging.info("loading data from MySQL")
    entries_info = mysql.get_entries(uri)
    sets = mysql.get_sets(uri)
    databases = mysql.get_entry_databases(uri)

    accessions = sorted(entries_info.keys())
    bucket_size = int(math.ceil(len(accessions) / dir_size))

    attic = disk.Attic(
        workdir=tmpdir,
        # List of maximum `dir_size` elements,
        # to determine in which bucket an entry should be stored
        accessions=[
            accessions[i].upper()
            for i in range(0, len(accessions), bucket_size)
        ],
        max_xref=max_xref,
        persist=False
    )

    proteins = disk.Store(proteins_f)
    prot_matches = disk.Store(prot_matches_f)
    struct_matches = disk.Store(struct_matches_f)
    proteomes = disk.Store(proteomes_f)

    n_proteins = 0
    ts = time.time()
    logging.info("storing protein-related cross-references")
    for accession, protein in proteins.iter():
        tax_id = protein["taxon"]
        _prot_matches = prot_matches.get(accession, [])

        """
        Flatten structures (dict to list)
        {feature/prediction -> dbname -> acc} to [(dbname ,acc), ...]
        """
        _struct_matches = struct_matches.get(accession, {})
        _struct_matches = [
            (dbname.upper(), acc.upper())
            for _type in _struct_matches
            for dbname in _struct_matches[_type]
            for acc in _struct_matches[_type][dbname]
        ]

        _proteomes = set(map(str.upper, proteomes.get(accession, [])))

        entries = set()
        for m in _prot_matches:
            method_ac = m["method_ac"].upper()
            entry_ac = m["entry_ac"]

            entries.add(method_ac)
            if entry_ac:
                entries.add(entry_ac.upper())

        entries_xref = {}
        for entry_ac in entries:
            e = entries_xref[entry_ac] = {
                "TAXONOMY": {tax_id},
                "UNIPROT": {accession},
                "PROTEOME": _proteomes
            }

            for dbname, acc in _struct_matches:
                if dbname in e:
                    e[dbname].add(acc)
                else:
                    e[dbname] = {acc}

        attic.put(entries_xref)
        n_proteins += 1
        if not n_proteins % 1000000:
            logging.info(
                "{:>12} ({:.0f} proteins/sec)".format(
                    n_proteins, n_proteins // (time.time() - ts)
                )
            )

    attic.close()

    for store in (proteins, prot_matches, struct_matches, proteomes):
        store.close()

    logging.info(
        "{:>12} ({:.0f} proteins/sec)".format(
            n_proteins, n_proteins // (time.time() - ts)
        )
    )

    if not os.path.isdir(dst):
        os.makedirs(dst)

    chunk = []          # entries to write in the current chunk
    n_xref = 0          # number of cross-references in the current chunk
    total_xref = 0      # total number of cross-references
    logging.info("dumping entries to JSON files")
    for accession in sorted(entries_info):
        entry = entries_info[accession]
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

        xref = {
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
                    if dbname in xref:
                        xref[dbname].add(acc.upper())
                    else:
                        xref[dbname] = {acc.upper()}

            for dbname in entry["cross_references"]:
                for acc in entry["cross_references"][dbname]:
                    if dbname in xref:
                        xref[dbname].add(acc.upper())
                    else:
                        xref[dbname] = {acc.upper()}

            for pub in entry["citations"].values():
                if pub.get("PMID"):
                    xref["PUBMED"].add(pub["PMID"])

            for term in entry["go_terms"]:
                xref["GO"].add(term["identifier"])

            for acc in entry["relations"]:
                xref["INTERPRO"].add(acc.upper())
        else:
            # Member DB signature
            if entry["integrated"]:
                xref["INTERPRO"] = {entry["integrated"].upper()}

            for pub in entry["citations"].values():
                if pub.get("PMID"):
                    xref["PUBMED"].add(pub["PMID"])

        # Read protein-related cross-references
        for dbname, dbkeys in attic.get(accession.upper()):
            xref[dbname] = dbkeys

        # Transform into EBI Search JSON data format
        l_fields = []
        for name in fields:
            if isinstance(fields[name], (list, tuple)):
                for value in fields[name]:
                    l_fields.append({"name": name, "value": value})
            else:
                l_fields.append({"name": name, "value": fields[name]})

        l_xref = [
            {"dbname": dbname, "dbkey": dbkey}
            for dbname in xref
            for dbkey in xref[dbname]
        ]

        if n_xref + len(l_xref) >= json_size:
            """
            Cannot add this entry without exceeding ratio
            of number of cross-ref / JSON file
                -> dump chunk
            """
            writer.write(chunk)
            chunk = []
            n_xref = 0

        # Add entry to chunk
        chunk.append({
            "fields": l_fields,
            "cross_references": l_xref
        })

        # Update counters
        n_xref += len(l_xref)
        total_xref += len(l_xref)

    if chunk:
        writer.write(chunk)

    logging.info(
        "{} cross-references for {} entries".format(
            total_xref, len(entries_info)
        )
    )
