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

    logging.info("loading data from MySQL")
    entries_info = mysql.get_entries(uri)
    sets = mysql.get_sets(uri)
    databases = mysql.get_entry_databases(uri)

    accessions = sorted(entries_info.keys())
    bucket_size = math.ceil(len(accessions) / dir_size)

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

    cnt = 0
    total = 0
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

        _proteomes = list(map(str.upper, proteomes.get(accession, [])))

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
                "PROTEOME": set(_proteomes)
            }

            for dbname, acc in _struct_matches:
                if dbname in e:
                    e[dbname].add(acc)
                else:
                    e[dbname] = {acc}

        attic.put(entries_xref)
        cnt += 1
        total += 1
        if not total % 1000000:
            logging.info(
                "{:>12} ({:.0f} proteins/sec)".format(
                    total, cnt // (time.time() - ts)
                )
            )
            cnt = 0
            ts = time.time()

    attic.close()

    for store in (proteins, prot_matches, struct_matches, proteomes):
        store.close()

    logging.info(
        "{:>12} ({:.0f} proteins/sec)".format(
            total, cnt // (time.time() - ts)
        )
    )

    if not os.path.isdir(dst):
        os.makedirs(dst)

    chunk_dir = dst
    chunk = []          # entries to write in the current chunk
    dir_cnt = 0         # number of files in the current dir
    n_xref = 0          # number of cross-references in the current chunk
    total_xref = 0      # total number of cross-references
    logging.info("dumping entries to JSON files")
    for accession in sorted(entries_info):
        entry = entries_info[accession]

        fields = {
            "id": accession.upper(),
            "short_name": entry["short_name"],
            "name": entry["name"],
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
            chunk_dir, dir_cnt = _write_json(name, version, release_date,
                                             chunk, chunk_dir, dir_size,
                                             dir_cnt)

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
        chunk_dir, dir_cnt = _write_json(name, version, release_date,
                                         chunk, chunk_dir, dir_size,
                                         dir_cnt)

    logging.info(
        "{} cross-references for {} entries".format(
            total_xref, len(entries_info)
        )
    )


def _write_json(name, version, rel_date, entries, dst, dir_size, dir_cnt):
    # number of characters in files/dirs name
    n_chars = math.ceil(math.log10(dir_size)) + 1

    if dir_cnt == dir_size:
        # Too many files in current directory

        dst = os.path.join(dst, "{:0{}d}".format(dir_cnt + 1, n_chars))
        dir_cnt = 0
        os.mkdir(dst)

    filepath = os.path.join(dst, "{:0{}d}.json".format(dir_cnt + 1, n_chars))
    with open(filepath, "wt") as fh:
        json.dump({
            "name": name,
            "release": version,
            "release_date": rel_date,
            "entry_count": len(entries),
            "entries": entries
        }, fh, indent=4)

    return dst, dir_cnt + 1
