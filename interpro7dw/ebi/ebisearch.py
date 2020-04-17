# -*- coding: utf-8 -*-

import json
import os
import shutil
from typing import Tuple

import MySQLdb

from interpro7dw import logger
from interpro7dw.utils import DataDump, DirectoryTree, loadobj, url2dict


def _init_fields(entry) -> Tuple[list, list]:
    fields = [
        {
            "name": "id",
            "value": entry.accession
        },
        {
            "name": "name",
            "value": entry.name or entry.accession
        },
        {
            "name": "short_name",
            "value": entry.short_name or entry.accession
        },
        {
            "name": "type",
            "value": entry.type
        },
        {
            "name": "creation_date",
            "value": entry.creation_date.strftime("%Y-%m-%d")
        },
        {
            "name": "description",
            "value": " ".join(entry.description)
        }
    ]
    xrefs = []

    if entry.clan:
        # Add clan
        fields.append({
            "name": "set",
            "value": entry.clan["accession"]
        })

    if entry.database == "interpro":
        # InterPro entry

        for database, signatures in entry.integrates.items():
            # Add contributing database to fields
            fields.append({
                "name": "contributing_database",
                "value": database
            })

            for signature_acc in signatures:
                xrefs.append({
                    "dbname": database,
                    "dbkey": signature_acc
                })

        for database, references in entry.cross_references.items():
            for ref_id in references:
                xrefs.append({
                    "dbname": database,
                    "dbkey": ref_id
                })

        for citation in entry.literature.values():
            try:
                pmid = citation["PMID"]
            except KeyError:
                continue
            else:
                xrefs.append({
                    "dbname": "PUBMED",
                    "dbkey": pmid
                })

        for term in entry.go_terms:
            xrefs.append({
                "dbname": "GO",
                "dbkey": term["identifier"]
            })

        for rel_acc in entry.relations:
            xrefs.append({
                "dbname": "INTERPRO",
                "dbkey": rel_acc
            })
    else:
        # Member database signature
        if entry.integrated_in:
            xrefs.append({
                "dbname": "INTERPRO",
                "dbkey": entry.integrated_in
            })

        for citation in entry.literature.values():
            try:
                pmid = citation["PMID"]
            except KeyError:
                continue
            else:
                xrefs.append({
                    "dbname": "PUBMED",
                    "dbkey": pmid
                })

    return fields, xrefs


def export(url: str, p_entries: str, p_entry2xrefs: str, outdir: str,
           max_xrefs: int=100000):
    logger.info("preparing data")
    con = MySQLdb.connect(**url2dict(url))
    cur = con.cursor()
    cur.execute(
        """
        SELECT name, name_long, version, release_date
        FROM webfront_database
        WHERE type = 'entry'
        """
    )
    databases = {}
    release_version = release_date = None
    for name, full_name, version, date in cur:
        databases[name] = full_name

        if name == "interpro":
            release_version = version
            release_date = date.strftime("%Y-%m-%d")

    cur.execute(
        """
        SELECT accession, scientific_name
        FROM webfront_taxonomy
        """
    )
    sci_names = dict(cur.fetchall())
    cur.close()
    con.close()

    if release_version is None:
        raise RuntimeError("missing release version/date for InterPro")

    try:
        shutil.rmtree(outdir)
    except FileNotFoundError:
        pass
    finally:
        os.makedirs(outdir, mode=0o775)

    entries = loadobj(p_entries)

    logger.info("starting")
    i = 0
    types = {}
    num_xrefs = {}
    with DataDump(p_entry2xrefs, compress=True) as entry2xrefs:
        for accession, entry_xrefs in entry2xrefs:
            entry = entries.pop(accession)

            fields, xrefs = _init_fields(entry)

            fields.append({
                "name": "source_database",
                "value": databases[entry.database]
            })


            for uniprot_acc, uniprot_id in entry_xrefs["proteins"]:
                xrefs.append({
                    "dbname": "UNIPROT",
                    "dbkey": uniprot_acc
                })

                xrefs.append({
                    "dbname": "UNIPROT",
                    "dbkey": uniprot_id
                })

            for tax_id in entry_xrefs["taxa"]:
                xrefs.append({
                    "dbname": "TAXONOMY",
                    "dbkey": tax_id
                })

                xrefs.append({
                    "dbname": "TAXONOMY",
                    "dbkey": sci_names[tax_id]
                })

            for upid in entry_xrefs["proteomes"]:
                xrefs.append({
                    "dbname": "PROTEOMES",
                    "dbkey": upid
                })

            for pdbe_id in entry_xrefs["structures"]:
                xrefs.append({
                    "dbname": "PDB",
                    "dbkey": pdbe_id
                })

            entry_type = entry.type.lower()
            try:
                dt, items = types[entry_type]
            except KeyError:
                dt = DirectoryTree(outdir, entry_type)
                items = []
                types[entry_type] = (dt, items)
                num_xrefs[entry_type] = 0

            items.append({
                "fields": fields,
                "cross_references": xrefs
            })
            num_xrefs[entry_type] += len(xrefs)

            if num_xrefs[entry_type] >= max_xrefs:
                path = dt.mktemp(suffix=".json")
                with open(path, "wt") as fh:
                    json.dump({
                        "name": "InterPro",
                        "release": release_version,
                        "release_date": release_date,
                        "entry_count": len(items),
                        "entries": items
                    }, fh, indent=4)

                items.clear()
                num_xrefs[entry_type] = 0

            i += 1
            if not i % 10000:
                logger.info(f"{i:>12,}")

    # Export entries not matching any protein
    for entry in entries.values():
        if entry.is_deleted:
            continue

        fields, xrefs = _init_fields(entry)

        fields.append({
            "name": "source_database",
            "value": databases[entry.database]
        })

        entry_type = entry.type.lower()
        try:
            dt, items = types[entry_type]
        except KeyError:
            dt = DirectoryTree(outdir, entry_type)
            items = []
            types[entry_type] = (dt, items)
            num_xrefs[entry_type] = 0

        items.append({
            "fields": fields,
            "cross_references": xrefs
        })

        if num_xrefs[entry_type] >= max_xrefs:
            path = dt.mktemp(suffix=".json")
            with open(path, "wt") as fh:
                json.dump({
                    "name": "InterPro",
                    "release": release_version,
                    "release_date": release_date,
                    "entry_count": len(items),
                    "entries": items
                }, fh, indent=4)

            items.clear()
            num_xrefs[entry_type] = 0

        i += 1
        if not i % 10000:
            logger.info(f"{i:>12,}")

    logger.info(f"{i:>12,}")

    for entry_type, (dt, items) in types.items():
        if num_xrefs[entry_type]:
            path = dt.mktemp(suffix=".json")
            with open(path, "wt") as fh:
                json.dump({
                    "name": "InterPro",
                    "release": release_version,
                    "release_date": release_date,
                    "entry_count": len(items),
                    "entries": items
                }, fh, indent=4)

    logger.info("complete")


def publish(src: str, dst: str):
    try:
        shutil.rmtree(dst)
    except FileNotFoundError:
        shutil.move(src, dst)
