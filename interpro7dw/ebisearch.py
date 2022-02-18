import json
import os
import shutil
from typing import Tuple
from xml.sax.saxutils import escape

import MySQLdb

from interpro7dw.utils import logger
from interpro7dw.utils.mysql import url2dict
from interpro7dw.utils.store import Directory, SimpleStore, loadobj


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
            "value": entry.type.lower()
        },
        {
            "name": "creation_date",
            "value": entry.creation_date.strftime("%Y-%m-%d")
        },
        {
            "name": "description",
            "value": escape(' '.join(entry.descriptions))
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

        parent, children = entry.relations
        if parent:
            children.insert(0, parent)

        for rel_acc in children:
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


def export(uri: str, entries_file: str, entry2xrefs_file: str, taxa_file: str,
           outdir: str, max_xrefs: int = 100000):
    """Creates JSON files containing entries (InterPro + signatures) and
    cross-references to be ingested by EBISearch

    :param uri: InterPro MySQL connection string
    :param entries_file: File of InterPro entries
    and member DB signatures
    :param entry2xrefs_file: File of entries cross-references
    :param taxa_file: File of taxonomic information
    :param outdir: Output directory
    :param max_xrefs: Maximum number of cross-references in a JSON file
    """
    logger.info("loading database versions")
    con = MySQLdb.connect(**url2dict(uri))
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

    cur.close()
    con.close()

    if release_version is None:
        raise RuntimeError("missing release version/date for InterPro")

    logger.info("loading taxonomic info")
    sci_names = {}
    for taxon_id, taxon in loadobj(taxa_file).items():
        sci_names[taxon_id] = taxon["sci_name"]

    try:
        shutil.rmtree(outdir)
    except FileNotFoundError:
        pass

    entries = loadobj(entries_file)

    logger.info("starting")
    i = 0
    types = {}
    num_xrefs = {}
    with SimpleStore(entry2xrefs_file) as store:
        for accession, entry_xrefs in store:
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

            for tax_id in entry_xrefs["taxa"]["all"]:
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
                directory, items = types[entry_type]
            except KeyError:
                type_dir = os.path.join(outdir, entry_type)
                os.makedirs(type_dir, exist_ok=True)
                directory = Directory(type_dir)
                items = []
                types[entry_type] = (directory, items)
                num_xrefs[entry_type] = 0

            items.append({
                "fields": fields,
                "cross_references": xrefs
            })
            num_xrefs[entry_type] += len(xrefs)

            if num_xrefs[entry_type] >= max_xrefs:
                path = directory.mktemp(suffix=".json")
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

    for entry in entries.values():
        if not entry.is_public:
            continue

        fields, xrefs = _init_fields(entry)
        fields.append({
            "name": "source_database",
            "value": databases[entry.database]
        })

        entry_type = entry.type.lower()
        try:
            directory, items = types[entry_type]
        except KeyError:
            type_dir = os.path.join(outdir, entry_type)
            os.makedirs(type_dir, exist_ok=True)
            directory = Directory(type_dir)
            items = []
            types[entry_type] = (directory, items)
            num_xrefs[entry_type] = 0

        items.append({
            "fields": fields,
            "cross_references": xrefs
        })

        i += 1

    logger.info(f"{i:>12,}")

    for entry_type, (directory, items) in types.items():
        if num_xrefs[entry_type]:
            path = directory.mktemp(suffix=".json")
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
        pass
    finally:
        shutil.copytree(src, dst)
