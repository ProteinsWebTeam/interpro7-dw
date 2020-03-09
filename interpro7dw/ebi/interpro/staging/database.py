# -*- coding: utf-8 -*-

import json
from datetime import datetime
from typing import List

import MySQLdb

from interpro7dw import logger
from interpro7dw.ebi.interpro import production as ippro
from interpro7dw.ebi.interpro.utils import Table
from interpro7dw.utils import Store, dataload, url2dict


def get_entry_databases(url: str) -> List[str]:
    con = MySQLdb.connect(**url2dict(url))
    cur = con.cursor()
    cur.execute("SELECT name FROM webfront_database WHERE type='entry'")
    names = [row[0] for row in cur]
    cur.close()
    con.close()
    return names


def insert_databases(pro_url: str, stg_url: str, **kwargs):
    version = kwargs.get("version")
    date = kwargs.get("date")

    con = MySQLdb.connect(**url2dict(stg_url))
    cur = con.cursor()
    cur.execute("DROP TABLE IF EXISTS webfront_database")
    cur.execute(
        """
        CREATE TABLE webfront_database
        (
            name VARCHAR(10) NOT NULL PRIMARY KEY,
            name_long VARCHAR(25) NOT NULL,
            description LONGTEXT,
            type ENUM('protein', 'entry', 'other') NOT NULL,
            version VARCHAR(20),
            release_date DATETIME,
            prev_version VARCHAR(20),
            prev_release_date DATETIME
        ) CHARSET=utf8 DEFAULT COLLATE=utf8_unicode_ci
        """
    )
    cur.close()

    sql = """
        INSERT INTO webfront_database VALUES (%s, %s, %s, %s, %s, %s, %s, %s) 
    """
    with Table(con, sql) as table:
        for record in ippro.get_databases(pro_url):
            if record[0] == "interpro" and version and version != record[4]:
                # Oracle has not been updated yet for the upcoming release
                record = list(record)  # tuple -> list
                record[6] = record[4]  # Previous version = Oracle's current version
                record[7] = record[5]
                record[4] = version    # Current version
                record[5] = datetime.strptime(date, "%Y-%m-%d")
                record = tuple(record)

            table.insert(record)

    con.commit()
    con.close()


def make_release_notes(p_entries: str, p_proteins: str, p_proteomes: str,
                       p_structures: str, p_taxonomy: str,
                       p_uniprot2entries: str, p_uniprot2proteome: str,
                       rel_url: str, stg_url: str):
    logger.info("preparing data")
    uniprot2pdbe = {}
    for pdb_id, entry in dataload(p_structures).items():
        for uniprot_acc in entry["proteins"]:
            try:
                uniprot2pdbe[uniprot_acc].add(pdb_id)
            except KeyError:
                uniprot2pdbe[uniprot_acc] = {pdb_id}

    con = MySQLdb.connect(**url2dict(stg_url))
    cur = con.cursor()
    cur.execute(
        """
        SELECT name_long, version
        FROM webfront_database
        WHERE name_long IN ('UniProtKB', 'UniProtKB/Swiss-Prot', 'UniProtKB/TrEMBL')
        """
    )
    uniprot = {}
    for name, version in cur:
        uniprot[name] = {
            "version": version,
            "count": 0,
            "signatures": 0,
            "integrated_signatures": 0
        }
    cur.close()
    con.close()

    proteins = Store(p_proteins)
    u2entries = Store(p_uniprot2entries)
    u2proteome = Store(p_uniprot2proteome)

    # Entities found in InterPro
    integrated_proteomes = set()
    integrated_structures = set()
    integrated_taxonomy = set()

    i = 0
    for uniprot_acc, info in proteins.items():
        i += 1
        if not i % 10000000:
            logger.info(f"{i:>12,}")

        if info["reviewed"]:
            database = uniprot["UniProtKB/Swiss-Prot"]
        else:
            database = uniprot["UniProtKB/TrEMBL"]

        database["count"] += 1

        try:
            entries = u2entries[uniprot_acc]
        except KeyError:
            continue

        # Protein matched by at least one signature
        database["signatures"] += 1

        for entry_acc, database, clan, go_terms in entries:
            if database == "interpro":
                """
                Protein matched by at least one InterPro entry,
                i.e. at least one integrated signature
                """
                database["integrated_signatures"] += 1

                try:
                    proteome_id = u2proteome[uniprot_acc]
                except KeyError:
                    pass
                else:
                    integrated_proteomes.add(proteome_id)

                try:
                    pdb_ids = uniprot2pdbe[uniprot_acc]
                except KeyError:
                    pass
                else:
                    integrated_structures |= pdb_ids

                integrated_taxonomy.add(info["taxid"])

                break

    proteins.close()
    u2entries.close()
    u2proteome.close()

    logger.info(f"{i:>12,}")

    # Sum Swiss-Prot and TrEMBL counts
    for key in ["count", "signatures", "integrated_signatures"]:
        value_sp = uniprot["UniProtKB/Swiss-Prot"][key]
        value_tr = uniprot["UniProtKB/TrEMBL"][key]
        uniprot["UniProtKB"][key] = value_sp + value_tr

    logger.info("tracking changes since last releases")
    con = MySQLdb.connect(**url2dict(rel_url))
    cur = con.cursor()
    cur.execute(
        """
        SELECT accession, source_database, integrated_id
        FROM webfront_entry 
        WHERE is_alive = 1
        """
    )
    public_entries = set()
    public_integrated = set()
    for entry_acc, database, integrated_in in cur:
        if database == "interpro":
            public_entries.add(entry_acc)
        elif integrated_in:
            # Signature already integrated in the previous release
            public_integrated.add(entry_acc)

    cur.execute(
        """
        SELECT name, version 
        FROM webfront_database 
        WHERE type = 'entry'
        """
    )
    public_databases = dict(cur.fetchall())
    cur.execute("SELECT * FROM webfront_release_note")
    prev_releases = cur.fetchall()
    cur.close()
    con.close()

    con = MySQLdb.connect(**url2dict(stg_url))
    cur = con.cursor()
    cur.execute(
        """
        CREATE TABLE IF NOT EXISTS webfront_release_note
        (
            version VARCHAR(20) PRIMARY KEY NOT NULL,
            release_date DATETIME NOT NULL,
            content LONGTEXT NOT NULL
        ) CHARSET=utf8 DEFAULT COLLATE=utf8_unicode_ci
        """
    )
    cur.execute(
        """
        SELECT name, name_long, version, release_date
        FROM webfront_database 
        WHERE type = 'entry'
        """
    )
    staging_databases = {row[0]: (row[1], row[2], row[3]) for row in cur}
    cur.execute("SELECT COUNT(*) FROM webfront_release_note")
    if not cur.fetchone()[0]:
        # Table is empty: import previous release notes
        cur.executemany(
            """
            INSERT INTO webfront_release_note
            VALUES (%s, %s, %s)
            """, prev_releases
        )
        con.commit()
        prev_releases = None

    interpro_new = []
    interpro_types = {}
    member_databases = {}
    pubmed_citations = set()
    interpro2go = None
    latest_entry = None

    for entry in sorted(dataload(p_entries).values(), key=lambda e: e.creation_date):
        if entry.is_deleted:
            continue

        if entry.database == "interpro":
            for pub in entry.literature.values():
                if pub["PMID"] is not None:
                    pubmed_citations.add(pub["PMID"])

            try:
                interpro_types[entry.type] += 1
            except KeyError:
                interpro_types[entry.type] = 1

            if entry.accession not in public_entries:
                interpro_new.append(entry.accession)

            interpro2go += len(entry.go_terms)
            latest_entry = entry.accession
        else:
            try:
                obj = member_databases[entry.database]
            except KeyError:
                database, version, _ = staging_databases[entry.database]

                is_new = is_updated = False
                if entry.database not in public_databases:
                    is_new = True
                elif version != public_databases[entry.database]:
                    is_updated = True

                obj = member_databases[entry.database] = {
                    "name": database,
                    "version": version,
                    "signatures": 0,
                    "integrated_signatures": 0,
                    "recently_integrated": [],
                    "is_new": is_new,
                    "is_updated": is_updated,
                    "sets": set()
                }

            obj["signatures"] += 1
            if entry.integrated_in:
                obj["integrated_signatures"] += 1

                if entry.accession not in public_integrated:
                    # Recent integration
                    obj["recently_integrated"].append(entry.accession)

            if entry.clan:
                obj["sets"].add(entry.clan["accession"])

    # Transform sets of clans to counts:
    for obj in member_databases.values():
        obj["sets"] = len(obj["sets"])

    structures = list(dataload(p_structures).values())

    proteomes = set(dataload(p_proteomes).keys())
    errors = integrated_proteomes - proteomes
    if errors:
        raise RuntimeError(f"{len(errors)} invalid proteomes")

    taxa = set(dataload(p_taxonomy).keys())
    errors = integrated_taxonomy - taxa
    if errors:
        raise RuntimeError(f"{len(errors)} invalid taxa")

    content = {
        "interpro": {
            "entries": sum(interpro_types.values()),
            "new_entries": interpro_new,
            "latest_entry": latest_entry,
            "types": interpro_types,
            "go_terms": interpro2go
        },
        "member_databases": member_databases,
        "proteins": uniprot,
        "structures": {
            "total": len(structures),
            "integrated": len(integrated_structures),
            "version": max(entry["data"] for entry in structures)
        },
        "proteomes": {
            "total": len(proteomes),
            "integrated": len(integrated_proteomes),
            "version": uniprot["UniProtKB"]["version"]
        },
        "taxonomy": {
            "total": len(taxa),
            "integrated": len(integrated_taxonomy),
            "version": uniprot["UniProtKB"]["version"]
        },
        "citations": len(pubmed_citations)
    }

    _, version, date = staging_databases["interpro"]
    cur.execute(
        """
        SELECT COUNT(*)
        FROM webfront_release_note
        WHERE version = %s
        """, (version,)
    )
    n_rows, = cur.fetchone()

    if n_rows:
        cur.execute(
            """
            UPDATE webfront_release_note
            SET content = %s
            WHERE version = %s
            """, (json.dumps(content), version)
        )
    else:
        cur.execute(
            """
            INSERT INTO webfront_release_note
            VALUES (%s, %s, %s)
            """, (version, date, json.dumps(content) )
        )

    con.commit()
    cur.close()
    con.close()
    logger.info("complete")
