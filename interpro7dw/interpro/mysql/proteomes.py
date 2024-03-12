import pickle

import MySQLdb

from interpro7dw.utils import logger
from interpro7dw.utils.mysql import uri2dict
from interpro7dw.utils.store import BasicStore
from .utils import create_index, jsonify


def populate(uri: str, proteomes_file: str, xrefs_file: str):
    logger.info("loading proteomes")
    with open(proteomes_file, "rb") as fh:
        proteomes = pickle.load(fh)

    logger.info("creating webfront_proteome")
    con = MySQLdb.connect(**uri2dict(uri), charset="utf8mb4")
    cur = con.cursor()
    cur.execute("DROP TABLE IF EXISTS webfront_proteome")
    cur.execute(
        """
        CREATE TABLE webfront_proteome
        (
            accession VARCHAR(20) PRIMARY KEY NOT NULL,
            name VARCHAR(215) NOT NULL,
            is_reference TINYINT NOT NULL,
            strain VARCHAR(512),
            assembly VARCHAR(512),
            taxonomy_id VARCHAR(20) NOT NULL,
            num_proteins INT NOT NULL,
            counts LONGTEXT NOT NULL
        ) CHARSET=utf8mb4 DEFAULT COLLATE=utf8mb4_unicode_ci
        """
    )
    cur.execute("DROP TABLE IF EXISTS webfront_proteomeperentry")
    cur.execute(
        """
        CREATE TABLE webfront_proteomeperentry
        (
          id INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
          accession VARCHAR(20) NOT NULL,
          entry_acc VARCHAR(30) NOT NULL,
          num_proteins INT NOT NULL,
          counts LONGTEXT NULL NULL
        ) CHARSET=utf8mb4 DEFAULT COLLATE=utf8mb4_unicode_ci
        """
    )
    cur.execute("DROP TABLE IF EXISTS webfront_proteomeperentrydb")
    cur.execute(
        """
        CREATE TABLE webfront_proteomeperentrydb
        (
          id INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
          accession VARCHAR(20) NOT NULL,
          source_database VARCHAR(10) NOT NULL,
          num_proteins INT NOT NULL,
          counts LONGTEXT NOT NULL
        ) CHARSET=utf8mb4 DEFAULT COLLATE=utf8mb4_unicode_ci
        """
    )

    query1 = """
            INSERT INTO webfront_proteome
            VALUES (%s, %s, %s, %s, %s, %s, %s, %s)
        """
    params1 = []
    query2 = """
            INSERT INTO webfront_proteomeperentry
                (accession, entry_acc, num_proteins, counts)
            VALUES (%s, %s, %s, %s) 
        """
    params2 = []
    query3 = """
            INSERT INTO webfront_proteomeperentrydb 
                (accession, source_database, num_proteins, counts)
            VALUES (%s, %s, %s, %s) 
        """
    params3 = []

    with BasicStore(xrefs_file, mode="r") as store:
        for proteome_id, xrefs in store:
            proteome = proteomes[proteome_id]

            # Build dict of member databases
            databases = {}
            for database, obj in xrefs["proteins"]["databases"].items():
                db = databases[database.lower()] = {
                    "proteins": obj["count"],
                    "entries": {},
                    "structures": set()
                }

                for entry_acc, num_proteins in obj["entries"].items():
                    db["entries"][entry_acc] = {
                        "proteins": num_proteins,
                        "structures": set()
                    }

            # Add structures matched by entries
            structures = xrefs["structures"]["all"]
            for database, obj in xrefs["structures"]["databases"].items():
                try:
                    db = databases[database.lower()]
                except KeyError:
                    db = databases[database.lower()] = {
                        "proteins": 0,
                        "entries": {},
                        "structures": set()
                    }

                for entry_acc, entry_structures in obj.items():
                    try:
                        e = db["entries"][entry_acc]
                    except KeyError:
                        e = db["entries"][entry_acc] = {
                            "proteins": 0,
                            "structures": set()
                        }

                    e["structures"] = entry_structures
                    db["structures"] |= entry_structures
                    structures |= entry_structures

            # Track total number of entries across all databases
            entries_per_db = {"total": 0}
            for database, db in databases.items():
                entries_per_db[database] = 0
                for entry_acc, e in db["entries"].items():
                    entries_per_db["total"] += 1
                    entries_per_db[database] += 1
                    params2.append((
                        proteome_id,
                        entry_acc,
                        e["proteins"],
                        jsonify({
                            "proteomes": len(xrefs["proteomes"]),
                            "proteins": e["proteins"],
                            "structures": len(e["structures"])
                        })
                    ))

                    if len(params2) == 1000:
                        cur.executemany(query2, params2)
                        params2 = []

                params3.append((
                    proteome_id,
                    database,
                    db["proteins"],
                    jsonify({
                        "entries": entries_per_db[database],
                        "proteins": db["proteins"],
                        "structures": len(db["structures"])
                    })
                ))

                if len(params3) == 1000:
                    cur.executemany(query3, params3)
                    params3 = []

            params1.append((
                proteome_id,
                proteome["name"],
                1 if proteome["is_reference"] else 0,
                proteome["strain"],
                proteome["assembly"],
                proteome["taxon_id"],
                xrefs["proteins"]["all"],
                jsonify({
                    "entries": entries_per_db,
                    "proteins": xrefs["proteins"]["all"],
                    "structures": len(structures),
                })
            ))

            if len(params1) == 1000:
                cur.executemany(query1, params1)
                params1 = []

    for query, params in zip([query1, query2, query3],
                             [params1, params2, params3]):
        if params:
            cur.executemany(query, params)

    con.commit()
    cur.close()
    con.close()

    logger.info("done")


def index(uri: str):
    con = MySQLdb.connect(**uri2dict(uri), charset="utf8mb4")
    cur = con.cursor()
    logger.info("i_webfront_proteomeperentry_tax_entry")
    create_index(
        cur,
        """
        CREATE UNIQUE INDEX i_webfront_proteomeperentry_tax_entry 
        ON webfront_proteomeperentry (accession, entry_acc)
        """
    )
    logger.info("i_webfront_proteomeperentrydb_tax_db")
    create_index(
        cur,
        """
        CREATE INDEX i_webfront_proteomeperentrydb_tax_db
        ON webfront_proteomeperentrydb (accession, source_database)
        """
    )
    logger.info("i_webfront_proteomeperentrydb_tax")
    create_index(
        cur,
        """
        CREATE INDEX i_webfront_proteomeperentrydb_tax
        ON webfront_proteomeperentrydb (accession)
        """
    )
    logger.info("i_webfront_proteomeperentrydb_db")
    create_index(
        cur,
        """
        CREATE INDEX i_webfront_proteomeperentrydb_db
        ON webfront_proteomeperentrydb (source_database)
        """
    )

    cur.close()
    con.close()
    logger.info("done")
