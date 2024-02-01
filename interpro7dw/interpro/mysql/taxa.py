import pickle

import MySQLdb

from interpro7dw.utils import logger
from interpro7dw.utils.store import BasicStore
from interpro7dw.utils.mysql import uri2dict
from .utils import create_index, jsonify


def populate(uri: str, taxa_file: str, xrefs_file: str):
    logger.info("loading taxa")
    with open(taxa_file, "rb") as fh:
        taxa = pickle.load(fh)

    logger.info("creating taxonomy tables")
    con = MySQLdb.connect(**uri2dict(uri), charset="utf8mb4")
    cur = con.cursor()
    cur.execute("DROP TABLE IF EXISTS webfront_taxonomy")
    cur.execute(
        """
        CREATE TABLE webfront_taxonomy
        (
            accession VARCHAR(20) PRIMARY KEY NOT NULL,
            scientific_name VARCHAR(255) NOT NULL,
            full_name VARCHAR(512) NOT NULL,
            lineage LONGTEXT NOT NULL,
            parent_id VARCHAR(20),
            rank VARCHAR(20) NOT NULL,
            children LONGTEXT,
            counts LONGTEXT NOT NULL
        ) CHARSET=utf8mb4 DEFAULT COLLATE=utf8mb4_unicode_ci
        """
    )
    cur.execute("DROP TABLE IF EXISTS webfront_taxonomyperentry")
    cur.execute(
        """
        CREATE TABLE webfront_taxonomyperentry
        (
          id INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
          tax_id VARCHAR(20) NOT NULL,
          entry_acc VARCHAR(30) NOT NULL,
          counts LONGTEXT NULL NULL
        ) CHARSET=utf8mb4 DEFAULT COLLATE=utf8mb4_unicode_ci
        """
    )
    cur.execute("DROP TABLE IF EXISTS webfront_taxonomyperentrydb")
    cur.execute(
        """
        CREATE TABLE webfront_taxonomyperentrydb
        (
          id INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
          tax_id VARCHAR(20) NOT NULL,
          source_database VARCHAR(10) NOT NULL,
          counts LONGTEXT NOT NULL
        ) CHARSET=utf8mb4 DEFAULT COLLATE=utf8mb4_unicode_ci
        """
    )

    query1 = """
        INSERT INTO webfront_taxonomy 
        VALUES (%s, %s, %s, %s, %s, %s, %s, %s)
    """
    params1 = []
    query2 = """
        INSERT INTO webfront_taxonomyperentry (tax_id,entry_acc,counts)
        VALUES (%s, %s, %s) 
    """
    params2 = []
    query3 = """
        INSERT INTO webfront_taxonomyperentrydb (tax_id,source_database,counts)
        VALUES (%s, %s, %s) 
    """
    params3 = []

    with BasicStore(xrefs_file, mode="r") as store:
        for taxon_id, xrefs in store:
            taxon = taxa[taxon_id]

            # Adds total number of entries
            num_entries = {"total": 0}

            for database, obj in xrefs["proteins"]["databases"].items():
                num_entries[database.lower()] = len(obj["entries"])
                num_entries["total"] += len(obj["entries"])

            params1.append((
                taxon_id,
                taxon["sci_name"],
                taxon["full_name"],
                f" {' '.join(taxon['lineage'])} ",
                taxon["parent"],
                taxon["rank"],
                jsonify(taxon["children"]),
                jsonify({
                    "entries": num_entries,
                    "proteomes": len(xrefs["proteomes"]),
                    "proteins": xrefs["proteins"]["all"],
                    "structures": len(xrefs["structures"]["all"]),
                })
            ))

            if len(params1) == 1000:
                cur.executemany(query1, params1)
                params1 = []

            for database, obj in xrefs["proteins"]["databases"].items():
                structures_in_db = set()
                for entry_acc, num_proteins in obj["entries"].items():
                    if entry_acc in xrefs["structures"]["entries"]:
                        structures = xrefs["structures"]["entries"][entry_acc]
                        num_structures = len(structures)
                        structures_in_db |= structures
                    else:
                        num_structures = 0

                    params2.append((
                        taxon_id,
                        entry_acc,
                        jsonify({
                            "proteomes": len(xrefs["proteomes"]),
                            "proteins": num_proteins,
                            "structures": num_structures
                        })
                    ))

                    if len(params2) == 1000:
                        cur.executemany(query2, params2)
                        params2 = []

                params3.append((
                    taxon_id,
                    database.lower(),
                    jsonify({
                        "entries": len(obj["entries"]),
                        "proteomes": len(xrefs["proteomes"]),
                        "proteins": obj["count"],
                        "structures": len(structures_in_db)
                    })
                ))

                if len(params3) == 1000:
                    cur.executemany(query3, params3)
                    params3 = []

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
    logger.info("i_webfront_taxonomyperentry_tax_entry")
    create_index(
        cur,
        """
        CREATE UNIQUE INDEX i_webfront_taxonomyperentry_tax_entry 
        ON webfront_taxonomyperentry (tax_id, entry_acc)
        """
    )
    logger.info("i_webfront_taxonomyperentrydb_tax_db")
    create_index(
        cur,
        """
        CREATE INDEX i_webfront_taxonomyperentrydb_tax_db
        ON webfront_taxonomyperentrydb (tax_id, source_database)
        """
    )
    logger.info("i_webfront_taxonomyperentrydb_tax")
    create_index(
        cur,
        """
        CREATE INDEX i_webfront_taxonomyperentrydb_tax
        ON webfront_taxonomyperentrydb (tax_id)
        """
    )
    logger.info("i_webfront_taxonomyperentrydb_db")
    create_index(
        cur,
        """
        CREATE INDEX i_webfront_taxonomyperentrydb_db
        ON webfront_taxonomyperentrydb (source_database)
        """
    )

    cur.close()
    con.close()
    logger.info("done")
