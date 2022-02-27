import pickle

import MySQLdb

from interpro7dw.utils import logger
from interpro7dw.utils.mysql import uri2dict
from interpro7dw.utils.store import BasicStore

from .utils import jsonify


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
            counts LONGTEXT NOT NULL
        ) CHARSET=utf8mb4 DEFAULT COLLATE=utf8mb4_unicode_ci
        """
    )

    query = """
        INSERT INTO webfront_proteome VALUES (%s, %s, %s, %s, %s, %s, %s) 
    """
    params = []

    with BasicStore(xrefs_file, mode="r") as store:
        for proteome_id, xrefs in store:
            proteome = proteomes[proteome_id]

            # Adds total number of entries
            num_entries = {"total": 0}
            for database, entries in xrefs["entries"].items():
                num_entries[database] = len(entries)
                num_entries["total"] += len(entries)

            params.append((
                proteome_id,
                proteome["name"],
                1 if proteome["is_reference"] else 0,
                proteome["strain"],
                proteome["assembly"],
                proteome["taxon_id"],
                jsonify({
                    "domain_architectures": len(xrefs["dom_orgs"]),
                    "entries": entries,
                    "proteins": xrefs["proteins"],
                    "sets": len(xrefs["sets"]),
                    "structures": len(xrefs["structures"]),
                    "taxa": len(xrefs["taxa"])
                })
            ))

            if len(params) == 1000:
                cur.executemany(query, params)
                params = []

    if params:
        cur.executemany(query, params)

    con.commit()
    cur.close()
    con.close()

    logger.info("done")
