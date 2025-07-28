import pickle

from interpro7dw.utils import logger
from interpro7dw.utils.store import BasicStore
from .utils import connect, jsonify


def populate(uri: str, taxa_file: str, xrefs_file: str):
    logger.info("loading taxa")
    with open(taxa_file, "rb") as fh:
        taxa = pickle.load(fh)

    logger.info("creating taxonomy tables")
    con = connect(uri)
    cur = con.cursor()
    cur.execute("DROP TABLE IF EXISTS interpro.webfront_taxonomy")
    cur.execute(
        """
        CREATE TABLE interpro.webfront_taxonomy
        (
            accession VARCHAR(20) PRIMARY KEY NOT NULL,
            scientific_name VARCHAR(255) NOT NULL,
            full_name VARCHAR(512) NOT NULL,
            lineage TEXT NOT NULL,
            parent_id VARCHAR(20),
            "rank" VARCHAR(20) NOT NULL,
            children JSONB,
            num_proteins INTEGER NOT NULL,
            counts JSONB NOT NULL
        )
        """
    )
    cur.execute("DROP TABLE IF EXISTS interpro.webfront_taxonomyperentry")
    cur.execute(
        """
        CREATE TABLE interpro.webfront_taxonomyperentry
        (
          id SERIAL NOT NULL PRIMARY KEY,
          tax_id VARCHAR(20) NOT NULL,
          entry_acc VARCHAR(30) NOT NULL,
          num_proteins INTEGER NOT NULL,
          counts JSONB NULL NULL
        )
        """
    )
    cur.execute("DROP TABLE IF EXISTS interpro.webfront_taxonomyperentrydb")
    cur.execute(
        """
        CREATE TABLE interpro.webfront_taxonomyperentrydb
        (
          id SERIAL NOT NULL PRIMARY KEY,
          tax_id VARCHAR(20) NOT NULL,
          source_database VARCHAR(10) NOT NULL,
          num_proteins INTEGER NOT NULL,
          counts JSONB NOT NULL
        )
        """
    )

    query1 = """
        INSERT INTO interpro.webfront_taxonomy
        VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s)
    """
    params1 = []
    query2 = """
        INSERT INTO interpro.webfront_taxonomyperentry
            (tax_id, entry_acc, num_proteins, counts)
        VALUES (%s, %s, %s, %s) 
    """
    params2 = []
    query3 = """
        INSERT INTO interpro.webfront_taxonomyperentrydb 
            (tax_id, source_database, num_proteins, counts)
        VALUES (%s, %s, %s, %s) 
    """
    params3 = []

    with BasicStore(xrefs_file, mode="r") as store:
        for taxon_id, xrefs in store:
            taxon = taxa[taxon_id]

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
                        taxon_id,
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
                    taxon_id,
                    database,
                    db["proteins"],
                    jsonify({
                        "entries": entries_per_db[database],
                        "proteomes": len(xrefs["proteomes"]),
                        "proteins": db["proteins"],
                        "structures": len(db["structures"])
                    })
                ))

                if len(params3) == 1000:
                    cur.executemany(query3, params3)
                    params3 = []

            params1.append((
                taxon_id,
                taxon["sci_name"],
                taxon["full_name"],
                f" {' '.join(taxon['lineage'])} ",
                taxon["parent"],
                taxon["rank"],
                jsonify(taxon["children"]),
                xrefs["proteins"]["all"],
                jsonify({
                    "entries": entries_per_db,
                    "proteomes": len(xrefs["proteomes"]),
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
    con = connect(uri)
    cur = con.cursor()
    logger.info("i_webfront_taxonomyperentry_tax_entry")
    cur.execute(
        """
        CREATE UNIQUE INDEX IF NOT EXISTS i_interpro.webfront_taxonomyperentry_tax_entry 
        ON interpro.webfront_taxonomyperentry (tax_id, entry_acc)
        """
    )
    logger.info("i_webfront_taxonomyperentry_entry")
    cur.execute(
        """
        CREATE INDEX IF NOT EXISTS i_interpro.webfront_taxonomyperentry_entry 
        ON interpro.webfront_taxonomyperentry (entry_acc)
        """
    )
    logger.info("i_webfront_taxonomyperentrydb_tax_db")
    cur.execute(
        """
        CREATE INDEX IF NOT EXISTS i_interpro.webfront_taxonomyperentrydb_tax_db
        ON interpro.webfront_taxonomyperentrydb (tax_id, source_database)
        """
    )
    logger.info("i_webfront_taxonomyperentrydb_tax")
    cur.execute(
        """
        CREATE INDEX IF NOT EXISTS i_interpro.webfront_taxonomyperentrydb_tax
        ON interpro.webfront_taxonomyperentrydb (tax_id)
        """
    )
    logger.info("i_webfront_taxonomyperentrydb_db")
    cur.execute(
        """
        CREATE INDEX IF NOT EXISTS i_interpro.webfront_taxonomyperentrydb_db
        ON interpro.webfront_taxonomyperentrydb (source_database)
        """
    )
    con.commit()
    cur.close()
    con.close()
    logger.info("done")
