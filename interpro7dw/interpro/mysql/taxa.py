import MySQLdb

from interpro7dw.utils import logger
from interpro7dw.utils.store import loadobj, SimpleStore
from interpro7dw.utils.mysql import url2dict
from .utils import jsonify


def insert_taxa(url: str, entries_file: str, taxa_file: str, xrefs_file: str):
    logger.info("creating taxonomy tables")
    con = MySQLdb.connect(**url2dict(url), charset="utf8mb4")
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
          entry_acc VARCHAR(25) NOT NULL,
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
    args1 = []
    query2 = """
        INSERT INTO webfront_taxonomyperentry (tax_id,entry_acc,counts)
        VALUES (%s, %s, %s) 
    """
    args2 = []
    query3 = """
        INSERT INTO webfront_taxonomyperentrydb (tax_id,source_database,counts)
        VALUES (%s, %s, %s) 
    """
    args3 = []

    entry2db = {e.accession: e.database
                for e in loadobj(entries_file).values()}
    taxa = loadobj(taxa_file)
    with SimpleStore(xrefs_file) as store:
        for taxon_id, xrefs in store:
            taxon = taxa[taxon_id]

            # Adds total number of entries
            entries = {"total": 0}
            for db, accessions in xrefs["databases"].items():
                entries[db] = len(accessions)
                entries["total"] += entries[db]

            args1.append((
                taxon_id,
                taxon["sci_name"],
                taxon["full_name"],
                f" {' '.join(taxon['lineage'])} ",
                taxon["parent"],
                taxon["rank"],
                jsonify(taxon["children"]),
                jsonify({
                    "entries": entries,
                    "proteomes": len(xrefs["proteomes"]),
                    "proteins": xrefs["proteins"]["all"],
                    "structures": len(xrefs["structures"]["all"]),
                })
            ))

            if len(args1) == 1000:
                cur.executemany(query1, args1)
                args1.clear()

            struct_per_db = {}
            for acc, num_proteins in xrefs["proteins"]["entries"].items():
                try:
                    structures = xrefs["structures"]["entries"][acc]
                except KeyError:
                    num_structures = 0
                else:
                    num_structures = len(structures)

                    db = entry2db[acc]
                    try:
                        struct_per_db[db] |= structures
                    except KeyError:
                        struct_per_db[db] = structures.copy()

                args2.append((
                    taxon_id,
                    acc,
                    jsonify({
                        "proteomes": len(xrefs["proteomes"]),
                        "proteins": num_proteins,
                        "structures": num_structures
                    })
                ))

            if len(args2) >= 1000:
                cur.executemany(query2, args2)
                args2.clear()

            for db, num_proteins in xrefs["proteins"]["databases"].items():
                args3.append((
                    taxon_id,
                    db,
                    jsonify({
                        "entries": entries[db],
                        "proteomes": len(xrefs["proteomes"]),
                        "proteins": num_proteins,
                        "structures": len(struct_per_db.get(db, []))
                    })
                ))

            if len(args3) >= 1000:
                cur.executemany(query3, args3)
                args3.clear()

    if args1:
        cur.executemany(query1, args1)
        args1.clear()

    if args2:
        cur.executemany(query2, args2)
        args2.clear()

    if args3:
        cur.executemany(query3, args3)
        args3.clear()

    con.commit()

    logger.info("creating indexes")
    cur.execute(
        """
        CREATE INDEX i_webfront_taxonomyperentry_tax
        ON webfront_taxonomyperentry (tax_id)
        """
    )
    cur.execute(
        """
        CREATE INDEX i_webfront_taxonomyperentry_entry
        ON webfront_taxonomyperentry (entry_acc)
        """
    )
    cur.execute(
        """
        CREATE INDEX i_webfront_taxonomyperentrydb_tax
        ON webfront_taxonomyperentrydb (tax_id)
        """
    )
    cur.execute(
        """
        CREATE INDEX i_webfront_taxonomyperentrydb_database
        ON webfront_taxonomyperentrydb (source_database)
        """
    )

    cur.close()
    con.close()
    logger.info("done")
