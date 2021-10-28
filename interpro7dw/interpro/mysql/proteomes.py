import MySQLdb

from interpro7dw.utils import logger
from interpro7dw.utils.mysql import url2dict
from interpro7dw.utils.store import SimpleStore, loadobj

from .utils import jsonify


def insert_proteomes(url: str, proteomes_file: str, xrefs_file: str):
    """

    :param url:
    :param proteomes_file:
    :param xrefs_file:
    :return:
    """
    logger.info("creating webfront_proteome")
    con = MySQLdb.connect(**url2dict(url), charset="utf8mb4")
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
    args = []

    proteomes = loadobj(proteomes_file)
    with SimpleStore(xrefs_file) as store:
        for proteome_id, xrefs in store:
            proteome = proteomes[proteome_id]

            # Adds total number of entries
            entries = {}
            total = 0
            for db, accessions in xrefs["entries"].items():
                entries[db] = len(accessions)
                total += entries[db]

            entries["total"] = total

            args.append((
                proteome_id,
                proteome["name"],
                1 if proteome["is_reference"] else 0,
                proteome["strain"],
                proteome["assembly"],
                proteome["taxon_id"],
                jsonify({
                    "domain_architectures": len(xrefs["domain_architectures"]),
                    "entries": entries,
                    "proteins": xrefs["proteins"],
                    "sets": len(xrefs["sets"]),
                    "structures": len(xrefs["structures"]),
                    "taxa": len(xrefs["taxa"])
                })
            ))

            if len(args) == 1000:
                cur.executemany(query, args)
                args.clear()

    if args:
        cur.executemany(query, args)
        args.clear()

    con.commit()
    cur.close()
    con.close()

    logger.info("done")
