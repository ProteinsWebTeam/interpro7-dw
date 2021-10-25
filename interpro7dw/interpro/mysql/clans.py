import MySQLdb

from interpro7dw.utils import logger
from interpro7dw.utils.mysql import url2dict
from interpro7dw.utils.store import loadobj, SimpleStore
from .utils import jsonify


def insert_clans(url: str, clans_file: str, clanxrefs_file: str,
                 alignments_file: str):
    logger.info("loading clans")
    clans = loadobj(clans_file)

    logger.info("creating webfront_set")
    con = MySQLdb.connect(**url2dict(url), charset="utf8mb4")
    cur = con.cursor()
    cur.execute("DROP TABLE IF EXISTS webfront_set")
    cur.execute(
        """
        CREATE TABLE webfront_set
        (
            accession VARCHAR(20) PRIMARY KEY NOT NULL,
            name VARCHAR(400),
            description TEXT,
            source_database VARCHAR(10) NOT NULL,
            relationships LONGTEXT NOT NULL,
            authors TEXT,
            literature TEXT,
            counts LONGTEXT DEFAULT NULL
        ) CHARSET=utf8mb4 DEFAULT COLLATE=utf8mb4_unicode_ci
        """
    )

    query = """
        INSERT INTO webfront_set
        VALUES (%s, %s, %s, %s, %s, %s, %s, %s)
    """

    with SimpleStore(clanxrefs_file) as store:
        args = []

        for accession, xrefs in store:
            clan = clans[accession]
            record = (
                accession,
                clan["name"],
                clan["description"],
                clan["database"],
                jsonify(clan["relationships"], nullable=False),
                jsonify(clan.get("authors"), nullable=True),     # only Pfam
                jsonify(clan.get("literature"), nullable=True),  # only Pfam
                jsonify({
                    "domain_architectures": len(xrefs["dom_orgs"]),
                    "entries": {k: len(v) for k, v in xrefs["entries"].items()},
                    "proteins": len(xrefs["proteins"]),
                    "proteomes": len(xrefs["proteomes"]),
                    "structures": len(xrefs["structures"]),
                    "taxa": len(xrefs["taxa"])
                })
            )

            args.append(record)
            if len(args) == 1000:
                cur.executemany(query, args)
                args.clear()

        if args:
            cur.executemany(query, args)
            args.clear()

    logger.info("creating webfront_alignment")
    cur.execute("DROP TABLE IF EXISTS webfront_alignment")
    cur.execute(
        """
        CREATE TABLE webfront_alignment
        (
            id INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
            set_acc VARCHAR(20) NOT NULL,
            entry_acc VARCHAR(25) NOT NULL,
            target_acc VARCHAR(25) NOT NULL,
            target_set_acc VARCHAR(20),
            score DOUBLE NOT NULL,
            seq_length MEDIUMINT NOT NULL,
            domains TEXT NOT NULL
        ) CHARSET=utf8mb4 DEFAULT COLLATE=utf8mb4_unicode_ci
        """
    )

    query = """
        INSERT INTO webfront_alignment (
            set_acc, entry_acc, target_acc, target_set_acc, score,
            seq_length, domains
        )
        VALUES (%s, %s, %s, %s, %s, %s, %s)
    """

    with SimpleStore(alignments_file) as store:
        args = []

        for alignment in store:
            args.append(alignment)

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
