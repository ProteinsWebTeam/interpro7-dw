import MySQLdb

from interpro7dw import pfam
from interpro7dw.utils import logger
from interpro7dw.utils.mysql import url2dict
from interpro7dw.utils.store import loadobj
from .utils import jsonify


def insert_entries(ipr_url: str, pfam_url: str, entries_file: str):
    logger.info("fetching Wikipedia data for Pfam entries")
    wiki = pfam.get_wiki(pfam_url)

    logger.info("loading Pfam curation/family details")
    pfam_details = pfam.get_details(pfam_url)

    logger.info("populating webfront_entry")
    entries = loadobj(entries_file)

    con = MySQLdb.connect(**url2dict(ipr_url), charset="utf8mb4")
    cur = con.cursor()
    cur.execute("DROP TABLE IF EXISTS webfront_entry")
    cur.execute(
        """
        CREATE TABLE webfront_entry
        (
            entry_id VARCHAR(10) DEFAULT NULL,
            accession VARCHAR(25) PRIMARY KEY NOT NULL,
            type VARCHAR(50) NOT NULL,
            name LONGTEXT,
            short_name VARCHAR(100),
            source_database VARCHAR(10) NOT NULL,
            member_databases LONGTEXT,
            integrated_id VARCHAR(25),
            go_terms LONGTEXT,
            description LONGTEXT,
            wikipedia LONGTEXT,
            details LONGTEXT,
            literature LONGTEXT,
            hierarchy LONGTEXT,
            cross_references LONGTEXT,
            interactions LONGTEXT,
            pathways LONGTEXT,
            overlaps_with LONGTEXT,
            taxa LONGTEXT NOT NULL,
            is_featured TINYINT NOT NULL,
            is_alive TINYINT NOT NULL,
            history LONGTEXT,
            entry_date DATETIME NOT NULL,
            deletion_date DATETIME,
            counts LONGTEXT NOT NULL
        ) CHARSET=utf8mb4 DEFAULT COLLATE=utf8mb4_unicode_ci
        """
    )

    query = """
        INSERT INTO webfront_entry
        VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s,
          %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
    """

    args = []
    for entry in entries.values():
        args.append((
            None,
            entry.accession,
            entry.type.lower(),
            entry.name,
            entry.shot_name,
            entry.database,
            jsonify(entry.integrates, nullable=True),
            entry.integrated_in,
            jsonify(entry.go_terms, nullable=True),
            jsonify(entry.descriptions, nullable=True),
            jsonify(wiki.get(entry.accession), nullable=True),
            jsonify(pfam_details.get(entry.accession), nullable=True),
            jsonify(entry.literature, nullable=True),
            jsonify(entry.hierarchy, nullable=True),
            jsonify(entry.xrefs, nullable=True),
            jsonify(entry.ppi, nullable=True),
            jsonify(entry.pathways, nullable=True),
            jsonify(entry.overlaps_with, nullable=True),
            jsonify(entry.taxa, nullable=False),
            0,
            1 if entry.is_public else 0,
            jsonify(entry.history, nullable=True),
            entry.creation_date,
            entry.deletion_date,
            jsonify(entry.counts, nullable=False)
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

    logger.info("complete")
