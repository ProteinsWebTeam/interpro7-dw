import MySQLdb

from interpro7dw import pfam
from interpro7dw.utils import logger
from interpro7dw.utils.mysql import url2dict
from interpro7dw.utils.store import loadobj, SimpleStore
from .utils import jsonify


def insert_annotations(url: str, hmms_file: str, pfam_alignments: str):
    logger.info("creating webfront_entryannotation")
    con = MySQLdb.connect(**url2dict(url))
    cur = con.cursor()
    cur.execute("DROP TABLE IF EXISTS webfront_entryannotation")
    cur.execute(
        """
        CREATE TABLE webfront_entryannotation
        (
            annotation_id INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
            accession VARCHAR(25) NOT NULL,
            type VARCHAR(20) NOT NULL,
            value LONGBLOB NOT NULL,
            mime_type VARCHAR(32) NOT NULL,
            num_sequences INT
        ) CHARSET=utf8mb4 DEFAULT COLLATE=utf8mb4_unicode_ci
        """
    )

    for file in (hmms_file, pfam_alignments):
        with SimpleStore(file) as store:
            for accession, anno_type, anno_value, count in store:
                if anno_type == "logo":
                    mime_type = "application/json"
                else:
                    mime_type = "application/gzip"

                cur.execute(
                    """
                    INSERT INTO webfront_entryannotation (
                        accession, type, value, mime_type, num_sequences
                    ) VALUES (%s, %s, %s, %s, %s)
                    """, (accession, anno_type, anno_value, mime_type, count)
                )

    con.commit()

    logger.info("indexing")
    cur.execute(
        """
        CREATE INDEX i_entryannotation 
        ON webfront_entryannotation (accession)
        """
    )

    cur.close()
    con.close()

    logger.info("done")


def insert_databases(url: str, databases_file: str):
    logger.info("creating webfront_database")
    con = MySQLdb.connect(**url2dict(url), charset="utf8mb4")
    cur = con.cursor()
    cur.execute("DROP TABLE IF EXISTS webfront_database")
    cur.execute(
        """
        CREATE TABLE webfront_database
        (
            name VARCHAR(10) NOT NULL PRIMARY KEY,
            name_alt VARCHAR(10) NOT NULL,
            name_long VARCHAR(25) NOT NULL,
            description LONGTEXT,
            type ENUM('protein', 'entry', 'feature', 'other') NOT NULL,
            num_entries INTEGER,
            version VARCHAR(20),
            release_date DATETIME,
            prev_version VARCHAR(20),
            prev_release_date DATETIME
        ) CHARSET=utf8mb4 DEFAULT COLLATE=utf8mb4_unicode_ci
        """
    )

    query = """
        INSERT INTO webfront_database 
        VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s) 
    """
    args = []

    for database in loadobj(databases_file):
        args.append(database)

    cur.executemany(query, args)
    con.commit()
    cur.close()
    con.close()

    logger.info("done")


def insert_entries(ipr_url: str, pfam_url: str, entries_file: str,
                   entry2xrefs_file: str):
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

    with SimpleStore(entry2xrefs_file) as store:
        for accession, xrefs in store:
            entry = entries[accession]
            rec = (
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
                jsonify(xrefs["taxa"]["tree"], nullable=False),
                0,
                1 if entry.is_public else 0,
                jsonify(entry.history, nullable=True),
                entry.creation_date,
                entry.deletion_date,
                jsonify(entry.counts, nullable=False)
            )

            cur.execute(query, rec)

    con.commit()
    cur.close()
    con.close()

    logger.info("done")
