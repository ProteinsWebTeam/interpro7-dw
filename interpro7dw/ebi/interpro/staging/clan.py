# -*- coding: utf-8 -*-

import MySQLdb

from interpro7dw import logger
from interpro7dw.ebi.interpro.utils import Table
from interpro7dw.utils import DumpFile, DirectoryTree
from interpro7dw.utils import loadobj, deepupdate, url2dict
from interpro7dw.utils import merge_dumps
from .utils import jsonify, reduce


def insert_clans(stg_url: str, p_alignments: str, p_clans: str, p_entries: str,
                 p_entry2xrefs: str, **kwargs):
    max_xrefs = kwargs.get("max_xrefs", 1000000)
    tmpdir = kwargs.get("tmpdir")

    logger.info("aggregating clan cross-references")
    dt = DirectoryTree(tmpdir)
    entry2clan = {}
    for entry_acc, entry in loadobj(p_entries).items():
        if entry.clan:
            entry2clan[entry_acc] = entry.clan["accession"]

    clans = {}
    files = []
    num_xrefs = 0
    with DumpFile(p_entry2xrefs) as df:
        for entry_acc, entry_xrefs in df:
            try:
                clan_acc = entry2clan[entry_acc]
            except KeyError:
                continue

            try:
                clan_xrefs = clans[clan_acc]
            except KeyError:
                clan_xrefs = clans[clan_acc] = {}

            # We do not need the number of matches
            del entry_xrefs["matches"]

            cnt_before = sum(map(len, clan_xrefs.values()))
            deepupdate(entry_xrefs, clan_xrefs)
            cnt_after = sum(map(len, clan_xrefs.values()))
            num_xrefs += cnt_after - cnt_before

            if num_xrefs >= max_xrefs:
                file = dt.mktemp()
                with DumpFile(file, compress=True) as df2:
                    for clan_acc in sorted(clans):
                        df2.dump((clan_acc, clans[clan_acc]))

                files.append(file)
                clans = {}
                num_xrefs = 0

    file = dt.mktemp()
    with DumpFile(file, compress=True) as df2:
        for clan_acc in sorted(clans):
            df2.dump((clan_acc, clans[clan_acc]))

    files.append(file)

    logger.info("inserting clans")
    clans = loadobj(p_clans)
    con = MySQLdb.connect(**url2dict(stg_url), charset="utf8mb4")
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
    cur.close()

    sql = """
        INSERT INTO webfront_set
        VALUES (%s, %s, %s, %s, %s, %s, %s, %s)
    """
    with Table(con, sql) as table:
        for clan_acc, xrefs in merge_dumps(files):
            clan = clans[clan_acc]
            counts = reduce(xrefs)
            counts["entries"] = {
                clan["database"]: len(clan["members"]),
                "total": len(clan["members"])
            }

            table.insert((
                clan_acc,
                clan["name"],
                clan["description"],
                clan["database"],
                jsonify(clan["relationships"], nullable=False),
                jsonify(clan.get("authors")),
                jsonify(clan.get("literature")),
                jsonify(counts)
            ))

    logger.info(f"temporary files: {dt.size / 1024 / 1024:.0f} MB")
    dt.remove()

    logger.info("inserting alignments")
    cur = con.cursor()
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
    cur.close()

    sql = """
        INSERT INTO webfront_alignment (
            set_acc, entry_acc, target_acc, target_set_acc, score, 
            seq_length, domains
        )
        VALUES (%s, %s, %s, %s, %s, %s, %s)
    """
    with DumpFile(p_alignments) as df, Table(con, sql) as table:
        for alignments in df:
            for aln in alignments:
                table.insert(aln)

    con.commit()
    con.close()

    logger.info("complete")
