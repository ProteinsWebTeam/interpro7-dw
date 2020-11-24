# -*- coding: utf-8 -*-

import json
import os
from typing import Optional

import MySQLdb

from interpro7dw import logger
from interpro7dw.ebi.interpro import production as ippro
from interpro7dw.ebi.interpro.utils import Table
from interpro7dw.ebi.pfam import get_clans
from interpro7dw.utils import DumpFile, DirectoryTree
from interpro7dw.utils import dumpobj, loadobj, deepupdate, url2dict
from interpro7dw.utils import merge_dumps
from .utils import jsonify, reduce


def init_clans(pro_url: str, stg_url: str, output: str,
               threshold: float = 1e-2):
    logger.info("loading clans")
    clans = ippro.get_clans(pro_url)
    entry2clan = {}
    for accession, clan in clans.items():
        for entry_acc, score, seq_length in clan["members"]:
            entry2clan[entry_acc] = (accession, seq_length)

    logger.info("inserting profile-profile alignments")
    con = MySQLdb.connect(**url2dict(stg_url), charset="utf8mb4")
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
        INSERT INTO webfront_alignment (set_acc, entry_acc, target_acc,
                                        target_set_acc, score, seq_length,
                                        domains)
        VALUES (%s, %s, %s, %s, %s, %s, %s)
    """
    clan_links = {clan_acc: {} for clan_acc in clans}
    with Table(con, sql) as table:
        gen = ippro.iter_clan_alignments(pro_url)
        i = 0
        for query, target, score, domains in gen:
            i += 1
            if not i % 10000000:
                logger.info(f"{i:>12,}")

            try:
                clan_acc, seq_length = entry2clan[query]
            except KeyError:
                continue

            if score > threshold:
                continue

            try:
                target_clan_acc, _ = entry2clan[target]
            except KeyError:
                target_clan_acc = None

            table.insert((clan_acc, query, target, target_clan_acc,
                          score, seq_length, json.dumps(domains)))

            if clan_acc == target_clan_acc:
                # Query and target from the same clan: update the clan's links
                links = clan_links[clan_acc]

                if query > target:
                    query, target = target, query

                try:
                    targets = links[query]
                except KeyError:
                    links[query] = {target: score}
                else:
                    if target not in targets or score < targets[target]:
                        targets[target] = score

        logger.info(f"{i:>12,}")

    con.commit()
    cur = con.cursor()
    cur.execute("CREATE INDEX i_alignment "
                "ON webfront_alignment (set_acc)")
    cur.close()
    con.close()

    for clan_acc, clan in clans.items():
        nodes = []
        for accession, score, seq_length in clan["members"]:
            nodes.append({
                "accession": accession,
                "type": "entry",
                "score": score
            })

        links = []
        for query_acc, targets in clan_links[clan_acc].items():
            for target_acc, score in targets.items():
                links.append({
                    "source": query_acc,
                    "target": target_acc,
                    "score": score
                })

        clan["relationships"] = {
            "nodes": nodes,
            "links": links
        }

    dumpobj(output, clans)
    logger.info("complete")


def dump_xrefs(xrefs: dict, output: str):
    with DumpFile(output) as f:
        for clan_acc in sorted(xrefs):
            f.dump((clan_acc, xrefs[clan_acc]))


def insert_clans(pfam_url: str, stg_url: str, p_clans: str, p_entries: str,
                 p_entry2xrefs: str, dir: Optional[str] = None,
                 max_xrefs: int = 1000000):
    dt = DirectoryTree(dir)
    entry2clan = {}
    for entry_acc, entry in loadobj(p_entries).items():
        if entry.clan:
            entry2clan[entry_acc] = entry.clan["accession"]

    clans = {}
    files = []
    num_xrefs = 0
    with DumpFile(p_entry2xrefs) as entry2xrefs:
        for entry_acc, entry_xrefs in entry2xrefs:
            try:
                clan_acc = entry2clan[entry_acc]
            except KeyError:
                continue

            try:
                clan_xrefs = clans[clan_acc]
            except KeyError:
                clan_xrefs = clans[clan_acc] = {}

            cnt_before = sum(map(len, clan_xrefs.values()))
            deepupdate(entry_xrefs, clan_xrefs)
            cnt_after = sum(map(len, clan_xrefs.values()))
            num_xrefs += cnt_after - cnt_before

            if num_xrefs >= max_xrefs:
                output = dt.mktemp()
                dump_xrefs(clans, output)
                files.append(output)
                clans = {}
                num_xrefs = 0

    if clans:
        output = dt.mktemp()
        dump_xrefs(clans, output)
        files.append(output)

    logger.info(f"temporary files: "
                f"{sum(map(os.path.getsize, files))/1024/1024:.0f} MB")

    clans = loadobj(p_clans)

    logger.info("loading additional details for Pfam clans")
    pfam_clans = get_clans(pfam_url)

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

            try:
                pfam_clan = pfam_clans[clan_acc]
            except KeyError:
                pass
            else:
                """
                Replace `description`
                Add `authors` and `literature`
                """
                clan.update(pfam_clan)

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

    con.commit()
    con.close()
    dt.remove()
    logger.info("complete")
