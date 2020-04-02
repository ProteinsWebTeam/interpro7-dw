# -*- coding: utf-8 -*-

import json
from typing import Optional

import MySQLdb

from interpro7dw import logger
from interpro7dw.ebi.interpro import production as ippro
from interpro7dw.ebi.interpro.utils import Table
from interpro7dw.utils import DataDump, DirectoryTree
from interpro7dw.utils import datadump, dataload, deepupdate, url2dict
from .utils import jsonify, merge_xrefs, reduce


def init_clans(pro_url: str, stg_url: str, output: str, threshold: float=1e-2):
    logger.info("loading clans")
    clans = ippro.get_clans(pro_url)
    entry2clan = {}
    for accession, clan in clans.items():
        for entry_acc, score, seq_length in clan["members"]:
            entry2clan[entry_acc] = (accession, seq_length)

    logger.info("inserting profile-profile alignments")
    con = MySQLdb.connect(**url2dict(stg_url))
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
        ) CHARSET=utf8 DEFAULT COLLATE=utf8_unicode_ci
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

    datadump(output, clans)
    logger.info("complete")


def dump_xrefs(xrefs: dict, output: str):
    with DataDump(output) as f:
        for clan_acc in sorted(xrefs):
            f.dump((clan_acc, xrefs[clan_acc]))


def insert_clans(url: str, p_clans: str, p_entries: str, p_entry2xrefs: str,
                 dir: Optional[str]=None, max_xrefs: int=1000000):
    dt = DirectoryTree(dir)
    entry2clan = {}
    for entry_acc, entry in dataload(p_entries).items():
        if entry.clan:
            entry2clan[entry_acc] = entry.clan["accession"]

    clans = {}
    files = []
    num_xrefs = 0
    with DataDump(p_entry2xrefs, compress=True) as entry2xrefs:
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

    clans = dataload(p_clans)

    con = MySQLdb.connect(**url2dict(url))
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
            counts LONGTEXT DEFAULT NULL
        ) CHARSET=utf8 DEFAULT COLLATE=utf8_unicode_ci
        """
    )
    cur.close()

    sql = """
        INSERT INTO webfront_set
        VALUES (%s, %s, %s, %s, %s, %s)
    """

    with Table(con, sql) as table:
        for clan_acc, xrefs in merge_xrefs(files):
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
                jsonify(clan["relationships"]),
                jsonify(counts)
            ))

    con.commit()
    con.close()
    logger.info("complete")