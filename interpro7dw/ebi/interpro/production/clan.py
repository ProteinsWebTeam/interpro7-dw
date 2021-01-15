# -*- coding: utf-8 -*-

import json
from typing import Dict

import cx_Oracle

from interpro7dw import logger
from interpro7dw.ebi import pfam
from interpro7dw.utils import dumpobj, DumpFile


def get_clans(cur: cx_Oracle.Cursor) -> Dict[str, dict]:
    cur.execute(
        """
        SELECT
          C.CLAN_AC, C.NAME, C.DESCRIPTION, LOWER(D.DBSHORT), M.MEMBER_AC, 
          M.LEN, M.SCORE
        FROM INTERPRO.CLAN C
        INNER JOIN INTERPRO.CV_DATABASE D
          ON C.DBCODE = D.DBCODE
        INNER JOIN INTERPRO.CLAN_MEMBER M
          ON C.CLAN_AC = M.CLAN_AC
        """
    )

    clans = {}
    for row in cur:
        accession = row[0]
        name = row[1]
        descr = row[2]
        database = row[3]
        member_acc = row[4]
        seq_length = row[5]
        score = row[6]

        try:
            c = clans[accession]
        except KeyError:
            c = clans[accession] = {
                "accession": accession,
                "name": name,
                "description": descr,
                "database": database,
                "members": []
            }
        finally:
            c["members"].append((member_acc, score, seq_length))

    return clans


def iter_alignments(cur: cx_Oracle.Cursor):
    cur.execute(
        """
        SELECT QUERY_AC, TARGET_AC, EVALUE, DOMAINS
        FROM INTERPRO.CLAN_MATCH
        """
    )

    for query, target, evalue, clob in cur:
        # DOMAINS is a LOB object: need to call read()

        domains = []
        for start, end in json.loads(clob.read()):
            domains.append({
                "start": start,
                "end": end
            })

        yield query, target, evalue, domains


def export_clans(ipr_url: str, pfam_url: str, p_clans: str, p_alignments: str,
                 **kwargs):
    buffer_size = kwargs.get("buffer_size", 1000000)
    threshold = kwargs.get("threshold", 1e-2)

    logger.info("loading clans")
    con = cx_Oracle.connect(ipr_url)
    cur = con.cursor()
    clans = get_clans(cur)

    clan_links = {}
    entry2clan = {}
    for accession, clan in clans.items():
        clan_links[accession] = {}
        for member_acc, score, seq_length in clan["members"]:
            entry2clan[member_acc] = (accession, seq_length)

    logger.info("exporting alignments")
    with DumpFile(p_alignments, compress=True) as df:
        i = 0
        alignments = []
        for query_acc, target_acc, evalue, domains in iter_alignments(cur):
            i += 1
            if not i % 10000000:
                logger.info(f"{i:>12,}")

            try:
                query_clan_acc, seq_length = entry2clan[query_acc]
            except KeyError:
                continue

            if evalue > threshold:
                continue

            try:
                target_clan_acc, _ = entry2clan[target_acc]
            except KeyError:
                target_clan_acc = None

            alignments.append((
                query_clan_acc,
                query_acc,
                target_acc,
                target_clan_acc,
                evalue,
                seq_length,
                json.dumps(domains)
            ))

            if len(alignments) == buffer_size:
                df.dump(alignments)
                alignments = []

            if query_clan_acc == target_clan_acc:
                # Query and target from the same clan: update the clan's links
                links = clan_links[query_clan_acc]

                if query_acc > target_acc:
                    query_acc, target_acc = target_acc, query_acc

                try:
                    targets = links[query_acc]
                except KeyError:
                    links[query_acc] = {target_acc: evalue}
                else:
                    if target_acc not in targets or evalue < targets[target_acc]:
                        targets[target_acc] = evalue

        df.dump(alignments)
        alignments = []
        logger.info(f"{i:>12,}")

    cur.close()
    con.close()

    logger.info("loading additional details for Pfam clans")
    pfam_clans = pfam.get_clans(pfam_url)

    logger.info("finalizing")
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

        if clan_acc in pfam_clans:
            # Replace `description`, add `authors` and `literature`
            clan.update(pfam_clans[clan_acc])

    dumpobj(p_clans, clans)
    logger.info("complete")
