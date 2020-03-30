# -*- coding: utf-8 -*-

import json
from typing import Dict

import cx_Oracle


def get_clans(url: str) -> Dict[str, dict]:
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    cur.execute(
        """
        SELECT
          C.CLAN_AC, C.NAME, C.DESCRIPTION, LOWER(D.DBSHORT), M.METHOD_AC, 
          LENGTH(M.SEQ), M.SCORE
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

    cur.close()
    con.close()

    return clans


def iter_clan_alignments(url: str):
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    cur.execute(
        """
        SELECT QUERY_AC, TARGET_AC, EVALUE, DOMAINS
        FROM INTERPRO.CLAN_MEMBER_ALN
        """
    )

    for query, target, evalue, clob in cur:
        # DOMAINS is a LOB object: need to call read()
        obj = json.loads(clob.read())
        domains = []

        for dom in obj:
            # Do not use query/target sequences and iEvalue
            domains.append({
                "start": dom["start"],
                "end": dom["end"]
            })

        yield query, target, evalue, domains

    cur.close()
    con.close()
