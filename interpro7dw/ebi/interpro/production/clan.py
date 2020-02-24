# -*- coding: utf-8 -*-

import json
from typing import Dict

import cx_Oracle


class Clan(object):
    def __init__(self, accession: str, name: str, desc: str, database: str):
        self.accession = accession
        self.name = name
        self.description = desc
        self.database = database
        self.members = []
        self.links = {}

    def add_link(self, query_acc: str, target_acc: str, score: float):
        if query_acc > target_acc:
            query_acc, target_acc = target_acc, query_acc

        try:
            links = self.links[query_acc]
        except KeyError:
            self.links[query_acc] = {target_acc: score}
        else:
            if target_acc not in links or score < links[target_acc]:
                links[target_acc] = score

    def astuple(self) -> tuple:
        nodes = []
        for accession, score, seq_length in self.members:
            nodes.append({
                "accession": accession,
                "type": "entry",
                "score": score
            })

        links = []
        for query_acc, targets in self.links.items():
            for target_acc, score in targets.items():
                links.append({
                    "source": query_acc,
                    "target": target_acc,
                    "score": score
                })

        return (
            self.accession,
            self.name,
            self.description,
            self.database,
            1,
            json.dumps({
                "nodes": nodes,
                "links": links
            })
        )


def get_clans(url: str) -> Dict[str, Clan]:
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
            c = clans[accession] = Clan(accession, name, descr, database)

        c.members.append((member_acc, score, seq_length))

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
