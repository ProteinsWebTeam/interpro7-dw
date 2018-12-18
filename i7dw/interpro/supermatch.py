#!/usr/bin/env python
# -*- coding: utf-8 -*-

import json
import logging
import re
import time

from . import mysql
from .. import dbms, io


class Supermatch(object):
    def __init__(self, entry_ac, entry_root, start, end):
        self.entries = {(entry_ac, entry_root)}
        self.start = start
        self.end = end

    def in_same_hierarchy(self, other):
        for acc_1, root_1 in self.entries:
            for acc_2, root_2 in other.entries:
                if root_1 != root_2:
                    return False

        return True

    def merge_if_overlap(self, other, min_overlap):
        overlap = min(self.end, other.end) - max(self.start, other.start) + 1
        shortest = min(self.end - self.start, other.end - other.start) + 1

        if (overlap / shortest * 100) >= min_overlap:
            self.entries = self.entries | other.entries
            self.start = min(self.start, other.start)
            self.end = max(self.end, other.end)
            return self
        else:
            return None

    def format(self):
        return '{}:{:.0f}-{:.0f}'.format(
            self.format_entries(),
            self.start,
            self.end
        )

    def format_entries(self):
        return '&'.join(sorted([re.sub(r'IPR0*', '', accession) for accession in self.get_entries()]))

    def get_entries(self):
        return [entry_ac for entry_ac, entry_root in self.entries]

    def __eq__(self, other):
        return (
                isinstance(other, Supermatch) and
                self.start == other.start and
                self.end == other.end and
                self.entries == other.entries
        )


class SupermatchSet(object):
    def __init__(self, supermatch):
        self.supermatches = [supermatch]

    def add(self, candidate, min_overlap):
        if not self.supermatches[0].in_same_hierarchy(candidate):
            return False

        merged = None
        for sm in self.supermatches:
            merged = sm.merge_if_overlap(candidate, min_overlap)

            if merged is not None:
                break

        if merged is None:
            self.supermatches.append(candidate)
        else:
            # Merged supermatch: we now need to remove overlaps between the newly merged supermatch and others
            indexes_ok = set()

            while True:
                index = None

                for i, sm in enumerate(self.supermatches):
                    if sm == merged or i in indexes_ok:
                        continue

                    if merged.merge_if_overlap(sm, min_overlap):
                        # Overlap so merged, we now have to remove the merged supermatch (sm)
                        index = i
                        break
                    else:
                        # No overlap, might be skipped during next iteration
                        indexes_ok.add(i)

                if index is None:
                    # No move overlaps
                    break
                else:
                    self.supermatches.pop(index)

        return True


def merge_supermatches(supermatches, min_overlap=20):
    sets = []

    for sm in supermatches:
        in_set = False

        for s in sets:
            in_set = s.add(sm, min_overlap)

            if in_set:
                break

        if not in_set:
            sets.append(SupermatchSet(sm))

    return sets


def calculate_relationships(my_uri: str, src_proteins: str, src_matches: str,
                            threshold: float, min_overlap: int=20,
                            ora_uri: str=None):
    logging.info("starting")
    entries = mysql.get_entries(my_uri)
    proteins = io.Store(src_proteins)
    protein2matches = io.Store(src_matches)
    types = ("homologous_superfamily", "domain", "family", "repeat")

    if ora_uri:
        con, cur = dbms.connect(ora_uri)

        try:
            cur.execute(
                """
                DROP TABLE INTERPRO.SUPERMATCH2
                CASCADE CONSTRAINTS
                """
            )
        except:
            pass

        cur.execute(
            """
            CREATE TABLE INTERPRO.SUPERMATCH2
            (
                PROTEIN_AC VARCHAR2(15) NOT NULL,
                DBCODE CHAR(1) NOT NULL,
                ENTRY_AC VARCHAR2(9) NOT NULL,
                POS_FROM NUMBER(5) NOT NULL,
                POS_TO NUMBER(5) NOT NULL
            ) NOLOGGING
            """
        )

        cur.close()
        con.close()

        query = """
            INSERT /*+APPEND*/ INTO INTERPRO.SUPERMATCH2 (
              PROTEIN_AC, DBCODE, ENTRY_AC, POS_FROM, POS_TO
            )
            VALUES (:1, :2, :3, :4, :5)
        """
        table = dbms.Populator(ora_uri, query, autocommit=True)
    else:
        table = None

    ts = time.time()
    n_proteins = 0
    sets = {}
    overlaps = {}
    for acc, protein in proteins:
        matches = protein2matches.get(acc)
        if not matches:
            continue

        dbcode = 'S' if protein["isReviewed"] else 'T'
        supermatches = []
        for m in matches:
            method_ac = m["method_ac"]
            entry_ac = entries[method_ac]["integrated"]

            if entry_ac:
                pos_start = None
                pos_end = None
                for f in m["fragments"]:
                    if pos_start is None or f["start"] < pos_start:
                        pos_start = f["start"]
                    if pos_end is None or f["end"] > pos_end:
                        pos_end = f["end"]

                supermatches.append(
                    Supermatch(
                        entry_ac,
                        entries[entry_ac]["root"],
                        pos_start,
                        pos_end
                    )
                )

        # Merge overlapping supermatches
        sm_sets = merge_supermatches(supermatches, min_overlap)
        supermatches = {}
        for s in sm_sets:
            for sm in s.supermatches:
                for entry_ac in sm.get_entries():
                    if table:
                        table.insert((acc, dbcode, entry_ac.upper(),
                                      sm.start, sm.end))

                    # Current implementation: leftmost match only
                    if entry_ac not in supermatches:
                        supermatches[entry_ac] = [(sm.start, sm.end)]

        intersect(supermatches, sets, overlaps)

        n_proteins += 1
        if not n_proteins % 1000000:
            logging.info("{:>12,} ({:.0f} proteins/sec)".format(
                n_proteins, n_proteins / (time.time() - ts)
            ))

    logging.info("{:>12,} ({:.0f} proteins/sec)".format(
        n_proteins, n_proteins / (time.time() - ts)
    ))

    proteins.close()
    protein2matches.close()

    if table:
        table.close()
        logging.info("{} supermatches inserted".format(table.count))

        logging.info("indexing SUPERMATCH2")
        con, cur = dbms.connect(ora_uri)

        cur.execute(
            """
            ALTER TABLE INTERPRO.SUPERMATCH2
            ADD CONSTRAINT PK_SUPERMATCH2
            PRIMARY KEY (PROTEIN_AC, ENTRY_AC, POS_FROM, POS_TO, DBCODE)
            """
        )

        try:
            cur.execute(
                """
                ALTER TABLE INTERPRO.SUPERMATCH2
                ADD CONSTRAINT FK_SUPERMATCH2$PROTEIN_AC
                FOREIGN KEY (PROTEIN_AC)
                REFERENCES INTERPRO.PROTEIN (PROTEIN_AC)
                ON DELETE CASCADE
                """
            )
        except:
            pass

        try:
            cur.execute(
                """
                ALTER TABLE INTERPRO.SUPERMATCH2
                ADD CONSTRAINT FK_SUPERMATCH2$ENTRY_AC
                FOREIGN KEY (ENTRY_AC)
                REFERENCES INTERPRO.ENTRY (ENTRY_AC)
                ON DELETE CASCADE
                """
            )
        except:
            pass

        # Indexes
        cur.execute(
            """
            CREATE INDEX I_SUPERMATCH2$PROTEIN
            ON INTERPRO.SUPERMATCH2 (PROTEIN_AC)
            NOLOGGING
            """
        )
        cur.execute(
            """
            CREATE INDEX I_SUPERMATCH2$ENTRY
            ON INTERPRO.SUPERMATCH2 (ENTRY_AC)
            NOLOGGING
            """
        )
        cur.execute(
            """
            CREATE INDEX I_SUPERMATCH2$DBCODE$ENTRY
            ON INTERPRO.SUPERMATCH2 (DBCODE, ENTRY_AC)
            NOLOGGING
            """
        )

        # Statistics
        cur.execute(
            """
                BEGIN
                    DBMS_STATS.GATHER_TABLE_STATS(:1, :2, cascade => TRUE);
                END;
            """,
            ("INTERPRO", "SUPERMATCH2")
        )

        # Privileges
        cur.execute(
            """
            GRANT SELECT
            ON INTERPRO.SUPERMATCH2
            TO INTERPRO_SELECT
            """
        )

        cur.close()
        con.close()

    # Compute Jaccard coefficients
    overlapping = {}
    for acc1 in overlaps:
        s1 = sets[acc1]

        for acc2 in overlaps[acc1]:
            s2 = sets[acc2]
            o1, o2 = overlaps[acc1][acc2]

            # Independent coefficients
            coef1 = o1 / (s1 + s2 - o1)
            coef2 = o2 / (s1 + s2 - o2)

            # Final coefficient: average of independent coefficients
            coef = (coef1 + coef2) * 0.5

            # Containment indices
            c1 = o1 / s1
            c2 = o2 / s2

            if any([item >= threshold for item in (coef, c1, c2)]):
                e1 = entries[acc1]
                e2 = entries[acc2]
                t1 = e1["type"]
                t2 = e2["type"]

                if t1 == "homologous_superfamily":
                    if t2 not in types:
                        continue
                elif t2 == "homologous_superfamily":
                    if t1 not in types:
                        continue
                else:
                    continue

                e1 = {
                    "accession": e1["accession"],
                    "name": e1["name"],
                    "type": e1["type"]
                }

                e2 = {
                    "accession": e2["accession"],
                    "name": e2["name"],
                    "type": e2["type"]
                }

                if acc1 in overlapping:
                    overlapping[acc1].append(e2)
                else:
                    overlapping[acc1] = [e2]

                if acc2 in overlapping:
                    overlapping[acc2].append(e1)
                else:
                    overlapping[acc2] = [e1]

    logging.info("updating table")
    con, cur = dbms.connect(my_uri)
    for acc in overlapping:
        cur.execute(
            """
            UPDATE webfront_entry
            SET overlaps_with = %s
            WHERE accession = %s
            """,
            (json.dumps(overlapping[acc]), acc)
        )

    con.commit()
    cur.close()
    con.close()

    logging.info("complete")


def intersect(matches: dict, sets: dict, intersections: dict):
    for acc1 in matches:
        if acc1 in sets:
            sets[acc1] += 1
        else:
            sets[acc1] = 1

        for acc2 in matches:
            if acc1 >= acc2:
                continue
            elif acc1 not in intersections:
                intersections[acc1] = {acc2: [0, 0]}
            elif acc2 not in intersections[acc1]:
                intersections[acc1][acc2] = [0, 0]

            m1 = matches[acc1][0]
            m2 = matches[acc2][0]
            o = min(m1[1], m2[1]) - max(m1[0], m2[0]) + 1

            l1 = m1[1] - m1[0] + 1
            l2 = m2[1] - m2[0] + 1

            if o > l1 * 0.5:
                # acc1 is in acc2 (because it overlaps acc2 at least 50%)
                intersections[acc1][acc2][0] += 1

            if o > l2 * 0.5:
                # acc2 is in acc1
                intersections[acc1][acc2][1] += 1
