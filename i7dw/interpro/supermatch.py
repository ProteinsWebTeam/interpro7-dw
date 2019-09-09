import json
from typing import Optional

from . import repr_frag, mysql
from .. import dbms, io, logger


class Supermatch(object):
    def __init__(self, entry_ac: str, entry_root: Optional[str],
                 fragments: list):
        self.members = {(entry_ac, entry_root)}
        self.fragments = fragments
        self.start = fragments[0]["start"]
        self.end = fragments[-1]["end"]

    @property
    def entries(self):
        for entry_ac, entry_root in self.members:
            yield entry_ac

    @property
    def fragments_str(self):
        return ','.join(["{start}-{end}".format(**frag)
                         for frag in self.fragments])

    @staticmethod
    def overlaps(start1, stop1, start2, stop2, min_overlap):
        overlap = min(stop1, stop2) - max(start1, start2) + 1
        shortest = min(stop1 - start1, stop2 - start2) + 1
        return overlap >= shortest * min_overlap / 100

    def merge_if_overlap(self, other, min_overlap):
        for acc1, root1 in self.members:
            for acc2, root2 in other.members:
                if root1 != root2:
                    return False

        if self.overlaps(self.start, self.end, other.start, other.end,
                         min_overlap):
            self.members |= other.members
            self.start = min(self.start, other.start)
            self.end = max(self.end, other.end)

            fragments = []
            for f1 in sorted(self.fragments + other.fragments,
                             key=repr_frag):
                for f2 in fragments:
                    if self.overlaps(f1["start"], f1["end"], f2["start"],
                                     f2["end"], min_overlap):
                        f2["start"] = min(f1["start"], f2["start"])
                        f2["end"] = max(f1["end"], f2["end"])
                        break
                else:
                    fragments.append(f1)
            self.fragments = fragments

            return True
        return False

    def __eq__(self, other):
        return self.start == other.start and self.end == other.end

    def __ne__(self, other):
        return not self == other

    def __lt__(self, other):
        return self.start < other.start or self.end < other.end

    def __le__(self, other):
        return self == other or self < other

    def __gt__(self, other):
        return self.start > other.start or self.end > other.end

    def __ge__(self, other):
        return self == other or self > other


def merge_supermatches(supermatches, min_overlap=20):
    merged_supermatches = []

    for sm in sorted(supermatches):
        for other in merged_supermatches:
            if other.merge_if_overlap(sm, min_overlap):
                break
        else:
            merged_supermatches.append(sm)

    return merged_supermatches


def intersect(entries: dict, counts: dict, intersections: dict):
    for acc1 in entries:
        try:
            counts[acc1] += 1
        except KeyError:
            counts[acc1] = 1

        for acc2 in entries:
            if acc1 >= acc2:
                continue
            elif acc1 in intersections:
                try:
                    overlaps = intersections[acc1][acc2]
                except KeyError:
                    overlaps = intersections[acc1][acc2] = [0, 0]
            else:
                intersections[acc1] = {acc2: [0, 0]}
                overlaps = intersections[acc1][acc2]

            flag = 0
            for f1 in entries[acc1]:
                len1 = f1["end"] - f1["start"] + 1

                for f2 in entries[acc2]:
                    len2 = f2["end"] - f2["start"] + 1
                    o = min(f1["end"], f2["end"]) - max(f1["start"], f2["start"]) + 1

                    if not flag & 1 and o >= len1 * 0.5:
                        flag |= 1
                        overlaps[0] += 1

                    if not flag & 2 and o >= len2 * 0.5:
                        flag |= 2
                        overlaps[1] += 1

                    if flag == 3:
                        break

                if flag == 3:
                    break


def calculate_relationships(my_uri: str, src_proteins: str, src_matches: str,
                            threshold: float, min_overlap: int=20,
                            ora_uri: str=None):
    logger.info("starting")
    entries = mysql.entry.get_entries(my_uri)
    proteins = io.Store(src_proteins)
    protein2matches = io.Store(src_matches)
    supfam = "homologous_superfamily"
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
                ENTRY_AC VARCHAR2(9) NOT NULL,
                FRAGMENTS VARCHAR(400) NOT NULL
            ) NOLOGGING
            """
        )

        cur.close()
        con.close()

        query = """
            INSERT /*+APPEND*/ INTO INTERPRO.SUPERMATCH2
            VALUES (:1, :2, :3)
        """
        con, cur = dbms.connect(ora_uri)
        cur.close()
        table = dbms.Populator(con, query, autocommit=True)
    else:
        table = None

    n_proteins = 0
    counts = {}
    interesctions = {}
    for acc, protein in proteins:
        matches = protein2matches.get(acc)
        if not matches:
            continue

        supermatches = []
        for m in matches:
            method_ac = m["method_ac"]
            entry_ac = entries[method_ac]["integrated"]

            if entry_ac:
                supermatches.append(
                    Supermatch(
                        entry_ac,
                        entries[entry_ac]["root"],
                        m["fragments"]
                    )
                )

        # Merge overlapping supermatches
        entry_matches = {}
        for sm in merge_supermatches(supermatches, min_overlap):
            for entry_ac in sm.entries:
                if table:
                    table.insert((acc, entry_ac, sm.fragments_str))

                entry_matches[entry_ac] = sm.fragments

        intersect(entry_matches, counts, interesctions)

        n_proteins += 1
        if not n_proteins % 10000000:
            logger.info("{:>12,}".format(n_proteins))

    logger.info("{:>12,}".format(n_proteins))

    proteins.close()
    protein2matches.close()

    if table:
        table.close()
        logger.info("{} supermatches inserted".format(table.count))

        logger.info("indexing SUPERMATCH2")
        cur = con.cursor()
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
    for acc1 in interesctions:
        cnt1 = counts[acc1]

        for acc2, (o1, o2) in interesctions[acc1].items():
            cnt2 = counts[acc2]

            # Independent coefficients
            coef1 = o1 / (cnt1 + cnt2 - o1)
            coef2 = o2 / (cnt1 + cnt2 - o2)

            # Final coefficient: average of independent coefficients
            coef = (coef1 + coef2) * 0.5

            # Containment indices
            c1 = o1 / cnt1
            c2 = o2 / cnt2

            if any([item >= threshold for item in (coef, c1, c2)]):
                e1 = entries[acc1]
                e2 = entries[acc2]
                t1 = e1["type"]
                t2 = e2["type"]

                if (t1 == supfam and t2 in types) or (t2 == supfam and t1 in types):
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

                    for k, v in [(acc1, e2), (acc2, e1)]:
                        try:
                            overlapping[k].append(v)
                        except KeyError:
                            overlapping[k] = [v]

    logger.info("updating table")
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

    logger.info("complete")
