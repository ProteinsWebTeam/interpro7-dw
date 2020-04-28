# -*- coding: utf-8 -*-

import os
from typing import Optional, Sequence

import cx_Oracle
import MySQLdb

from interpro7dw import logger
from interpro7dw.ebi import pfam
from interpro7dw.ebi.interpro.utils import Table
from interpro7dw.ebi.interpro.utils import overlaps_pdb_chain, repr_fragment
from interpro7dw.utils import DumpFile, DirectoryTree, Store, merge_dumps
from interpro7dw.utils import loadobj, url2dict
from .utils import jsonify, reduce


def init_xrefs() -> dict:
    return {
        "domain_architectures": set(),
        "matches": 0,
        "proteins": set(),
        "proteomes": set(),
        "structures": set(),
        "taxa": set()
    }


def dump_xrefs(xrefs: dict, output: str):
    with DumpFile(output) as f:
        for entry_acc in sorted(xrefs):
            f.dump((entry_acc, xrefs[entry_acc]))


def insert_entries(stg_url: str, p_entries: str, p_proteins: str,
                   p_structures: str, p_uniprot2ida: str,
                   p_uniprot2matches: str, p_uniprot2proteome: str,
                   p_entry2xrefs: str, **kwargs):
    dir = kwargs.get("dir")
    min_overlap = kwargs.get("overlap", 0.2)
    min_similarity = kwargs.get("similarity", 0.75)
    pro_url = kwargs.get("pro_url")

    logger.info("preparing data")
    dt = DirectoryTree(dir)
    uniprot2pdbe = {}
    for pdb_id, entry in loadobj(p_structures).items():
        for uniprot_acc, chains in entry["proteins"].items():
            try:
                uniprot2pdbe[uniprot_acc][pdb_id] = chains
            except KeyError:
                uniprot2pdbe[uniprot_acc] = {pdb_id: chains}

    entries = loadobj(p_entries)
    proteins = Store(p_proteins)
    u2ida = Store(p_uniprot2ida)
    u2matches = Store(p_uniprot2matches)
    u2proteome = Store(p_uniprot2proteome)

    logger.info("exporting cross-references")
    i = 0
    xrefs = {}
    files = []
    entry_counts = {}
    entry_intersections = {}
    for uniprot_acc, matches in u2matches.items():
        info = proteins[uniprot_acc]

        try:
            dom_members, dom_arch, dom_arch_id = u2ida[uniprot_acc]
        except KeyError:
            dom_members = []

        try:
            proteome_id = u2proteome[uniprot_acc]
        except KeyError:
            proteome_id = None

        pdb_entries = uniprot2pdbe.get(uniprot_acc, {})

        supermatches = []
        for entry_acc, locations in matches.items():
            try:
                entry_xref = xrefs[entry_acc]
            except KeyError:
                entry_xref = xrefs[entry_acc] = init_xrefs()

            if entry_acc in dom_members:
                entry_xref["domain_architectures"].add(dom_arch_id)

            entry_xref["matches"] += len(locations)
            entry_xref["proteins"].add((uniprot_acc, info["identifier"]))

            if proteome_id:
                entry_xref["proteomes"].add(proteome_id)

            for pdb_id, chains in pdb_entries.items():
                for chain_id, segments in chains.items():
                    if overlaps_pdb_chain(locations, segments):
                        entry_xref["structures"].add(pdb_id)
                        break  # Skip other chains

            entry_xref["taxa"].add(info["taxid"])

            # Create a Supermatch for each integrated signature match
            entry = entries[entry_acc]
            if entry.integrated_in:
                # Integrated member DB signature
                interpro_acc = entry.integrated_in
                root = entries[interpro_acc].hierarchy["accession"]
                for loc in locations:
                    sm = Supermatch(interpro_acc, loc["fragments"], root)
                    supermatches.append(sm)

        # Merge overlapping supermatches
        merged = []
        for sm_to_merge in sorted(supermatches):
            for sm_merged in merged:
                if sm_merged.overlaps(sm_to_merge, min_overlap):
                    """
                    Supermatches overlap
                        (sm_to_merge has been merged into sm_merged)
                    """
                    break
            else:
                # sm_to_merge does not overlap with any other supermatches
                merged.append(sm_to_merge)

        # Group by entry and populate table
        merged_grouped = {}
        for sm in merged:
            # fragments = sm.stringify_fragments()
            for interpro_acc in sm.entries:
                # table.insert((uniprot_acc, interpro_acc, fragments))

                try:
                    merged_grouped[interpro_acc] += sm.fragments
                except KeyError:
                    merged_grouped[interpro_acc] = list(sm.fragments)

        # Evaluate how entries overlap
        for interpro_acc, fragments1 in merged_grouped.items():
            try:
                entry_counts[interpro_acc] += 1
            except KeyError:
                entry_counts[interpro_acc] = 1

            for other_acc, fragments2 in merged_grouped.items():
                if other_acc >= interpro_acc:
                    continue

                try:
                    obj = entry_intersections[interpro_acc]
                except KeyError:
                    obj = entry_intersections[interpro_acc] = {
                        other_acc: [0, 0]
                    }

                try:
                    overlaps = obj[other_acc]
                except KeyError:
                    overlaps = obj[other_acc] = [0, 0]

                flag = 0
                for f1 in fragments1:
                    start1 = f1["start"]
                    end1 = f1["end"]
                    length1 = end1 - start1 + 1

                    for f2 in fragments2:
                        start2 = f2["start"]
                        end2 = f2["end"]
                        length2 = end2 - start2 + 1
                        overlap = min(end1, end2) - max(start1,
                                                        start2) + 1

                        if not flag & 1 and overlap >= length1 * 0.5:
                            # 1st time fragments overlap >= 50% of f1
                            flag |= 1
                            overlaps[0] += 1

                        if not flag & 2 and overlap >= length2 * 0.5:
                            # 1st time fragments overlap >= 50% of f2
                            flag |= 2
                            overlaps[1] += 1

                    if flag == 3:
                        # Both cases already happened, no need to continue iterating
                        break

        i += 1
        if not i % 1000000:
            output = dt.mktemp()
            dump_xrefs(xrefs, output)
            files.append(output)
            xrefs = {}

            if not i % 10000000:
                logger.info(f"{i:>12,}")

    if xrefs:
        output = dt.mktemp()
        dump_xrefs(xrefs, output)
        files.append(output)
        xrefs = {}

    logger.info(f"{i:>12,}")
    logger.info(f"temporary files: "
                f"{sum(map(os.path.getsize, files))/1024/1024:.0f} MB")

    proteins.close()
    u2ida.close()
    u2matches.close()
    u2proteome.close()

    logger.info("calculating Jaccard coefficients")
    supfam = "homologous_superfamily"
    types = (supfam, "domain", "family", "repeat")
    overlapping = {}
    for entry_acc, overlaps in entry_intersections.items():
        entry1 = entries[entry_acc]
        entry_cnt = entry_counts[entry_acc]

        for other_acc, (o1, o2) in overlaps.items():
            other_cnt = entry_counts[other_acc]

            # Independent coefficients
            coef1 = o1 / (entry_cnt + other_cnt - o1)
            coef2 = o2 / (entry_cnt + other_cnt - o2)

            # Final coefficient: average of independent coefficients
            coef = (coef1 + coef2) * 0.5

            # Containment indices
            c1 = o1 / entry_cnt
            c2 = o2 / other_cnt

            if all([item < min_similarity for item in (coef, c1, c2)]):
                continue

            # Entries are similar enough

            entry2 = entries[other_acc]
            if ((entry1.type == supfam and entry2.type in types)
                    or (entry2.type == supfam and entry1.type in types)):
                # e1 -> e2 relationship
                try:
                    obj = overlapping[entry_acc]
                except KeyError:
                    obj = overlapping[entry_acc] = []
                finally:
                    obj.append({
                        "accession": other_acc,
                        "name": entry2.name,
                        "type": entry2.type
                    })

                # e2 -> e1 relationship
                try:
                    obj = overlapping[other_acc]
                except KeyError:
                    obj = overlapping[other_acc] = []
                finally:
                    obj.append({
                        "accession": entry_acc,
                        "name": entry1.name,
                        "type": entry1.type
                    })

    logger.info("populating webfront_entry")
    con = MySQLdb.connect(**url2dict(stg_url))
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
            literature LONGTEXT,
            hierarchy LONGTEXT,
            cross_references LONGTEXT,
            interactions LONGTEXT,
            pathways LONGTEXT,
            overlaps_with LONGTEXT,
            is_featured TINYINT NOT NULL,
            is_alive TINYINT NOT NULL,
            history LONGTEXT,
            entry_date DATETIME NOT NULL,
            deletion_date DATETIME,
            counts LONGTEXT
        ) CHARSET=utf8 DEFAULT COLLATE=utf8_unicode_ci
        """
    )
    cur.close()

    sql = """
        INSERT INTO webfront_entry
        VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s,
          %s, %s, %s, %s, %s, %s, %s, %s)
    """

    with Table(con, sql) as table:
        with DumpFile(p_entry2xrefs, compress=True) as f:
            for accession, xrefs in merge_dumps(files):
                entry = entries.pop(accession)
                counts = reduce(xrefs)
                counts.update({
                    "interactions": len(entry.ppi),
                    "pathways": sum([len(v) for v in entry.pathways.values()]),
                    "sets": 1 if entry.clan else 0
                })

                table.insert((
                    None,
                    accession,
                    entry.type,
                    entry.name,
                    entry.short_name,
                    entry.database,
                    jsonify(entry.integrates),
                    entry.integrated_in,
                    jsonify(entry.go_terms),
                    jsonify(entry.description),
                    jsonify(entry.wikipedia),
                    jsonify(entry.literature),
                    jsonify(entry.hierarchy),
                    jsonify(entry.cross_references),
                    jsonify(entry.ppi),
                    jsonify(entry.pathways),
                    jsonify(overlapping.get(accession)),
                    0,
                    0 if entry.is_deleted else 1,
                    jsonify(entry.history),
                    entry.creation_date,
                    entry.deletion_date,
                    jsonify(counts)
                ))

                # We don't need matches in cross-references
                del xrefs["matches"]
                f.dump((accession, xrefs))

        # Entries not matching any proteins
        for entry in entries.values():
            counts = reduce(init_xrefs())
            counts["interactions"] = 0
            counts["pathways"] = 0
            counts["sets"] = 1 if entry.clan else 0
            table.insert((
                None,
                entry.accession,
                entry.type,
                entry.name,
                entry.short_name,
                entry.database,
                jsonify(entry.integrates),
                entry.integrated_in,
                jsonify(entry.go_terms),
                jsonify(entry.description),
                jsonify(entry.wikipedia),
                jsonify(entry.literature),
                jsonify(entry.hierarchy),
                jsonify(entry.cross_references),
                jsonify(entry.ppi),
                None,  # No pathways (because no proteins/enzymes)
                jsonify(overlapping.get(entry.accession)),
                0,
                0 if entry.is_deleted else 1,
                jsonify(entry.history),
                entry.creation_date,
                entry.deletion_date,
                jsonify(counts)
            ))

    con.commit()
    con.close()
    dt.remove()
    logger.info("complete")


def insert_annotations(pfam_url: str, stg_url: str):
    con = MySQLdb.connect(**url2dict(stg_url))
    cur = con.cursor()
    cur.execute("DROP TABLE IF EXISTS webfront_entryannotation")
    cur.execute(
        """
        CREATE TABLE webfront_entryannotation
        (
            annotation_id VARCHAR(255) PRIMARY KEY NOT NULL,
            accession_id VARCHAR(25) NOT NULL,
            type VARCHAR(32) NOT NULL,
            value LONGBLOB NOT NULL,
            mime_type VARCHAR(32) NOT NULL
        ) CHARSET=utf8 DEFAULT COLLATE=utf8_unicode_ci
        """
    )
    cur.close()

    sql = """
        INSERT INTO webfront_entryannotation
        VALUES (%s, %s, %s, %s, %s)
    """

    with Table(con, sql) as table:
        for acc, anno_type, value, mime in pfam.get_annotations(pfam_url):
            anno_id = f"{acc}--{anno_type}"
            table.insert((anno_id, acc, anno_type, value, mime))

    con.commit()
    con.close()


class Supermatch(object):
    def __init__(self, acc: str, frags: Sequence[dict], root: Optional[str]):
        self.members = {(acc, root)}
        self.fragments = frags
        """
        frags is sorted by (start, end):
          - start of the first frag is guaranteed to be the leftmost one
          - end of the last frag is NOT guaranteed to be the rightmost one
             (e.g. [(5, 100), (6, 80)])
        """
        self.start = frags[0]["start"]
        self.end = max(f["end"] for f in frags)

    def __eq__(self, other) -> bool:
        return self.start == other.start and self.end == other.end

    def __ne__(self, other) -> bool:
        return not self == other

    def __lt__(self, other) -> bool:
        return self.start < other.start or self.end < other.end

    def __le__(self, other) -> bool:
        return self == other or self < other

    def __gt__(self, other) -> bool:
        return self.start > other.start or self.end > other.end

    def __ge__(self, other) -> bool:
        return self == other or self > other

    @property
    def entries(self):
        return [acc for acc, root in self.members]

    def stringify_fragments(self):
        return ','.join(["{start}-{end}".format(**f) for f in self.fragments])

    def overlaps(self, other, min_overlap: float) -> bool:
        for acc1, root1 in self.members:
            for acc2, root2 in other.members:
                if root1 != root2:
                    return False

        # All members are in the same hierarchy
        overlap = min(self.end, other.end) - max(self.start, other.start) + 1
        shortest = min(self.end - self.start, other.end - other.start) + 1

        if overlap < shortest * min_overlap:
            # Supermatches do not significantly overlap
            return False

        # Merge supermatches
        self.members |= other.members
        self.start = min(self.start, other.start)
        self.end = max(self.end, other.end)

        # Merge fragments
        fragments = []
        for f1 in sorted(self.fragments+other.fragments, key=repr_fragment):
            start1 = f1["start"]
            end1 = f1["end"]

            for f2 in fragments:
                start2 = f2["start"]
                end2 = f2["end"]
                overlap = min(end1, end2) - max(start1, start2) + 1
                shortest = min(end1 - start1, end2 - start2) + 1

                if overlap >= shortest * min_overlap:
                    # Merge f1 into f2
                    f2["start"] = min(start1, start2)
                    f2["end"] = max(end1, end2)
                    break
            else:
                # f1 does not overlap with any other fragments
                fragments.append(f1)

        self.fragments = fragments
        return True


class SupermatchTable(object):
    def __init__(self, url: Optional[str]):
        if url:
            con = cx_Oracle.connect(url)
            self.create(con)

            sql = """
                INSERT /*+APPEND*/ INTO INTERPRO.SUPERMATCH2
                VALUES (:1, :2, :3)
            """
            self.table = Table(con, sql, autocommit=True)
            self.insert = self._insert
        else:
            self.table = None
            self.insert = self._do_nothing

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __del__(self):
        self.close()

    def _insert(self, record):
        self.table.insert(record)

    def _do_nothing(self, record):
        return

    @staticmethod
    def create(con: cx_Oracle.Connection):
        cur = con.cursor()

        try:
            cur.execute("DROP TABLE INTERPRO.SUPERMATCH2 PURGE")
        except cx_Oracle.DatabaseError:
            pass

        cur.execute(
            """
            CREATE TABLE INTERPRO.SUPERMATCH2 (
                PROTEIN_AC VARCHAR2(15) NOT NULL,
                ENTRY_AC VARCHAR2(9) NOT NULL,
                FRAGMENTS VARCHAR(400) NOT NULL
            ) NOLOGGING
            """
        )
        cur.close()

    def close(self):
        if self.table is None:
            return

        self.table.close()
        logger.info(f"{self.table.count:,} supermatches inserted")

        logger.info("indexing table")
        cur = self.table.con.cursor()
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
        cur.execute("GRANT SELECT ON INTERPRO.SUPERMATCH2 TO INTERPRO_SELECT")
        cur.close()
        self.table.con.close()
        self.table = None
