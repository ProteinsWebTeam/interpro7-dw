# -*- coding: utf-8 -*-

import json
import os
from typing import Any, Dict, Generator, List, Optional

import cx_Oracle
import MySQLdb
import MySQLdb.cursors

from i7dw import cdd, logger, io, pfam
from i7dw.interpro import Table, extract_frag, oracle
from .utils import parse_url


def to_json(obj: Any) -> Optional[str]:
    return json.dumps(obj) if obj else None


def from_json(string: str, default: Optional[Any]=None):
    return json.loads(string) if isinstance(string, str) else default


def insert_entries(my_url: str, ora_url: str, pfam_url: str):
    wiki = pfam.get_wiki(pfam_url)

    query = """
        INSERT INTO webfront_entry (accession, type, name, short_name, 
                                    source_database, member_databases, 
                                    integrated_id, go_terms, description, 
                                    wikipedia, literature, hierarchy, 
                                    cross_references, is_alive, entry_date, 
                                    deletion_date)
        VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)    
    """

    con = MySQLdb.connect(**parse_url(my_url), charset="utf8")
    with Table(con, query) as table:
        for e in oracle.get_entries(ora_url):
            table.insert((
                e["accession"],
                e["type"],
                e["name"],
                e["short_name"],
                e["database"],
                to_json(e["member_databases"]),
                e["integrated"],
                to_json(e["go_terms"]),
                to_json(e["descriptions"]),
                to_json(wiki.get(e["accession"])),
                to_json(e["citations"]),
                to_json(e["hierarchy"]),
                to_json(e["cross_references"]),
                1,  # is alive
                e["date"],
                None  # deletion date
            ))

        for e in oracle.get_deleted_entries(ora_url):
            table.insert((
                e["accession"],
                e["type"],
                e["name"],
                e["short_name"],
                e["database"],
                None, None, None, None, None, None, None, None,
                0,
                e["creation_date"],
                e["deletion_date"]
            ))

    con.commit()
    con.close()


def insert_annotations(my_url: str, pfam_url: str):
    query = """
        INSERT INTO webfront_entryannotation (annotation_id, accession_id, 
                                              type, value, mime_type ) 
        VALUES (%s, %s, %s, %s, %s)    
    """

    con = MySQLdb.connect(**parse_url(my_url), charset="utf8")
    with Table(con, query) as table:
        for entry_acc, annotations in pfam.get_annotations(pfam_url).items():
            for annotation in annotations:
                table.insert((
                    f"{entry_acc}--{annotation['type']}",
                    entry_acc,
                    annotation["type"],
                    annotation["value"],
                    annotation["mime_type"]
                ))

    con.commit()
    con.close()


def insert_sets(my_url: str, ora_url: str, pfam_url: str):
    query1 = """
        INSERT INTO webfront_set (accession, name, description, 
                                  source_database, is_set, relationships
        ) 
        VALUES (%s, %s, %s, %s, %s, %s)
    """

    query2 = """
        INSERT INTO webfront_alignment (set_acc, entry_acc, target_acc, 
                                        target_set_acc, score, seq_length, 
                                        domains)
        VALUES (%s, %s, %s, %s, %s, %s, %s)
    """

    sets = pfam.get_clans(pfam_url)
    sets.update(cdd.get_superfamilies())

    con = MySQLdb.connect(**parse_url(my_url), charset="utf8")

    with Table(con, query1, buffer_size=1) as t1, Table(con, query2) as t2:
        for obj in oracle.get_profile_alignments(ora_url):
            set_acc = obj[0]
            source_database = obj[1]
            relationships = obj[2]
            alignments = obj[3]

            if source_database in ("cdd", "pfam"):
                try:
                    s = sets[set_acc]
                except KeyError:
                    logger.warning(f"unknown CDD/Pfam set: {set_acc}")
                    continue

                if source_database == "pfam":
                    # Use nodes from Pfam DB (they have a 'real' score)
                    relationships["nodes"] = s["relationships"]["nodes"]

                name = s["name"] or set_acc
                desc = s["description"]
            else:
                name = set_acc
                desc = None

            t1.insert((set_acc, name, desc, source_database, 1,
                       json.dumps(relationships)))

            for entry_acc, targets in alignments.items():
                for target in targets:
                    t2.insert((set_acc, entry_acc, *target))

    con.commit()
    con.close()


def find_node(node, accession, relations=list()):
    """
    Find a entry (node) in its hierarchy tree
    and store its relations (ancestors + direct children)
    """
    if node["accession"] == accession:
        relations += [child["accession"] for child in node["children"]]
        return node

    for child in node["children"]:
        child = find_node(child, accession, relations)

        if child:
            relations.append(node["accession"])
            return child

    return None


def get_entries(url: str) -> Dict[str, Dict]:
    con = MySQLdb.connect(**parse_url(url), charset="utf8")
    cur = con.cursor()
    cur.execute(
        """
            SELECT accession, source_database, entry_date, description,
                   integrated_id, name, type, short_name, member_databases,
                   go_terms, literature, cross_references, hierarchy
            FROM webfront_entry
            WHERE is_alive = 1
        """
    )

    entries = {}
    for row in cur:
        accession = row[0]
        relations = []
        hierarchy = from_json(row[12])
        if hierarchy:
            find_node(hierarchy, accession, relations)
            _hierarchy = hierarchy.get("accession")
        else:
            _hierarchy = None

        entries[accession] = {
            "accession": accession,
            "database": row[1],
            "date": row[2],
            "descriptions": from_json(row[3], default=[]),
            "integrated": row[4],
            "name": row[5],
            "type": row[6],
            "short_name": row[7],
            "member_databases": from_json(row[8], default={}),
            "go_terms": from_json(row[9], default=[]),
            "citations": from_json(row[10], default={}),
            "cross_references": from_json(row[11], default={}),
            "root": _hierarchy,
            "relations": relations
        }

    cur.close()
    con.close()

    return entries


def iter_sets(url: str) -> Generator[Dict, None, None]:
    con = MySQLdb.connect(**parse_url(url), charset="utf8")
    cur = MySQLdb.cursors.SSCursor(con)
    cur.execute(
        """
        SELECT accession, name, description, source_database, relationships
        FROM webfront_set
        """
    )

    for row in cur:
        yield {
            "accession": row[0],
            "name": row[1],
            "description": row[2],
            "database": row[3],
            "members": [n["accession"] for n in json.loads(row[4])["nodes"]]
        }

    cur.close()
    con.close()


class Supermatch(object):
    def __init__(self, entry_acc: str, entry_root: Optional[str],
                 fragments: List[Dict]):
        self.members = {(entry_acc, entry_root)}
        self.fragments = fragments
        """
        `fragments` are sorted by (start, end):
          - `start` of the first fragm is guaranteed to be the leftmost one
          - `end` of the last frag is NOT guaranteed to be the rightmost one
             (e.g. [(5, 100), (6, 80)])
        """
        self.start = fragments[0]["start"]
        self.end = max(f["end"] for f in fragments)

    @property
    def entries(self) -> Generator[str, None, None]:
        for entry_acc, entry_root in self.members:
            yield entry_acc

    @property
    def fragments_str(self) -> str:
        return ','.join(["{start}-{end}".format(**f) for f in self.fragments])

    def overlaps(self, other, min_overlap: float) -> bool:
        for acc1, root1 in self.members:
            for acc2, root2 in other.members:
                if root1 != root2:
                    return False

        overlap = min(self.end, other.end) - max(self.start, other.start) + 1
        shortest = min(self.end - self.start, other.end - other.start) + 1

        if overlap < shortest * min_overlap:
            # Supermatches do not significantly overlap
            return False

        # Supermatches significantly overlap
        self.members |= other.members
        self.start = min(self.start, other.start)
        self.end = max(self.end, other.end)

        # Merge fragments
        fragments = []
        for f1 in sorted(self.fragments+other.fragments, key=extract_frag):
            s1 = f1["start"]
            e1 = f1["end"]
            for f2 in fragments:
                s2 = f2["start"]
                e2 = f2["end"]
                overlap = min(e1, e2) - max(s1, s2) + 1
                shortest = min(e1 - s1, e2 - s2) + 1

                if overlap >= shortest * min_overlap:
                    # Merge `f1` into `f2`
                    f2["start"] = min(s1, s2)
                    f2["end"] = max(e1, e2)
                    break
            else:
                # `f1` did not overlap with any fragments
                fragments.append(f1)

        self.fragments = fragments
        return True

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


def merge_supermatches(supermatches: List[Supermatch],
                       min_overlap: float) -> List[Supermatch]:
    merged = []
    for sm in sorted(supermatches):
        for other in merged:
            if other.overlaps(sm, min_overlap):
                break
        else:
            merged.append(sm)

    return merged


def find_overlapping_entries(url: str, src_matches: str,
                             min_overlap: float=0.2, threshold: float=0.75,
                             ora_url: Optional[str]=None):
    logger.info("starting")
    if ora_url:
        con = cx_Oracle.connect(ora_url)
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

        table = Table(con, "INSERT /*+APPEND*/ INTO INTERPRO.SUPERMATCH2 "
                           "VALUES (:1, :2, :3)", autocommit=True)
    else:
        con = table = None

    logger.info("loading MySQL data")
    entries = {}
    signatures = {}
    for acc, info in get_entries(url).items():
        if info["database"] == "interpro":
            entries[acc] = {
                "name": info["root"],
                "type": info["type"],
                "root": info["root"]
            }
        elif info["integrated"]:
            signatures[acc] = info["integrated"]

    counts = {}
    intersections = {}
    cnt = 0
    with io.Store(src_matches) as matches:
        for protein_acc, protein_entries in matches:
            supermatches = []
            for signature_acc, locations in protein_entries.items():
                try:
                    entry_acc = signatures[signature_acc]
                except KeyError:
                    continue
                else:
                    root = entries[entry_acc]["root"]
                    for loc in locations:
                        supermatches.append(Supermatch(entry_acc, root,
                                                       loc["fragments"]))

            # Merge overlapping supermatches
            merged = {}
            for sm in merge_supermatches(supermatches, min_overlap):
                for entry_acc in sm.entries:
                    if table:
                        table.insert((protein_acc, entry_acc, sm.fragments_str))

                    if entry_acc in merged:
                        # TODO: check if that happens
                        print(protein_acc)
                        logger.error(entry_acc)
                        logger.error(merged[entry_acc])
                        logger.error(sm.fragments)
                        raise RuntimeError()
                    else:
                        merged[entry_acc] = sm.fragments

            intersect(merged, counts, intersections)

            cnt += 1
            if not cnt % 10000000:
                logger.info(f"{cnt:>12,}")

        logger.info(f"{cnt:>12,}")

    if table:
        table.close()
        logger.info(f"{table.count} supermatches inserted")

        logger.info("indexing SUPERMATCH2")
        cur = con.cursor()
        cur.execute(
            """
            CREATE UNIQUE INDEX PK_SUPERMATCH2
            ON INTERPRO.SUPERMATCH2 (PROTEIN_AC, ENTRY_AC, POS_FROM, POS_TO, 
                                     DBCODE)
            NOLOGGING
            """
        )
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
        cur.execute("GRANT SELECT ON INTERPRO.SUPERMATCH2 TO INTERPRO_SELECT")
        cur.close()
        con.close()

    # Compute Jaccard coefficients
    overlapping = {}
    supfam = "homologous_superfamily"
    types = ("homologous_superfamily", "domain", "family", "repeat")
    for acc1 in intersections:
        cnt1 = counts[acc1]

        for acc2, (o1, o2) in intersections[acc1].items():
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
                        "accession": acc1,
                        "name": e1["name"],
                        "type": e1["type"]
                    }

                    e2 = {
                        "accession": acc2,
                        "name": e2["name"],
                        "type": e2["type"]
                    }

                    for k, v in [(acc1, e2), (acc2, e1)]:
                        try:
                            overlapping[k].append(v)
                        except KeyError:
                            overlapping[k] = [v]

    logger.info("updating table")
    con = MySQLdb.connect(**parse_url(url), charset="utf8")
    cur = con.cursor()
    for acc, obj in overlapping.items():
        cur.execute(
            """
            UPDATE webfront_entry
            SET overlaps_with = %s
            WHERE accession = %s
            """, (json.dumps(obj), acc)
        )

    con.commit()
    cur.close()
    con.close()

    logger.info("complete")


def intersect(entries: Dict[str, List[Dict]], counts: Dict[str, int],
              intersections: Dict[str, Dict[str, List[int, int]]]):
    for acc1 in entries:
        try:
            counts[acc1] += 1
        except KeyError:
            counts[acc1] = 1

        for acc2 in entries:
            if acc2 >= acc1:
                continue
            elif acc1 in intersections:
                try:
                    overlaps = intersections[acc1][acc2]
                except KeyError:
                    overlaps = intersections[acc1][acc2] = [0, 0]
            else:
                intersections[acc1] = {acc2: [0, 0]}
                overlaps = intersections[acc1] [acc2]

            flag = 0
            for f1 in entries[acc1]:
                len1 = f1["end"] - f1["start"] + 1

                for f2 in entries[acc2]:
                    len2 = f2["end"] - f2["start"] + 1
                    o = min(f1["end"] - f2["end"]) - max(f1["start"] - f2["start"]) + 1

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



def _is_structure_overlapping(start: int, end: int, chains: Dict[str, List[dict]]) -> bool:
    for locations in chains.values():
        for loc in locations:
            if is_overlapping(start, end, loc["protein_start"], loc["protein_end"]):
                return True
    return False


def _export(my_url: str, src_proteins: str, src_proteomes:str,
            src_matches: str, src_ida: str, dst_xrefs: str,
            sync_frequency: int, tmpdir: Optional[str]) -> io.Store:
    # Get required MySQL data
    entries = get_entries(my_url)
    accessions = set(entries)
    entry2set = get_entry2set(my_url)
    protein2structures = {}
    for pdb_id, s in structure.get_structures(my_url).items():
        for protein_acc, chains in s["proteins"].items():
            try:
                protein = protein2structures[protein_acc]
            except KeyError:
                protein = protein2structures[protein_acc] = {}
            finally:
                protein[pdb_id] = chains

    # Open existing stores containing protein-related info
    proteins = io.Store(src_proteins)
    protein2proteome = io.Store(src_proteomes)
    protein2matches = io.Store(src_matches)
    protein2ida = io.Store(src_ida)

    store = io.Store(dst_xrefs, io.Store.chunk_keys(accessions, 10), tmpdir)
    entry_match_counts = {
        "mobidb-lite": 0
    }
    mobidb_proteins = 0
    cnt_proteins = 0
    for protein_acc, p in proteins:
        cnt_proteins += 1
        if not cnt_proteins % 10000000:
            logger.info(f"{cnt_proteins:>12}")

        xrefs = {
            "domain_architectures": set(),
            "proteins": {(protein_acc, p["identifier"])},
            "proteomes": set(),
            "structures": set(),
            "taxa": {p["taxon"]}
        }

        try:
            ida, ida_id = protein2ida[protein_acc]
        except KeyError:
            pass
        else:
            xrefs["domain_architectures"].add(ida)

        try:
            upid = protein2proteome[protein_acc]
        except KeyError:
            pass
        else:
            xrefs["proteomes"].add(upid)

        structures = protein2structures.get(protein_acc, {})

        matches = {}
        match_counts = {}
        for match in protein2matches.get(protein_acc, []):
            method_acc = match["method_ac"]
            try:
                matches[method_acc] += match["fragments"]
            except KeyError:
                matches[method_acc] = list(match["fragments"])
                match_counts[method_acc] = 0
            finally:
                match_counts[method_acc] += 1

            entry_acc = entries[method_acc]["integrated"]
            if entry_acc:
                try:
                    match_counts[entry_acc] += 1
                except KeyError:
                    match_counts[entry_acc] = 1

        """
        As MobiDB-lite is not included in EBISearch,
        we do not need to keep track of matched proteins.
        We only need the number of proteins matched.
        """
        try:
            cnt = match_counts.pop("mobidb-lite")
        except KeyError:
            pass
        else:
            mobidb_proteins += 1
            entry_match_counts["mobidb-lite"] += cnt
            del matches["mobidb-lite"]

        for entry_acc, cnt in match_counts.items():
            try:
                entry_match_counts[entry_acc] += cnt
            except KeyError:
                entry_match_counts[entry_acc] = cnt

        for method_acc, fragments in matches.items():
            _xrefs = xrefs.copy()

            fragments.sort(key=repr_frag)
            start = fragments[0]["start"]
            end = fragments[-1]["end"]

            for pdb_id, chains in structures.items():
                if _is_structure_overlapping(start, end, chains):
                    _xrefs["structures"].add(pdb_id)

            store.update(method_acc, _xrefs)

            entry_acc = entries[method_acc]["integrated"]
            if entry_acc:
                store.update(entry_acc, _xrefs)

        if not cnt_proteins % sync_frequency:
            store.sync()

    proteins.close()
    protein2proteome.close()
    protein2matches.close()
    protein2ida.close()
    logger.info(f"{cnt_proteins:>12}")

    # Add protein count for MobiDB-lite
    store.update("mobidb-lite", {"proteins": mobidb_proteins})

    # Add match counts and set for entries *with* protein matches
    for entry_acc, cnt in entry_match_counts.items():
        xrefs = {"matches": cnt}
        try:
            set_acc = entry2set[entry_acc]
        except KeyError:
            xrefs["sets"] = set()
        else:
            xrefs["sets"] = {set_acc}
        finally:
            store.update(entry_acc, xrefs)
            accessions.remove(entry_acc)

    # Remaining entries without protein matches
    for entry_acc in accessions:
        xrefs = {
            "domain_architectures": set(),
            "matches": 0,
            "proteins": set(),
            "proteomes": set(),
            "structures": set(),
            "taxa": set()
        }

        try:
            set_acc = entry2set[entry_acc]
        except KeyError:
            xrefs["sets"] = set()
        else:
            xrefs["sets"] = {set_acc}
        finally:
            store.update(entry_acc, xrefs)

    return store


def export_xrefs(my_url: str, src_proteins: str, src_proteomes:str,
                 src_matches: str, src_ida: str, dst: str, processes: int=1,
                 sync_frequency: int=100000, tmpdir: Optional[str]=None):
    logger.info("starting")
    if tmpdir:
        os.makedirs(tmpdir, exist_ok=True)

    store = _export(my_url, src_proteins, src_proteomes, src_matches, src_ida,
                    dst, sync_frequency, tmpdir)
    size = store.merge(processes=processes)
    logger.info("disk usage: {:.0f}MB".format(size/1024**2))


def update_counts(my_url: str, src_entries: str, tmpdir: Optional[str]=None):
    logger.info("starting")
    if tmpdir:
        os.makedirs(tmpdir, exist_ok=True)

    entry2set = get_entry2set(my_url)

    con = MySQLdb.connect(**parse_url(my_url), use_unicode=True, charset="utf8")
    cur = con.cursor()
    cur.close()
    with io.KVdb(dir=tmpdir) as kvdb:
        logger.info("updating webfront_entry")
        query = "UPDATE webfront_entry SET counts = %s WHERE accession = %s"
        with Table(con, query) as table, io.Store(src_entries) as store:
            for entry_acc, xrefs in store:
                table.update((json.dumps(reduce(xrefs)), entry_acc))

                if entry_acc in entry2set:
                    kvdb[entry_acc] = xrefs

        logger.info("updating webfront_set")
        query = "UPDATE webfront_set SET counts = %s WHERE accession = %s"
        with Table(con, query) as table:
            for set_acc, s in get_sets(my_url).items():
                set_xrefs = {
                    "domain_architectures": set(),
                    "entries": {
                        s["database"]: len(s["members"]),
                        "total": len(s["members"])
                    },
                    "proteins": set(),
                    "proteomes": set(),
                    "structures": set(),
                    "taxa": set()
                }

                for entry_acc in s["members"]:
                    try:
                        entry_xrefs = kvdb[entry_acc]
                    except KeyError:
                        continue
                    else:
                        for key in entry_xrefs:
                            try:
                                set_xrefs[key] |= entry_xrefs[key]
                            except KeyError:
                                pass

                table.update((json.dumps(reduce(set_xrefs)), set_acc))

        logger.info("disk usage: {:.0f}MB".format(kvdb.size/1024**2))

    con.commit()
    con.close()
    logger.info("complete")


def get_entry2set(my_url: str) -> Dict[str, str]:
    entry2set = {}
    for set_acc, s in get_sets(my_url).items():
        for entry_acc in s["members"]:
            entry2set[entry_acc] = set_acc

    return entry2set
