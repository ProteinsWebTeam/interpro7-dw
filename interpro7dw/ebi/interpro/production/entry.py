# -*- coding: utf-8 -*-

import bisect
import hashlib
from multiprocessing import Process, Queue
from typing import List, Mapping, Optional, Sequence

import cx_Oracle

from interpro7dw import logger, metacyc
from interpro7dw.ebi import goa, intact, uniprot
from interpro7dw.ebi.interpro.utils import Table, blob_as_str
from interpro7dw.ebi.interpro.utils import overlaps_pdb_chain, repr_fragment
from interpro7dw.utils import DumpFile, DirectoryTree, Store
from interpro7dw.utils import deepupdate, dumpobj, loadobj, merge_dumps


PATHWAY_DATABASE = {
    "metacyc": 't',
    "reactome": 'r'
}


class Entry:
    def __init__(self, accession: str, sig_type: str, name: str,
                 short_name: str, database: str):
        self.accession = accession
        self.type = sig_type
        self.name = name
        self.short_name = short_name
        self.database = database
        self.integrates = {}
        self.integrated_in = None
        self.go_terms = []
        self.description = []
        self.literature = {}
        self.hierarchy = {}
        self.cross_references = {}
        self.pathways = {}
        self.ppi = []  # protein-protein interactions
        self.is_deleted = False
        self.history = {}
        self.counts = {}
        self.creation_date = None
        self.deletion_date = None
        self.clan = None
        self.overlaps_with = []

    @property
    def relations(self) -> tuple:
        if not self.hierarchy:
            return None, []

        parent, children = self.traverse_hierarchy(self.hierarchy,
                                                   self.accession)
        return parent, children

    def add_contributing_signature(self, accession: str, database: str,
                                   name: str, description: Optional[str]):
        try:
            obj = self.integrates[database]
        except KeyError:
            obj = self.integrates[database] = {}
        finally:
            obj[accession] = description or name

    def set_hierarchy(self, entries: dict, parent_of: dict, children_of: dict):
        # Find root
        accession = self.accession
        parent_acc = parent_of.get(accession)

        while parent_acc:
            accession = parent_acc
            parent_acc = parent_of.get(accession)

        self.hierarchy = Entry.format_node(entries, children_of, accession)

    @staticmethod
    def traverse_hierarchy(node, accession) -> tuple:
        if node["accession"] == accession:
            return None, [child["accession"] for child in node["children"]]

        for child in node["children"]:
            parent, children = Entry.traverse_hierarchy(child, accession)

            if parent:
                return parent, children
            elif parent is None:
                return node["accession"], children

        return False, []

    @staticmethod
    def format_node(entries: dict, children_of: dict, accession: str) -> dict:
        children = []
        for child_acc in sorted(children_of.get(accession, [])):
            children.append(Entry.format_node(entries, children_of, child_acc))

        try:
            e = entries[accession]
        except KeyError:
            logger.warning(f"{accession}")
            return {
                "accession": accession,
                "name": None,
                "type": None,
                "children": children
            }
        else:
            return {
                "accession": accession,
                "name": e.name,
                "type": e.type.lower(),
                "children": children
            }


def _get_interpro_releases(cur: cx_Oracle.Cursor) -> tuple:
    cur.execute(
        """
        SELECT VERSION, FILE_DATE
        FROM INTERPRO.DB_VERSION_AUDIT
        WHERE DBCODE = 'I'
        ORDER BY FILE_DATE
        """
    )
    rel_numbers = []
    rel_dates = []
    for row in cur:
        rel_numbers.append(row[0])
        rel_dates.append(row[1])

    return rel_numbers, rel_dates


def _get_name_history(cur: cx_Oracle.Cursor) -> dict:
    release_numbers, release_dates = _get_interpro_releases(cur)

    cur.execute(
        """
        SELECT ENTRY_AC, TRIM(NAME) AS NAME, TIMESTAMP
        FROM INTERPRO.ENTRY_AUDIT
        WHERE NAME IS NOT NULL
        ORDER BY TIMESTAMP
        """
    )

    # Get all names assigned to each entry
    entries = {}
    for acc, name, timestamp in cur:
        try:
            entries[acc].append((name, timestamp))
        except KeyError:
            entries[acc] = [(name, timestamp)]

    for acc, names in entries.items():
        # Select the last name given to an entry before each release
        releases = {}
        for name, timestamp in names:
            i = bisect.bisect_left(release_dates, timestamp)
            rel = release_numbers[i - 1]

            if rel not in releases or timestamp > releases[rel]["time"]:
                releases[rel] = {
                    "name": name,
                    "time": timestamp
                }

        names = []
        for rel in sorted(releases.values(), key=lambda x: x["time"]):
            last_name = rel["name"]

            if last_name not in names:
                names.append(last_name)

        entries[acc] = names

    return entries


def _get_integration_history(cur: cx_Oracle.Cursor) -> dict:
    release_numbers, release_dates = _get_interpro_releases(cur)

    # Get all past integrations
    cur.execute(
        """
        SELECT ENTRY_AC, METHOD_AC, TIMESTAMP, ACTION
        FROM INTERPRO.ENTRY2METHOD_AUDIT
        ORDER BY TIMESTAMP
        """
    )

    entries = {}
    for entry_acc, signature_acc, timestamp, action in cur:
        try:
            e = entries[entry_acc]
        except KeyError:
            e = entries[entry_acc] = {}

        i = bisect.bisect_left(release_dates, timestamp)
        rel = release_numbers[i - 1]
        try:
            signatures = e[rel]
        except KeyError:
            signatures = e[rel] = set()

        if action in ('I', 'U'):
            # Add a signature for Insert/Update actions
            signatures.add(signature_acc)
        elif signature_acc in signatures:  # Delete action
            signatures.remove(signature_acc)

    # Get signatures and their current integration
    cur.execute(
        """
        SELECT M.METHOD_AC, LOWER(DB.DBSHORT), E.ENTRY_AC
        FROM INTERPRO.METHOD M
        INNER JOIN INTERPRO.CV_DATABASE DB
          ON M.DBCODE = DB.DBCODE
        LEFT OUTER JOIN (
            SELECT EM.METHOD_AC, EM.ENTRY_AC
            FROM INTERPRO.ENTRY2METHOD EM
            INNER JOIN INTERPRO.ENTRY E
            ON EM.ENTRY_AC = E.ENTRY_AC
            AND E.CHECKED = 'Y'
        ) E ON M.METHOD_AC = E.METHOD_AC
        """
    )
    signatures = {}
    for signature_acc, source_database, entry_acc in cur:
        signatures[signature_acc] = (source_database, entry_acc)

    for entry_acc_then, releases in entries.items():
        member_databases = {}

        for release_signatures in releases.values():
            for signature_acc in release_signatures:
                try:
                    database, entry_acc_now = signatures[signature_acc]
                except KeyError:
                    database = "deleted"
                    entry_acc_now = None

                try:
                    member_databases[database][signature_acc] = entry_acc_now
                except KeyError:
                    member_databases[database] = {signature_acc: entry_acc_now}

        entries[entry_acc_then] = member_databases

    return entries


def _get_citations(cur: cx_Oracle.Cursor) -> dict:
    citations = {}
    cur.execute(
        """
        SELECT PUB_ID, PUBMED_ID, ISBN, VOLUME, ISSUE,
          YEAR, TITLE, URL, RAWPAGES, MEDLINE_JOURNAL,
          ISO_JOURNAL, AUTHORS, DOI_URL
        FROM INTERPRO.CITATION
        """
    )
    for row in cur:
        authors = []
        if row[11]:
            for name in row[11].split(','):
                authors.append(name.strip())

        citations[row[0]] = {
            "PMID": row[1],
            "ISBN": row[2],
            "volume": row[3],
            "issue": row[4],
            "year": row[5],
            "title": row[6],
            "URL": row[7],
            "raw_pages": row[8],
            "medline_journal": row[9],
            "ISO_journal": row[10],
            "authors": authors,
            "DOI_URL": row[12]
        }
    return citations


def _get_retired_interpro_entries(cur: cx_Oracle.Cursor) -> List[Entry]:
    release_numbers, release_dates = _get_interpro_releases(cur)

    cur.execute(
        """
        SELECT E.ENTRY_AC, T.ABBREV, E.NAME, E.SHORT_NAME,
          E.TIMESTAMP, E.ACTION, E.CHECKED
        FROM INTERPRO.ENTRY_AUDIT E
        LEFT OUTER JOIN INTERPRO.CV_ENTRY_TYPE T
          ON E.ENTRY_TYPE = T.CODE
        WHERE E.ENTRY_AC NOT IN (
            SELECT ENTRY_AC
            FROM INTERPRO.ENTRY
            WHERE CHECKED = 'Y'
        )
        ORDER BY E.TIMESTAMP
        """
    )

    entries = {}
    releases = {}
    for row in cur:
        accession = row[0]
        entry_type = row[1]
        name = row[2]
        short_name = row[3]
        timestamp = row[4]
        is_deleted = row[5] == 'D'
        is_public = row[6] == 'Y'

        i = bisect.bisect_left(release_dates, timestamp)
        rel = release_numbers[i - 1]

        try:
            e = releases[accession]
        except KeyError:
            e = releases[accession] = {}
        finally:
            # public in a release if checked and not deleted
            e[rel] = not is_deleted and is_public

        try:
            e = entries[accession]
        except KeyError:
            e = Entry(accession, entry_type, name, short_name, "interpro")
            e.creation_date = timestamp
            e.deletion_date = timestamp
            e.is_deleted = True  # No necessarily 'deleted', but not public
            entries[accession] = e
        else:
            e.type = entry_type
            e.name = name
            e.short_name = short_name
            e.deletion_date = timestamp

    public_entries = []
    for e in entries.values():
        if any(releases[e.accession].values()):
            # Entry was public in at least one release
            public_entries.append(e)

    return public_entries


def _get_interpro_entries(cur: cx_Oracle.Cursor) -> List[Entry]:
    cur.execute(
        """
        SELECT
          E.ENTRY_AC, ET.ABBREV, E.NAME, E.SHORT_NAME,
          E.CREATED, CA.TEXT
        FROM INTERPRO.ENTRY E
        INNER JOIN INTERPRO.CV_ENTRY_TYPE ET
          ON E.ENTRY_TYPE = ET.CODE
        LEFT OUTER JOIN INTERPRO.ENTRY2COMMON E2C
          ON E.ENTRY_AC = E2C.ENTRY_AC
        LEFT OUTER JOIN INTERPRO.COMMON_ANNOTATION CA
          ON E2C.ANN_ID = CA.ANN_ID
        WHERE E.CHECKED = 'Y'
        """
    )

    entries = {}
    for row in cur:
        accession = row[0]
        entry_type = row[1]
        name = row[2]
        short_name = row[3]
        creation_date = row[4]
        description = row[5]

        try:
            e = entries[accession]
        except KeyError:
            e = entries[accession] = Entry(accession, entry_type, name,
                                           short_name, "interpro")
            e.creation_date = creation_date

        if description:
            e.description.append(description)

    # Contributing signatures
    cur.execute(
        """
        SELECT EM.ENTRY_AC, M.METHOD_AC, LOWER(DB.DBSHORT), M.NAME,
          M.DESCRIPTION
        FROM INTERPRO.ENTRY2METHOD EM
        INNER JOIN INTERPRO.METHOD M
          ON EM.METHOD_AC = M.METHOD_AC
        INNER JOIN INTERPRO.CV_DATABASE DB
          ON M.DBCODE = DB.DBCODE
        """
    )

    for row in cur:
        accession = row[0]

        try:
            e = entries[accession]
        except KeyError:
            continue
        else:
            e.add_contributing_signature(accession=row[1],
                                         database=row[2],
                                         name=row[3],
                                         description=row[4])

    _entries = {}
    for accession, e in entries.items():
        if e.integrates:
            _entries[accession] = e
        else:
            logger.warning(f"{accession} has no contributing signatures")

    entries = _entries

    # GO terms
    cur.execute(
        """
        SELECT I2G.ENTRY_AC, GT.GO_ID, GT.NAME, GT.CATEGORY
        FROM INTERPRO.INTERPRO2GO I2G
        INNER JOIN GO.TERMS@GOAPRO GT
          ON I2G.GO_ID = GT.GO_ID
        ORDER BY GC.SORT_ORDER, I2G.GO_ID
        """
    )

    for accession, go_id, name, category in cur:
        try:
            e = entries[accession]
        except KeyError:
            continue
        else:
            e.go_terms.append({
                "identifier": go_id,
                "name": name,
                "category": {
                    "code": category,
                    "name": goa.CATEGORIES[category]
                }
            })

    # Hierarchy
    cur.execute(
        """
        SELECT ENTRY_AC, PARENT_AC
        FROM INTERPRO.ENTRY2ENTRY
        WHERE RELATION = 'TY'
        """
    )
    parent_of = {}
    children_of = {}
    for child_acc, parent_acc in cur:
        parent_of[child_acc] = parent_acc

        try:
            children_of[parent_acc].append(child_acc)
        except KeyError:
            children_of[parent_acc] = [child_acc]

    for e in entries.values():
        e.set_hierarchy(entries, parent_of, children_of)

    # Literature references
    citations = _get_citations(cur)
    cur.execute(
        """
        SELECT ENTRY_AC, PUB_ID
        FROM INTERPRO.ENTRY2PUB
        UNION ALL
        SELECT ENTRY_AC, PUB_ID
        FROM INTERPRO.SUPPLEMENTARY_REF
        """
    )

    for accession, pub_id in cur:
        try:
            e = entries[accession]
        except KeyError:
            continue
        else:
            e.literature[pub_id] = citations[pub_id]

    # Cross-references
    cur.execute(
        """
        SELECT X.ENTRY_AC, X.AC, LOWER(D.DBSHORT)
        FROM INTERPRO.ENTRY_XREF X
        INNER JOIN INTERPRO.CV_DATABASE D ON X.DBCODE = D.DBCODE
        """
    )

    for accession, ref_id, ref_db in cur:
        try:
            e = entries[accession]
        except KeyError:
            continue

        try:
            e.cross_references[ref_db].append(ref_id)
        except KeyError:
            e.cross_references[ref_db] = [ref_id]

    return list(entries.values())


def _get_signatures(cur: cx_Oracle.Cursor) -> List[Entry]:
    cur.execute(
        """
        SELECT
          M.METHOD_AC, M.NAME, M.DESCRIPTION, M.ABSTRACT, M.ABSTRACT_LONG,
          M.METHOD_DATE, ET.ABBREV, LOWER(DB.DBSHORT), E2M.ENTRY_AC
        FROM INTERPRO.METHOD M
        INNER JOIN INTERPRO.CV_ENTRY_TYPE ET
          ON M.SIG_TYPE = ET.CODE
        INNER JOIN INTERPRO.CV_DATABASE DB
          ON M.DBCODE = DB.DBCODE
        LEFT OUTER JOIN INTERPRO.ENTRY2METHOD E2M
          ON M.METHOD_AC = E2M.METHOD_AC
          AND E2M.ENTRY_AC IN (
            SELECT ENTRY_AC
            FROM INTERPRO.ENTRY
            WHERE CHECKED='Y'
          )
        WHERE M.DBCODE != 'g'  -- discarding MobiDB-Lite
        """
    )

    entries = {}
    for row in cur:
        accession = row[0]
        short_name = row[1] if row[1] != accession else None
        name = row[2]
        creation_date = row[5]
        entry_type = row[6]
        database = row[7]
        integrated_in = row[8]

        e = Entry(accession, entry_type, name, short_name, database)
        e.creation_date = creation_date
        e.integrated_in = integrated_in

        if row[4]:
            e.description.append(row[4].read().lstrip("<p>").rstrip("</p>"))
        elif row[3]:
            e.description.append(row[3].lstrip("<p>").rstrip("</p>"))

        entries[accession] = e

    # Literature references
    citations = _get_citations(cur)
    cur.execute(
        """
        SELECT METHOD_AC, PUB_ID
        FROM INTERPRO.METHOD2PUB
        """
    )

    for accession, pub_id in cur:
        try:
            e = entries[accession]
        except KeyError:
            continue
        else:
            e.literature[pub_id] = citations[pub_id]

    return list(entries.values())


class EntryXrefs:
    def __init__(self):
        self.ida = set()
        self.matches = 0
        self.proteins = set()
        self.proteomes = set()
        self.structures = set()
        self.taxa = set()

    def asdict(self):
        return {
            "domain_architectures": self.ida,
            "matches": self.matches,
            "proteins": self.proteins,
            "proteomes": self.proteomes,
            "structures": self.structures,
            "taxa": self.taxa
        }


class Supermatch:
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


def _process_proteins(inqueue: Queue, entries: Mapping[str, Entry],
                      min_overlap: bool, dt: DirectoryTree, outqueue: Queue):
    xrefs = {}                  # temporary dict accession->xrefs
    xref_files = []             # files containing xrefs
    entries_with_xrefs = set()  # accession of entries having xrefs
    entry_counts = {}           # number of matches
    entry_intersections = {}    # number of overlapping matches
    interpro2enzyme = {}        # InterPro-ENZYME mapping
    interpro2reactome = {}      # InterPro-Reactome mapping

    ida_file = dt.mktemp()
    with DumpFile(ida_file, compress=True) as ida_df:
        i = 0
        for obj in iter(inqueue.get, None):
            uniprot_acc = obj[0]     # str
            protein_info = obj[1]    # dict
            matches = obj[2]         # dict
            proteome_id = obj[3]     # str or None
            pdb_entries = obj[4]     # dict
            enzymes = obj[5]         # set
            pathways = obj[6]        # set

            supermatches = []
            all_locations = []
            for entry_acc, locations in matches.items():
                entry = entries[entry_acc]
                if entry.database == "interpro":
                    # Adding EC / Reactome mapping

                    if enzymes:
                        try:
                            interpro2enzyme[entry_acc] |= enzymes
                        except KeyError:
                            interpro2enzyme[entry_acc] = enzymes.copy()

                    if pathways:
                        try:
                            interpro2reactome[entry_acc] |= pathways
                        except KeyError:
                            interpro2reactome[entry_acc] = pathways.copy()
                elif entry.database == "pfam":
                    # Storing matches for IDA
                    for loc in locations:
                        all_locations.append({
                            "pfam": entry_acc,
                            "interpro": entry.integrated_in,
                            # We do not consider fragmented locations
                            "start": loc["fragments"][0]["start"],
                            "end": max(f["end"] for f in loc["fragments"])
                        })

                # Adding cross-references (except IDA, still being calculated)
                try:
                    entry_xrefs = xrefs[entry_acc]
                except KeyError:
                    entry_xrefs = xrefs[entry_acc] = EntryXrefs()
                    entries_with_xrefs.add(entry_acc)

                entry_xrefs.matches += len(locations)
                entry_xrefs.proteins.add((
                    uniprot_acc,
                    protein_info["identifier"]
                ))

                if proteome_id:
                    entry_xrefs.proteomes.add(proteome_id)

                for pdb_id, chains in pdb_entries.items():
                    for chain_id, segments in chains.items():
                        if overlaps_pdb_chain(locations, segments):
                            entry_xrefs.structures.add(pdb_id)
                            break  # Skip other chains

                entry_xrefs.taxa.add(protein_info["taxid"])

                # Create a Supermatch for each integrated signature match
                if entry.integrated_in:
                    # Integrated member database signature
                    interpro_acc = entry.integrated_in
                    root = entries[interpro_acc].hierarchy["accession"]
                    for loc in locations:
                        sm = Supermatch(interpro_acc, loc["fragments"], root)
                        supermatches.append(sm)

            # Finishing IDA
            domains = []
            dom_members = set()
            for loc in sorted(all_locations, key=repr_fragment):
                if loc["interpro"]:
                    domains.append(f"{loc['pfam']}:{loc['interpro']}")
                    dom_members.add(loc["interpro"])
                else:
                    domains.append(loc["pfam"])

                dom_members.add(loc["pfam"])

            if domains:
                # Flush IDA
                dom_str = '-'.join(domains)
                dom_id = hashlib.sha1(dom_str.encode("utf-8")).hexdigest()
                ida_df.dump((uniprot_acc, dom_members, dom_str, dom_id))

                # Adding cross-references now
                for key in dom_members:
                    xrefs[key].ida.add(dom_id)

            # Merging overlapping supermatches
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

            # Group by entry
            merged_grouped = {}
            for sm in merged:
                for interpro_acc in sm.entries:
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
                        obj = entry_intersections[interpro_acc] = {}

                    try:
                        overlaps = obj[other_acc]
                    except KeyError:
                        """
                        Use a dict rather than a list (or tuple)
                        because deepupdate() would concatenate the lists
                        created by different workers
                        """
                        overlaps = obj[other_acc] = {
                            "1": 0,
                            "2": 0,
                        }

                    flag = 0
                    for f1 in fragments1:
                        start1 = f1["start"]
                        end1 = f1["end"]
                        length1 = end1 - start1 + 1

                        for f2 in fragments2:
                            start2 = f2["start"]
                            end2 = f2["end"]
                            length2 = end2 - start2 + 1
                            overlap = min(end1, end2) - max(start1, start2) + 1

                            if not flag & 1 and overlap >= length1 * 0.5:
                                # 1st time fragments overlap >= 50% of f1
                                flag |= 1
                                overlaps["1"] += 1

                            if not flag & 2 and overlap >= length2 * 0.5:
                                # 1st time fragments overlap >= 50% of f2
                                flag |= 2
                                overlaps["2"] += 1

                        if flag == 3:
                            """
                            Both cases already happened
                              -> no need to keep iterating
                            """
                            break

            i += 1
            if not i % 100000:
                # Flush Xrefs
                file = dt.mktemp()
                with DumpFile(file, compress=True) as xref_df:
                    for entry_acc in sorted(xrefs):
                        xref_df.dump((entry_acc, xrefs[entry_acc].asdict()))

                xrefs = {}
                xref_files.append(file)

    # Remaining xrefs
    file = dt.mktemp()
    with DumpFile(file, compress=True) as df:
        for entry_acc in sorted(xrefs):
            df.dump((entry_acc, xrefs[entry_acc].asdict()))

    xref_files.append(file)

    # Merge files (each worker will produce one merged file)
    xref_file = dt.mktemp()
    with DumpFile(xref_file, compress=True) as df:
        for entry_acc, xrefs in merge_dumps(xref_files):
            df.dump((entry_acc, xrefs))

    outqueue.put((
        xref_file,
        entries_with_xrefs,
        ida_file,
        entry_counts,
        entry_intersections,
        interpro2enzyme,
        interpro2reactome
    ))


def export_entries(url: str, p_metacyc: str, p_clans: str,
                   p_proteins: str, p_structures: str,
                   p_uniprot2matches: str, p_uniprot2proteome: str,
                   p_uniprot2ida: str, p_entry2xrefs: str, p_entries: str,
                   **kwargs):
    min_overlap = kwargs.get("overlap", 0.2)
    processes = kwargs.get("processes", 1)
    min_similarity = kwargs.get("similarity", 0.75)
    tmpdir = kwargs.get("tmpdir")

    con = cx_Oracle.connect(url)
    cur = con.cursor()

    entries = {}
    logger.info("loading active InterPro entries")
    for entry in _get_interpro_entries(cur):
        entries[entry.accession] = entry

    logger.info("enriching entries with IntAct data")
    for accession, interactions in intact.get_interactions(cur).items():
        try:
            entry = entries[accession]
        except KeyError:
            continue
        else:
            entry.ppi = interactions

    logger.info("loading deleted InterPro entries")
    for entry in _get_retired_interpro_entries(cur):
        if entry.accession in entries:
            cur.close()
            con.close()
            raise RuntimeError(f"entry cannot be active "
                               f"and deleted {entry.accession}")

        entries[entry.accession] = entry

    logger.info("loading member database signatures")
    for entry in _get_signatures(cur):
        if entry.integrated_in and entry.integrated_in not in entries:
            cur.close()
            con.close()
            raise RuntimeError(f"{entry.accession} integrated "
                               f"in missing entry ({entry.integrated_in})")

        entries[entry.accession] = entry

    logger.info("loading past entry names")
    past_names = _get_name_history(cur)

    logger.info("loading past signature integrations")
    past_integrations = _get_integration_history(cur)

    logger.info("loading ENZYME")
    u2enzyme = uniprot.get_swissprot2enzyme(cur)

    logger.info("loading Reactome pathways")
    u2reactome = uniprot.get_swissprot2reactome(cur)
    cur.close()
    con.close()

    logger.info("loading MetaCyc pathways")
    ec2metacyc = metacyc.get_ec2pathways(p_metacyc)

    # Updating entry history
    for entry in entries.values():
        try:
            names = past_names[entry.accession]
        except KeyError:
            pass
        else:
            entry.history["names"] = names

        try:
            signatures = past_integrations[entry.accession]
        except KeyError:
            pass
        else:
            entry.history["signatures"] = signatures

    # Updating entry clan info
    for clan in loadobj(p_clans).values():
        for entry_acc, score, seq_length in clan["members"]:
            try:
                entry = entries[entry_acc]
            except:
                continue
            else:
                entry.clan = {
                    "accession": clan["accession"],
                    "name": clan["name"]
                }

    inqueue = Queue(maxsize=processes)
    outqueue = Queue()
    workers = []
    for _ in range(max(1, processes - 1)):
        dt = DirectoryTree(tmpdir)
        p = Process(target=_process_proteins,
                    args=(inqueue, entries, min_overlap, dt, outqueue))
        p.start()
        workers.append((p, dt))

    logger.info("processing")
    uniprot2pdbe = {}
    for pdb_id, entry in loadobj(p_structures).items():
        for uniprot_acc, chains in entry["proteins"].items():
            try:
                uniprot2pdbe[uniprot_acc][pdb_id] = chains
            except KeyError:
                uniprot2pdbe[uniprot_acc] = {pdb_id: chains}

    proteins = Store(p_proteins)
    u2matches = Store(p_uniprot2matches)
    u2proteome = Store(p_uniprot2proteome)
    i = 0
    for uniprot_acc, matches in u2matches.items():
        inqueue.put((
            uniprot_acc,
            proteins[uniprot_acc],
            matches,
            u2proteome.get(uniprot_acc),
            uniprot2pdbe.get(uniprot_acc, {}),
            set(u2enzyme.get(uniprot_acc, [])),
            set(u2reactome.get(uniprot_acc, []))
        ))

        i += 1
        if not i % 10000000:
            logger.info(f"{i:>15,}")

    proteins.close()
    u2matches.close()
    u2proteome.close()
    logger.info(f"{i:>15,}")

    # Send sentinel
    for _ in workers:
        inqueue.put(None)

    # Merge results from workers
    logger.info("exporting domain architectures")
    entries_with_xrefs = set()
    xref_files = []
    entry_counts = {}
    entry_intersections = {}
    interpro2enzyme = {}
    interpro2reactome = {}
    with Store(p_uniprot2ida, u2matches.get_keys(), tmpdir) as u2ida:
        for _ in workers:
            obj = outqueue.get()
            xref_files.append(obj[0])                               # str
            entries_with_xrefs |= obj[1]                            # set
            ida_file = obj[2]                                       # str
            deepupdate(obj[3], entry_counts, replace=False)         # dict
            deepupdate(obj[4], entry_intersections, replace=False)  # dict
            deepupdate(obj[5], interpro2enzyme)                     # dict
            deepupdate(obj[6], interpro2reactome)                   # dict

            with DumpFile(ida_file) as df:
                i = 0
                for uniprot_acc, dom_members, dom_str, dom_id in df:
                    u2ida[uniprot_acc] = (
                        dom_members,
                        dom_str,
                        dom_id
                    )
                    i += 1

                    if not i % 1000:
                        u2ida.sync()

            u2ida.sync()

        size = u2ida.merge(processes=processes)

    # Adding empty EntryXrefs objects for entries without xrefs
    xref_files.append(workers[0][1].mktemp())
    with DumpFile(xref_files[-1], compress=True) as df:
        for entry_acc in sorted(set(entries.keys()) - entries_with_xrefs):
            df.dump((entry_acc, EntryXrefs().asdict()))

    logger.info("exporting cross-references")
    with DumpFile(p_entry2xrefs, compress=True) as df:
        for entry_acc, xrefs in merge_dumps(xref_files):
            df.dump((entry_acc, xrefs))

            entry = entries[entry_acc]

            # Reactome pathways
            if entry_acc in interpro2reactome:
                pathways = interpro2reactome[entry_acc]
                entry.pathways["reactome"] = [
                    dict(zip(("id", "name"), pthw))
                    for pthw in sorted(pathways)
                ]

            # EC numbers
            if entry_acc in interpro2enzyme:
                ecnos = sorted(interpro2enzyme[entry_acc])
                entry.cross_references["ec"] = ecnos

                # MetaCyc pathways
                pathways = set()
                for ecno in ecnos:
                    pathways |= set(ec2metacyc.get(ecno, []))

                if pathways:
                    entry.pathways["metacyc"] = [
                        dict(zip(("id", "name"), pthw))
                        for pthw in sorted(pathways)
                    ]

    for p, dt in workers:
        size += dt.size
        dt.remove()

    logger.info(f"temporary files: {size / 1024 / 1024:.0f} MB")

    logger.info("calculating overlapping relationships")
    supfam = "homologous_superfamily"
    types = (supfam, "domain", "family", "repeat")
    for entry_acc, overlaps in entry_intersections.items():
        entry1 = entries[entry_acc]
        entry_cnt = entry_counts[entry_acc]
        type1 = entry1.type.lower()

        for other_acc, overlap_counts in overlaps.items():
            o1 = overlap_counts["1"]
            o2 = overlap_counts["2"]
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
            type2 = entry2.type.lower()
            if ((type1 == supfam and type2 in types)
                    or (type1 in types and type2 == supfam)):
                # e1 -> e2 relationship
                entry1.overlaps_with.append({
                    "accession": other_acc,
                    "name": entry2.name,
                    "type": type2
                })

                # e2 -> e1 relationship
                entry2.overlaps_with.append({
                    "accession": entry_acc,
                    "name": entry1.name,
                    "type": type1
                })

    dumpobj(p_entries, entries)

    logger.info("populating ENTRY2PATHWAY")
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    cur.execute("TRUNCATE TABLE INTERPRO.ENTRY2PATHWAY")
    cur.close()
    sql = "INSERT INTO INTERPRO.ENTRY2PATHWAY VALUES (:1, :2, :3, :4)"
    with Table(con, sql) as table:
        for e in entries.values():
            for database, pathways in e.pathways.items():
                code = PATHWAY_DATABASE[database]
                for pthw in pathways:
                    table.insert((
                        e.accession,
                        code,
                        pthw["id"],
                        pthw["name"]
                    ))

    con.commit()
    con.close()
    logger.info("complete")


def get_features(cur: cx_Oracle.Cursor) -> dict:
    cur.execute(
        """
        SELECT M.METHOD_AC, M.NAME, D.DBSHORT, V.VERSION, EVI.ABBREV
        FROM INTERPRO.FEATURE_METHOD M
        INNER JOIN INTERPRO.CV_DATABASE D
          ON M.DBCODE = D.DBCODE
        INNER JOIN INTERPRO.DB_VERSION V
          ON D.DBCODE = V.DBCODE
        INNER JOIN INTERPRO.IPRSCAN2DBCODE I2D
          ON D.DBCODE = I2D.DBCODE
        INNER JOIN INTERPRO.CV_EVIDENCE EVI
          ON I2D.EVIDENCE = EVI.CODE
        """
    )
    features = {}
    for row in cur:
        features[row[0]] = {
            "accession": row[0],
            "name": row[1],
            "database": row[2],
            "version": row[3],
            "evidence": row[4]
        }
    return features


def get_hmms(url: str, multi_models: bool = True):
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    cur.outputtypehandler = blob_as_str

    if multi_models:
        sql = """
            SELECT METHOD_AC, MODEL_AC, HMM
            FROM INTERPRO.METHOD_HMM
        """
    else:
        # Ignore databases with signatures having more than one model
        sql = """
            SELECT METHOD_AC, MODEL_AC, HMM
            FROM INTERPRO.METHOD_HMM
            WHERE METHOD_AC NOT IN (
                SELECT METHOD_AC
                FROM INTERPRO.METHOD
                WHERE DBCODE IN (
                    SELECT DISTINCT DBCODE
                    FROM INTERPRO.METHOD
                    WHERE METHOD_AC IN (
                        SELECT METHOD_AC
                        FROM INTERPRO.METHOD_HMM
                        GROUP BY METHOD_AC
                        HAVING COUNT(*) > 1
                    )
                )
            )
        """

    cur.execute(sql)
    for signature_acc, model_acc, hmm_bytes in cur:
        yield signature_acc, model_acc, hmm_bytes

    cur.close()
    con.close()


def get_signatures(cur: cx_Oracle.Cursor) -> dict:
    cur.execute(
        """
        SELECT M.METHOD_AC, M.NAME, DB.DBSHORT, EVI.ABBREV,
               E2M.ENTRY_AC, E2M.NAME, E2M.ABBREV, E2M.PARENT_AC
        FROM INTERPRO.METHOD M
        INNER JOIN  INTERPRO.CV_DATABASE DB
          ON M.DBCODE = DB.DBCODE
        INNER JOIN  INTERPRO.IPRSCAN2DBCODE I2D
          ON M.DBCODE = I2D.DBCODE
        INNER JOIN INTERPRO.CV_EVIDENCE EVI
          ON I2D.EVIDENCE = EVI.CODE
        LEFT OUTER JOIN (
          SELECT E2M.METHOD_AC, E.ENTRY_AC, E.NAME, ET.ABBREV, E2E.PARENT_AC
          FROM INTERPRO.ENTRY E
          INNER JOIN INTERPRO.ENTRY2METHOD E2M
            ON E.ENTRY_AC = E2M.ENTRY_AC
          INNER JOIN INTERPRO.CV_ENTRY_TYPE ET
            ON E.ENTRY_TYPE = ET.CODE
          LEFT OUTER JOIN INTERPRO.ENTRY2ENTRY E2E
            ON E.ENTRY_AC = E2E.ENTRY_AC
          WHERE E.CHECKED = 'Y'
        ) E2M
          ON M.METHOD_AC = E2M.METHOD_AC
        """
    )
    signatures = {}
    for row in cur:
        signatures[row[0]] = {
            "accession": row[0],
            "name": row[1] or row[0],
            "database": row[2],
            "evidence": row[3],
            "interpro": [
                ("id", row[4]),
                ("name", row[5]),
                ("type", row[6]),
                ("parent_id", row[7])
            ] if row[4] else None
        }
    return signatures
