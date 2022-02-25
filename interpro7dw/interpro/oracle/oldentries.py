import bisect
import hashlib
from datetime import datetime
from typing import Dict, List, Sequence

from interpro7dw import intact, uniprot
from interpro7dw.utils import logger
from interpro7dw.utils.oracle import lob_as_str
from interpro7dw.utils.store import SimpleStore, Store
from interpro7dw.utils.store import dumpobj, loadobj

import cx_Oracle


class Entry:
    def __init__(self, accession: str, short_name: str, name: str,
                 entry_type: str, database: str):
        self.accession = accession
        self.short_name = short_name
        self.name = name
        self.type = entry_type
        self.source_database = database
        self.is_public = True       # Only false for retired InterPro entries
        self.clan = None            # all entries (if any: dict)
        self.counts = {}            # all entries
        self.creation_date = None
        self.deletion_date = None   # Only not None for retired InterPro entries
        self.descriptions = []      # all entries
        self.go_terms = []          # InterPro only
        self.hierarchy = {}         # InterPro only
        self.history = {            # InterPro only
            "names": [],
            "signatures": {}
        }
        self.integrated_in = None   # signatures only
        self.evidence = None        # signatures only
        self.integrates = {}        # InterPro only
        self.literature = {}        # all entries
        self.overlaps_with = []     # InterPro only
        self.pathways = {}          # InterPro only
        self.ppi = []               # prot-prot interactions (InterPro only)
        self.xrefs = {}             # InterPro only

    @property
    def database(self) -> str:
        return self.source_database.lower()

    @property
    def relations(self) -> tuple:
        if not self.hierarchy:
            return None, []

        parent, children = self.traverse_hierarchy(self.hierarchy,
                                                   self.accession)
        return parent, children

    @staticmethod
    def traverse_hierarchy(node: dict, accession: str) -> tuple:
        if node["accession"] == accession:
            return None, [child["accession"] for child in node["children"]]

        for child in node["children"]:
            parent, children = Entry.traverse_hierarchy(child, accession)

            if parent:
                return parent, children
            elif parent is None:
                return node["accession"], children

        return False, []


def _make_hierarchy(accession: str,
                    entries: Dict[str, Entry],
                    child2parent: Dict[str, str],
                    parent2children: Dict[str, Sequence[str]]) -> dict:
    # Find root
    parent_acc = child2parent.get(accession)

    while parent_acc:
        accession = parent_acc
        parent_acc = child2parent.get(accession)

    return _format_node(accession, entries, parent2children)


def _format_node(accession: str,
                 entries: Dict[str, Entry],
                 parent2children: Dict[str, Sequence[str]]) -> dict:
    children = []
    for child_acc in sorted(parent2children.get(accession, [])):
        children.append(_format_node(child_acc, entries, parent2children))

    entry = entries[accession]
    return {
        "accession": accession,
        "name": entry.name,
        "type": entry.type,
        "children": children
    }


def _get_active_interpro_entries(cur: cx_Oracle.Cursor) -> Dict[str, Entry]:
    entries = {}
    cur.execute(
        """
        SELECT
          E.ENTRY_AC, ET.ABBREV, E.NAME, E.SHORT_NAME,
          E.CREATED, E2C.ORDER_IN, CA.TEXT
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

    for interpro_acc, _type, name, short_name, date, descr_id, text in cur:
        try:
            e = entries[interpro_acc]
        except KeyError:
            e = entries[interpro_acc] = Entry(interpro_acc, short_name, name,
                                              _type, "InterPro")
            e.creation_date = date

        if text:
            e.descriptions.append((descr_id, text))

    # Sorts descriptions
    for e in entries.values():
        if not e.descriptions:
            raise ValueError(f"{e.accession}: no descriptions")

        descriptions = []
        for descr_id, text in sorted(e.descriptions):
            descriptions.append(text)

        e.descriptions = descriptions

    return entries


def _add_go_terms(cur: cx_Oracle.Cursor, goa_url: str,
                  entries: Dict[str, Entry]):
    interpro2go = {}
    cur.execute("SELECT ENTRY_AC, GO_ID FROM INTERPRO.INTERPRO2GO")
    for interpro_acc, go_id in cur:
        if interpro_acc not in entries:
            continue

        try:
            interpro2go[interpro_acc].append(go_id)
        except KeyError:
            interpro2go[interpro_acc] = [go_id]

    # Gets GO terms from GOA.
    go_terms = uniprot.goa.get_terms(goa_url)

    while interpro2go:
        interpro_acc, term_ids = interpro2go.popitem()
        terms = []

        for go_id in term_ids:
            try:
                name, aspect, aspect_full, order = go_terms[go_id]
            except KeyError:
                logger.error(f"{interpro_acc}: term {go_id} not found")
                continue

            terms.append((order, go_id, name, aspect, aspect_full))

        # Sort terms
        entry = entries[interpro_acc]
        for _, go_id, name, aspect, aspect_full in sorted(terms):
            entry.go_terms.append({
                "identifier": go_id,
                "name": name,
                "category": {
                    "code": aspect,  # e.g. F
                    "name": aspect_full  # e.g. molecular_function
                }
            })


def _add_hierarchies(cur: cx_Oracle.Cursor, entries: Dict[str, Entry]):
    child2parent = {}
    parent2children = {}
    cur.execute(
        """
        SELECT ENTRY_AC, PARENT_AC
        FROM INTERPRO.ENTRY2ENTRY
        WHERE RELATION = 'TY'
        """
    )
    for child, parent in cur:
        if child not in entries:
            raise KeyError(f"{child}: unchecked entry in hierarchy")
        elif parent not in entries:
            raise KeyError(f"{parent}: unchecked entry in hierarchy")

        child2parent[child] = parent
        try:
            parent2children[parent].append(child)
        except KeyError:
            parent2children[parent] = [child]

    for entry in entries.values():
        entry.hierarchy = _make_hierarchy(entry.accession, entries,
                                          child2parent, parent2children)


def _add_xrefs(cur: cx_Oracle.Cursor, entries: Dict[str, Entry]):
    cur.execute(
        """
        SELECT X.ENTRY_AC, X.AC, LOWER(D.DBSHORT)
        FROM INTERPRO.ENTRY_XREF X
        INNER JOIN INTERPRO.CV_DATABASE D ON X.DBCODE = D.DBCODE
        """
    )
    for interpro_acc, xref_id, xref_db in cur:
        try:
            entry = entries[interpro_acc]
        except KeyError:
            continue

        try:
            entry.xrefs[xref_db].append(xref_id)
        except KeyError:
            entry.xrefs[xref_db] = [xref_id]


def _get_freeze_dates(cur: cx_Oracle.Cursor) -> tuple:
    cur.execute(
        """
        SELECT VERSION, TIMESTAMP
        FROM INTERPRO.DB_VERSION_AUDIT
        WHERE DBCODE = 'I'
        ORDER BY TIMESTAMP
        """
    )
    versions = []
    dates = []
    for version, date in cur:
        versions.append(version)
        dates.append(date)

    return versions, dates


def _get_retired_interpro_entries(cur: cx_Oracle.Cursor) -> List[Entry]:
    """Returns a list of InterPro entries that are not public anymore
    (i.e. deleted of checked=N, including in the upcoming release).

    Only entries that were public at least once are returned.
    We use production freeze times to evaluate if an entry was still existing
    and checked=Y for a release.

    :param cur: Oracle cursor.
    :return: A list of entries.
    """
    versions, dates = _get_freeze_dates(cur)

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
    entries_per_release = {}
    for acc, _type, name, short_name, timestamp, action, checked in cur:
        """
        Example:
        dates: [2021-04-01, 2021-06-01, 2021-08-01]
        timestamp: 2021-04-28
        Then, i will be 1
        and we will associate this edit to 2021-06-01
        """
        i = bisect.bisect_left(dates, timestamp)

        try:
            version = versions[i]
        except IndexError:
            # edit made after the most recent freeze time: ignore
            # (not for this/upcoming release but for the next one)
            continue

        try:
            releases = entries_per_release[acc]
        except KeyError:
            releases = entries_per_release[acc] = {}

        # public in a release if checked and not deleted
        releases[version] = action != "D" and checked == "Y"

        try:
            e = entries[acc]
        except KeyError:
            e = entries[acc] = Entry(acc, short_name, name, _type, "InterPro")
            e.creation_date = timestamp
            e.deletion_date = timestamp
            e.is_public = False
        else:
            # Update entry
            e.short_name = short_name
            e.name = name
            e.type = _type
            e.deletion_date = timestamp

    results = []
    for e in entries.values():
        if any(entries_per_release[e.accession].values()):
            # Entry public for at least once release
            results.append(e)

    return results


def _get_past_names(cur: cx_Oracle.Cursor) -> Dict[str, List[str]]:
    """Returns all the names that each InterPro entry ever had.
    Names are sorted chronologically.

    :param cur: Oracle connection cursor.
    :return: A dictionary (key: entry accession, value: list of names)
    """
    versions, dates = _get_freeze_dates(cur)

    # Gets all names assigned to entries
    cur.execute(
        """
        SELECT ENTRY_AC, TRIM(NAME) AS NAME, TIMESTAMP
        FROM INTERPRO.ENTRY_AUDIT
        WHERE NAME IS NOT NULL
        ORDER BY TIMESTAMP
        """
    )

    entry2names = {}
    for acc, name, timestamp in cur:
        try:
            entry2names[acc].append((name, timestamp))
        except KeyError:
            entry2names[acc] = [(name, timestamp)]

    for acc, names in entry2names.items():
        # Selects the last name given to an entry before each release
        releases = {}
        for name, timestamp in names:
            i = bisect.bisect_left(dates, timestamp)
            try:
                version = versions[i]
            except IndexError:
                # edit made after the most recent freeze time: ignore
                # (not for this/upcoming release but for the next one)
                continue

            if version not in releases or timestamp > releases[version][0]:
                releases[version] = (timestamp, name)

        # Sorts names by oldest to newest
        names = []
        for _, name in sorted(releases.values(), key=lambda x: x[0]):
            if name not in names:
                names.append(name)

        entry2names[acc] = names

    return entry2names


def _get_past_integrations(cur: cx_Oracle.Cursor) -> Dict[str, Dict]:
    """Returns all the signatures that each InterPro entry ever integrated.

    :param cur: Oracle connection cursor.
    :return: A dictionary (key: entry accession, value: dict of mem databases)
    """
    versions, dates = _get_freeze_dates(cur)

    # Gets all signatures that were ever integrated
    cur.execute(
        """
        SELECT ENTRY_AC, METHOD_AC, TIMESTAMP, ACTION
        FROM INTERPRO.ENTRY2METHOD_AUDIT
        ORDER BY TIMESTAMP
        """
    )

    entries = {}
    for interpro_acc, signature_acc, timestamp, action in cur:
        try:
            releases = entries[interpro_acc]
        except KeyError:
            releases = entries[interpro_acc] = {}

        i = bisect.bisect_left(dates, timestamp)
        try:
            version = versions[i]
        except IndexError:
            # edit made after the most recent freeze time: ignore
            # (not for this/upcoming release but for the next one)
            continue

        try:
            signatures = releases[version]
        except KeyError:
            signatures = releases[version] = set()

        if action in ('I', 'U'):
            # Add a signature for Insert/Update actions
            signatures.add(signature_acc)
        elif signature_acc in signatures:  # Delete action
            signatures.remove(signature_acc)

    # Gets signatures, with the entry they are currently integrated in
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

    now_integrated = {}
    for signature_acc, database, interpro_acc in cur:
        now_integrated[signature_acc] = (database, interpro_acc)

    for interpro_acc, releases in entries.items():
        mem_databases = {}

        for signatures in releases.values():
            for signature_acc in signatures:
                try:
                    database, now_interpro_acc = now_integrated[signature_acc]
                except KeyError:
                    database = "deleted"
                    now_interpro_acc = None

                try:
                    mem_databases[database][signature_acc] = now_interpro_acc
                except KeyError:
                    mem_databases[database] = {signature_acc: now_interpro_acc}

        entries[interpro_acc] = mem_databases

    return entries


def _add_signatures(cur: cx_Oracle.Cursor, entries: Dict[str, Entry]):
    cur.execute(
        """
        SELECT
          M.METHOD_AC, M.NAME, M.DESCRIPTION, M.ABSTRACT, M.ABSTRACT_LONG,
          M.METHOD_DATE, ET.ABBREV, DB.DBSHORT, E2M.ENTRY_AC,
          EVI.ABBREV
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
            WHERE CHECKED = 'Y'
          )
        LEFT OUTER JOIN INTERPRO.IPRSCAN2DBCODE I2D
          ON M.DBCODE = I2D.DBCODE
        LEFT OUTER JOIN INTERPRO.CV_EVIDENCE EVI
          ON I2D.EVIDENCE = EVI.CODE
        WHERE M.DBCODE != 'g'  -- discarding MobiDB-Lite
        """
    )

    for rec in cur:
        acc = rec[0]
        short_name = rec[1]
        name = rec[2]
        descr_text = rec[3] or rec[4]
        date = rec[5]
        _type = rec[6]
        database = rec[7]
        interpro_acc = rec[8]
        evidence = rec[9]
        if not evidence:
            raise ValueError(f"{acc}: no evidence")

        e = entries[acc] = Entry(acc, short_name, name, _type, database)
        e.creation_date = date

        if descr_text:
            e.descriptions.append(descr_text)

        e.integrated_in = interpro_acc
        e.evidence = evidence

        if interpro_acc:
            e = entries[interpro_acc]
            try:
                e.integrates[database][acc] = name or short_name
            except KeyError:
                e.integrates[database] = {acc: name or short_name}


def _add_citations(cur: cx_Oracle.Cursor, entries: Dict[str, Entry]):
    citations = {}
    cur.execute(
        """
        SELECT PUB_ID, PUBMED_ID, ISBN, VOLUME, ISSUE,
          YEAR, TITLE, URL, RAWPAGES, MEDLINE_JOURNAL,
          ISO_JOURNAL, AUTHORS, DOI_URL
        FROM INTERPRO.CITATION
        """
    )

    for rec in cur:
        authors = []
        if rec[11]:
            for name in rec[11].split(','):
                authors.append(name.strip())

        citations[rec[0]] = {
            "PMID": rec[1],
            "ISBN": rec[2],
            "volume": rec[3],
            "issue": rec[4],
            "year": rec[5],
            "title": rec[6],
            "URL": rec[7],
            "raw_pages": rec[8],
            "medline_journal": rec[9],
            "ISO_journal": rec[10],
            "authors": authors,
            "DOI_URL": rec[12]
        }

    cur.execute(
        """
        SELECT ENTRY_AC, PUB_ID
        FROM INTERPRO.ENTRY2PUB
        UNION ALL
        SELECT ENTRY_AC, PUB_ID
        FROM INTERPRO.SUPPLEMENTARY_REF
        UNION ALL
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


def export_entries(interpro_url: str, goa_url: str, intact_url: str,
                   clans_file: str, overlapping_file: str, xrefs_file: str,
                   entries_file: str, update: bool = False):
    """Export InterPro entries and member database signatures.

    :param interpro_url: InterPro Oracle connection string.
    :param goa_url: GOA Oracle connection string.
    :param intact_url: IntAct Oracle connection string.
    :param clans_file: Data file of clans and their members.
    :param overlapping_file: Data file of overlapping InterPro entries.
    :param xrefs_file: Store of entry cross-references.
    :param entries_file: Output file.
    :param update: If True, add pathways cross-reference to Oracle.
    """
    logger.info("loading from Oracle databases")

    con = cx_Oracle.connect(interpro_url)
    cur = con.cursor()
    cur.outputtypehandler = lob_as_str  # to fetch CLOB object as strings

    # Starts with InterPro entries
    entries = _get_active_interpro_entries(cur)

    # Adds GO terms
    _add_go_terms(cur, goa_url, entries)

    # Adds entry hierarchies
    _add_hierarchies(cur, entries)

    # Adds cross-references
    _add_xrefs(cur, entries)

    # Adds protein-protein interactions from IntAct
    for interpro_acc, ppi in intact.get_interactions(intact_url).items():
        if interpro_acc in entries:
            entries[interpro_acc].ppi = ppi

    # Adds retired entries (that were at least public in one release)
    for e in _get_retired_interpro_entries(cur):
        entries[e.accession] = e

    # Adds historical names
    for interpro_acc, names in _get_past_names(cur).items():
        if interpro_acc in entries:
            entries[interpro_acc].history["names"] = names

    # Adds historical signatures
    for interpro_acc, member_databases in _get_past_integrations(cur).items():
        if interpro_acc in entries:
            entries[interpro_acc].history["signatures"] = member_databases

    # Adds member database signatures
    _add_signatures(cur, entries)

    # Adds literature references
    _add_citations(cur, entries)

    cur.close()
    con.close()

    # Adds clans on signatures
    logger.info("loading clans")
    for clan in loadobj(clans_file).values():
        for acc, score, seq_length in clan["members"]:
            try:
                entry = entries[acc]
            except KeyError:
                continue
            else:
                entry.clan = {
                    "accession": clan["accession"],
                    "name": clan["name"]
                }

    # Add relationships between overlapping entries
    logger.info("loading overlapping entries")
    for interpro_acc, other_interpro_acc in loadobj(overlapping_file):
        e1 = entries[interpro_acc]
        e2 = entries[other_interpro_acc]

        e1.overlaps_with.append({
            "accession": e2.accession,
            "name": e2.name,
            "type": e2.type
        })

        e2.overlaps_with.append({
            "accession": e1.accession,
            "name": e1.name,
            "type": e1.type
        })

    # Adds cross-references
    logger.info("loading cross-references")
    with SimpleStore(xrefs_file, "r") as store:
        for acc, xrefs in store:
            entry = entries[acc]

            for key in ["metacyc", "reactome"]:
                if xrefs[key]:
                    entry.pathways[key] = list(xrefs[key])

            entry.counts = {
                "domain_architectures": len(xrefs["dom_orgs"]),
                "interactions": len(entry.ppi),
                "matches": xrefs["matches"],
                "pathways": sum([len(v) for v in entry.pathways.values()]),
                "proteins": len(xrefs["proteins"]),
                "proteomes": len(xrefs["proteomes"]),
                "sets": 1 if entry.clan else 0,
                "structural_models": xrefs["struct_models"],
                "structures": len(xrefs["structures"]),
                "taxa": len(xrefs["taxa"]["all"]),
            }

            if xrefs["enzymes"]:
                entry.xrefs["ec"] = sorted(xrefs["enzymes"])

    # Uses default counts for entries without cross-references
    for e in entries.values():
        if not e.counts:
            e.counts = {
                "domain_architectures": 0,
                "interactions": 0,
                "matches": 0,
                "pathways": 0,
                "proteins": 0,
                "proteomes": 0,
                "sets": 0,
                "structural_models": {},
                "structures": 0,
                "taxa": 0,
            }

    logger.info("writing file")
    dumpobj(entries, entries_file)

    if update:
        logger.info("updating ENTRY2PATHWAY")
        con = cx_Oracle.connect(interpro_url)
        cur = con.cursor()
        cur.execute("TRUNCATE TABLE INTERPRO.ENTRY2PATHWAY")
        cur.execute(
            """
            SELECT LOWER(DBSHORT), DBCODE
            FROM INTERPRO.CV_DATABASE
            """
        )
        id2dbcode = dict(cur.fetchall())

        for entry in entries.values():
            for database, pathways in entry.pathways.items():
                dbcode = id2dbcode[database]

                for pathway_id, name in pathways:
                    cur.execute(
                        """
                        INSERT INTO INTERPRO.ENTRY2PATHWAY
                        VALUES (:1, :2, :3, :4)
                        """, (entry.accession, dbcode, pathway_id, name)
                    )

        con.commit()
        cur.close()
        con.close()

    logger.info("done")
