# -*- coding: utf-8 -*-

import bisect
import json
from typing import Dict, List, Optional, Sequence

import cx_Oracle

from interpro7dw import logger
from interpro7dw.utils import Store, datadump
from .utils import DC_STATUSES, condense_locations, repr_fragment


ENTRY_DATABASES = [
    'B',    # SFLD
    'F',    # PRINTS
    'H',    # Pfam
    'I',    # InterPro
    'J',    # CDD
    'M',    # PROSITE profiles
    'N',    # TIGRFAMs
    'P',    # PROSITE patterns
    'Q',    # HAMAP
    'R',    # SMART
    'U',    # PIRSF
    'V',    # PANTHER
    'X',    # CATH-Gene3D
    'Y',    # SUPERFAMILY
    'g',    # MobiDB Lite
]


def chunk_proteins(url: str, output: str, chunk_size: int=50000):
    logger.info("loading")
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    cur.execute(
        """
        SELECT PROTEIN_AC
        FROM INTERPRO.PROTEIN
        """
    )

    accessions = [acc for acc, in cur]
    cur.close()
    con.close()

    logger.info("sorting")
    accessions.sort()

    keys = []
    for i in range(0, len(accessions), chunk_size):
        keys.append(accessions[i])

    Store.dump_keys(keys, output)
    logger.info("complete")


def export_features(url: str, input: str, output: str,
                    dir: Optional[str]=None, processes: int=1):
    logger.info("starting")
    with Store(output, input, dir) as store:
        con = cx_Oracle.connect(url)
        cur = con.cursor()
        cur.execute(
            """
            SELECT FM.PROTEIN_AC, FM.METHOD_AC, LOWER(DB.DBSHORT),
                   FM.POS_FROM, FM.POS_TO, FM.SEQ_FEATURE
            FROM INTERPRO.FEATURE_MATCH FM
            INNER JOIN INTERPRO.CV_DATABASE DB ON FM.DBCODE = DB.DBCODE
            """
        )

        i = 0
        for row in cur:
            protein_acc = row[0]
            signature_acc = row[1]
            database = row[2]
            pos_start = row[3]
            pos_end = row[4]
            seq_feature = row[5]

            if database == "mobidblt" and seq_feature is None:
                seq_feature = "Consensus Disorder Prediction"

            store.append(protein_acc, (signature_acc, database, pos_start,
                                       pos_end, seq_feature))

            i += 1
            if not i % 1000000:
                store.sync()

                if not i % 100000000:
                    logger.info(f"{i:>12,}")

        cur.close()
        con.close()

        logger.info(f"{i:>12,}")
        size = store.merge(fn=_post_features, processes=processes)
        logger.info(f"temporary files: {size/1024/1024:.0f} MB")


def _post_features(matches: Sequence[dict]) -> dict:
    entries = {}
    for acc, database, pos_start, pos_end, seq_feature in matches:
        try:
            obj = entries[acc]
        except KeyError:
            obj = entries[acc] = {
                "accession": acc,
                "source_database": database,
                "locations": []
            }
        finally:
            obj["locations"].append({
                "fragments": [{
                    "start": pos_start,
                    "end": pos_end,
                    "seq_feature": seq_feature
                }]
            })

    for obj in entries.values():
        # Only one fragment per location
        obj["locations"].sort(key=lambda l: repr_fragment(l["fragments"][0]))

    return entries


def export_matches(url: str, input: str, output: str,
                   dir: Optional[str]=None, processes: int=1):
    logger.info("starting")
    with Store(output, input, dir) as store:
        con = cx_Oracle.connect(url)
        cur = con.cursor()
        cur.execute(
            """
            SELECT EM.METHOD_AC, EM.ENTRY_AC 
            FROM INTERPRO.ENTRY2METHOD EM
            INNER JOIN INTERPRO.ENTRY E ON EM.ENTRY_AC = E.ENTRY_AC
            WHERE E.CHECKED = 'Y'
            """
        )
        integrated = dict(cur.fetchall())

        cur.execute(
            """
            SELECT PROTEIN_AC, METHOD_AC, MODEL_AC, POS_FROM, POS_TO,
                   FRAGMENTS
            FROM INTERPRO.MATCH
            """
        )

        i = 0
        for row in cur:
            protein_acc = row[0]
            signature_acc = row[1]
            model_acc = row[2] if signature_acc != row[2] else None
            pos_start = row[3]
            pos_end = row[4]
            fragments_str = row[5]

            if fragments_str is None:
                fragments = [{
                    "start": pos_start,
                    "end": pos_end,
                    "dc-status": DC_STATUSES['S']  # Continuous
                }]
            else:
                fragments = []
                for frag in fragments_str.split(','):
                    # Format: START-END-STATUS
                    s, e, t = frag.split('-')
                    fragments.append({
                        "start": int(s),
                        "end": int(e),
                        "dc-status": DC_STATUSES[t]
                    })
                fragments.sort(key=repr_fragment)

            store.append(protein_acc, {
                "accession": signature_acc,
                "condense": integrated.get(signature_acc),
                "fragments": fragments,
                "model": model_acc
            })

            i += 1
            if not i % 1000000:
                store.sync()

                if not i % 100000000:
                    logger.info(f"{i:>12,}")

        cur.close()
        con.close()

        logger.info(f"{i:>12,}")
        size = store.merge(fn=_post_matches, processes=processes)
        logger.info(f"temporary files: {size/1024/1024:.0f} MB")


def _post_matches(matches: Sequence[dict]) -> dict:
    entries = {}
    signatures = {}
    for match in matches:
        acc = match["accession"]

        try:
            signatures[acc].append({
                "fragments": match["fragments"],
                "model_acc": match["model"]
            })
        except KeyError:
            signatures[acc] = [{
                "fragments": match["fragments"],
                "model_acc": match["model"]
            }]

        acc = match["condense"]
        if acc:
            try:
                entries[acc].append(match["fragments"])
            except KeyError:
                entries[acc] = [match["fragments"]]

    for acc, locations in entries.items():
        condensed = []
        for start, end in condense_locations(locations):
            condensed.append({
                "fragments": [{
                    "start": start,
                    "end": end,
                    "dc-status": DC_STATUSES['S']
                }],
                "model_acc": None
            })

        entries[acc] = condensed

    for acc, locations in signatures.items():
        # Sort locations using their leftmost fragment (fragments are sorted)
        locations.sort(key=lambda l: repr_fragment(l["fragments"][0]))
        entries[acc] = locations

    return entries


def export_proteins(url: str, input: str, output: str,
                    dir: Optional[str]=None, processes: int=1):
    logger.info("starting")
    with Store(output, input, dir) as store:
        con = cx_Oracle.connect(url)
        cur = con.cursor()
        cur.execute(
            """
            SELECT 
              P.PROTEIN_AC, P.NAME, P.DBCODE, P.LEN, P.FRAGMENT, 
              TO_CHAR(P.TAX_ID)
            FROM INTERPRO.PROTEIN P
            INNER JOIN INTERPRO.ETAXI E ON P.TAX_ID = E.TAX_ID
            """
        )

        i = 0
        for row in cur:
            store[row[0]] = {
                "identifier": row[1],
                "reviewed": row[2] == 'S',
                "length": row[3],
                "fragment": row[4] == 'Y',
                "taxid": row[5]
            }

            i += 1
            if not i % 1000000:
                store.sync()

                if not i % 10000000:
                    logger.info(f"{i:>12,}")

        cur.close()
        con.close()

        logger.info(f"{i:>12,}")
        size = store.merge(processes=processes)
        logger.info(f"temporary files: {size/1024/1024:.0f} MB")


def export_residues(url: str, input: str, output: str,
                    dir: Optional[str]=None, processes: int=1):
    logger.info("starting")
    with Store(output, input, dir) as store:
        con = cx_Oracle.connect(url)
        cur = con.cursor()

        cur.execute(
            """
            SELECT S.PROTEIN_AC, S.METHOD_AC, M.NAME, LOWER(D.DBSHORT),
                   S.DESCRIPTION, S.RESIDUE, S.RESIDUE_START, S.RESIDUE_END
            FROM INTERPRO.SITE_MATCH S
            INNER JOIN INTERPRO.METHOD M ON S.METHOD_AC = M.METHOD_AC
            INNER JOIN INTERPRO.CV_DATABASE D ON M.DBCODE = D.DBCODE
            """
        )

        i = 0
        for row in cur:
            protein_acc = row[0]
            signature_acc = row[1]
            signature_name = row[2]
            database = row[3]
            description = row[4]
            residue = row[5]
            pos_start = row[6]
            pos_end = row[7]

            store.append(protein_acc, (signature_acc, signature_name,
                                       database, description, residue,
                                       pos_start, pos_end))

            i += 1
            if not i % 1000000:
                store.sync()

                if not i % 100000000:
                    logger.info(f"{i:>12,}")

        cur.close()
        con.close()

        logger.info(f"{i:>12,}")
        size = store.merge(fn=_post_residues, processes=processes)
        logger.info(f"temporary files: {size/1024/1024:.0f} MB")


def _post_residues(matches: Sequence[dict]) -> dict:
    entries = {}
    for acc, name, database, descr, residue, pos_start, pos_end in matches:
        try:
            obj = entries[acc]
        except KeyError:
            obj = entries[acc] = {
                "accession": acc,
                "name": name,
                "source_database": database,
                "locations": {}
            }

        try:
            d = obj["locations"][descr]
        except KeyError:
            d = obj["locations"][descr] = []
        finally:
            d.append({
                "residues": residue,
                "start": pos_start,
                "end": pos_end
            })

    return entries


def export_sequences(url: str, input: str, output: str,
                     dir: Optional[str]=None, processes: int=1):
    logger.info("starting")
    with Store(output, input, dir) as store:
        con = cx_Oracle.connect(url)
        cur = con.cursor()
        cur.execute(
            """
            SELECT /*+ PARALLEL */ UX.AC, UP.SEQ_SHORT, UP.SEQ_LONG
            FROM UNIPARC.XREF UX
            INNER JOIN UNIPARC.PROTEIN UP ON UX.UPI = UP.UPI
            WHERE UX.DBID IN (2, 3)
            AND UX.DELETED = 'N'
            """
        )

        i = 0
        for row in cur:
            store[row[0]] = row[2].read() if row[2] is not None else row[1]

            i += 1
            if not i % 1000000:
                store.sync()

                if not i % 10000000:
                    logger.info(f"{i:>12,}")

        cur.close()
        con.close()

        logger.info(f"{i:>12,}")
        size = store.merge(processes=processes)
        logger.info(f"temporary files: {size/1024/1024:.0f} MB")


def get_databases(url: str) -> List[tuple]:
    con = cx_Oracle.connect(url)
    cur = con.cursor()

    """
    Using RN=2 to join with the second most recent action
    (the most recent is the current record)
    """
    cur.execute(
        """
        SELECT
          DB.DBCODE, LOWER(DB.DBSHORT), DB.DBNAME, DB.DESCRIPTION,
          V.VERSION, V.FILE_DATE, VA.VERSION, VA.FILE_DATE
        FROM INTERPRO.CV_DATABASE DB
        LEFT OUTER JOIN INTERPRO.DB_VERSION V ON DB.DBCODE = V.DBCODE
        LEFT OUTER JOIN (
          SELECT
            DBCODE, VERSION, FILE_DATE,
            ROW_NUMBER() OVER (
              PARTITION BY DBCODE ORDER BY TIMESTAMP DESC
            ) RN
          FROM INTERPRO.DB_VERSION_AUDIT
          WHERE ACTION = 'U'
        ) VA ON DB.DBCODE = VA.DBCODE AND VA.RN = 2
        """
    )

    databases = []
    for row in cur:
        code = row[0]
        name_short = row[1]
        name_long = row[2]
        description = row[3]
        release_version = row[4]
        release_date = row[5]
        previous_version = row[6]
        previous_date = row[7]

        if code in ENTRY_DATABASES:
            db_type = "entry"
        elif code in ('S', 'T', 'u'):
            if code == 'S':
                name_short = "reviewed"
            elif code == 'T':
                name_short = "unreviewed"

            db_type = "protein"
        else:
            db_type = "other"

        databases.append((
            name_short,
            name_long,
            description,
            db_type,
            release_version,
            release_date,
            previous_version,
            previous_date
        ))

    cur.close()
    con.close()

    return databases


def _get_name_history(cur: cx_Oracle.Cursor) -> Dict:
    cur.execute(
        """
        SELECT VERSION, FILE_DATE
        FROM INTERPRO.DB_VERSION_AUDIT
        WHERE DBCODE = 'I'
        ORDER BY FILE_DATE
        """
    )
    releases = [row[1] for row in cur]

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
        last_changes = {}
        for name, timestamp in names:
            i = bisect.bisect_left(releases, timestamp)

            try:
                # Will raise an IndexError if timestamp > most recent release
                date = releases[i]
            except IndexError:
                continue

            if date not in last_changes or timestamp > last_changes[date]["time"]:
                last_changes[date] = {"name": name, "time": timestamp}

        names = []
        for date in sorted(last_changes):
            last_name = last_changes[date]["name"]

            if last_name not in names:
                names.append(last_name)

        entries[acc] = names

    return entries


def _get_integration_history(cur: cx_Oracle.Cursor) -> Dict:
    cur.execute(
        """
        SELECT VERSION, FILE_DATE
        FROM INTERPRO.DB_VERSION_AUDIT
        WHERE DBCODE = 'I'
        ORDER BY FILE_DATE
        """
    )
    releases = [row[1] for row in cur]

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
        i = bisect.bisect_left(releases, timestamp)

        try:
            # Will raise an IndexError if timestamp > most recent release
            date = releases[i]
        except IndexError:
            continue

        try:
            e = entries[entry_acc]
        except KeyError:
            e = entries[entry_acc] = {}

        try:
            signatures = e[date]
        except KeyError:
            signatures = e[date] = set()

        if action in ('I', 'U'):
            # Add a signature for Insert/Update actions
            signatures.add(signature_acc)
        elif signature_acc in signatures:  # Delete action
            signatures.remove(signature_acc)

    # Get signatures and their current integration
    cur.execute(
        """
        SELECT M.METHOD_AC, LOWER(DB.DBSHORT), EM.ENTRY_AC
        FROM INTERPRO.METHOD M
        INNER JOIN INTERPRO.CV_DATABASE DB
          ON M.DBCODE = DB.DBCODE
        LEFT OUTER JOIN INTERPRO.ENTRY2METHOD EM
          ON M.METHOD_AC = EM.METHOD_AC
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


class Entry(object):
    def __init__(self, accession: str, type: str, name: str, short_name: str,
                 database: str):
        self.accession = accession
        self.type = type
        self.name = name
        self.short_name = short_name
        self.database = database
        self.integrates = {}
        self.integrated_in = None
        self.go_terms = []
        self.description = []
        self.wikipedia = {}
        self.literature = {}
        self.hierarchy = {}
        self.cross_references = {}
        self.ppi = []  # protein-protein interactions
        self.pathways = {}
        self.overlaps_with = []
        self.is_featured = False
        self.is_alive = True
        self.is_public = False
        self.history = {}
        self.counts = {}
        self.creation_date = None
        self.deletion_date = None

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

        self.history = Entry.format_node(entries, children_of, accession)

    @staticmethod
    def format_node(entries: dict, children_of: dict, accession: str) -> Dict:
        children = []
        for child_acc in children_of.get(accession, []):
            children.append(Entry.format_node(entries, children_of, child_acc))

        e = entries[accession]
        return {
            "accession": accession,
            "name": e.name,
            "type": e.type,
            "children": children
        }


def _get_citations(cur: cx_Oracle.Cursor) -> Dict:
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


def _get_deleted_interpro_entries(cur: cx_Oracle.Cursor) -> List[Entry]:
    cur.execute(
        """
        SELECT E.ENTRY_AC, LOWER(T.ABBREV), E.NAME, E.SHORT_NAME, 
          E.TIMESTAMP, E.CHECKED
        FROM INTERPRO.ENTRY_AUDIT E
        LEFT OUTER JOIN INTERPRO.CV_ENTRY_TYPE T
          ON E.ENTRY_TYPE = T.CODE
        WHERE E.ENTRY_AC NOT IN (
            SELECT ENTRY_AC
            FROM INTERPRO.ENTRY
            WHERE CHECKED='Y'
        )
        ORDER BY E.TIMESTAMP
        """
    )

    entries = {}
    for row in cur:
        accession = row[0]
        entry_type = row[1]
        name = row[2]
        short_name = row[3]
        timestamp = row[4]
        is_public = row[5] == 'Y'

        try:
            e = entries[accession]
        except KeyError:
            e = Entry(accession, entry_type, name, short_name, "interpro")
            e.creation_date = timestamp
            e.deletion_date = timestamp
            e.is_alive = False
            e.is_public = is_public
        else:
            e.type = entry_type
            e.name = name
            e.short_name = short_name
            e.deletion_date = timestamp
            if is_public:
                """
                The entry was checked at least once.
                However, curators might have unchecked it right after, 
                but this becomes difficult to track.
                """
                e.is_public = True

    public_entries = []
    for e in entries.values():
        if e.is_public:
            # Only expose entries that were public at some point
            public_entries.append(e)

    return public_entries


def _get_interpro_entries(cur: cx_Oracle.Cursor) -> List[Entry]:
    cur.execute(
        """
        SELECT
          E.ENTRY_AC, LOWER(ET.ABBREV), E.NAME, E.SHORT_NAME,
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
        SELECT I2G.ENTRY_AC, GT.GO_ID, GT.NAME, GT.CATEGORY, GC.TERM_NAME
        FROM INTERPRO.INTERPRO2GO I2G
        INNER JOIN GO.TERMS@GOAPRO GT
          ON I2G.GO_ID = GT.GO_ID
        INNER JOIN GO.CV_CATEGORIES@GOAPRO GC
          ON GT.CATEGORY = GC.CODE
        """
    )

    for accession, go_id, name, category, term_name in cur:
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
                    "name": term_name
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
            children_of[parent_acc].add(child_acc)
        except KeyError:
            children_of[parent_acc] = {child_acc}

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

    """
    Cross-references
    Exclude the following databases:
        * C: PANDIT (outdated)
        * E: MSDsite (incorporated in PDB)
        * b: PDB (structures accessible from the "Structures" tab)
        * L: Blocks (outdated)
        * e: ENZYME (mapping ENZYME->UniProt->InterPro done later)
        * h, y: CATH, SCOP 
    """
    cur.execute(
        """
        SELECT X.ENTRY_AC, X.AC, LOWER(D.DBSHORT)
        FROM INTERPRO.ENTRY_XREF X
        INNER JOIN INTERPRO.CV_DATABASE D ON X.DBCODE = D.DBCODE
        WHERE X.DBCODE NOT IN ('C', 'E', 'b', 'L', 'e', 'h', 'y')
        """
    )

    for accession, ref_id, ref_db in cur:
        try:
            e = entries[accession]
        except KeyError:
            continue

        try:
            e.cross_references[ref_db].add(ref_id)
        except KeyError:
            e.cross_references[ref_db] = {ref_id}

    return list(entries.values())


def _get_signatures(cur: cx_Oracle.Cursor) -> List[Entry]:
    cur.execute(
        """
        SELECT
          M.METHOD_AC, M.NAME, M.DESCRIPTION, M.ABSTRACT, M.ABSTRACT_LONG, 
          M.METHOD_DATE, LOWER(ET.ABBREV), LOWER(DB.DBSHORT), E2M.ENTRY_AC
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


def export_entries(url: str, output: str):
    con = cx_Oracle.connect(url)
    cur = con.cursor()

    entries = {}
    logger.info("loading active InterPro entries")
    for entry in _get_interpro_entries(cur):
        entries[entry.accession] = entry

    logger.info("loading deleted InterPro entries")
    for entry in _get_deleted_interpro_entries(cur):
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

    cur.close()
    con.close()

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

    datadump(output, entries)


def export_taxonomy(url: str, output: str):
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    cur.execute(
        """
        SELECT TO_CHAR(TAX_ID), TO_CHAR(PARENT_ID), SCIENTIFIC_NAME, 
          FULL_NAME, RANK
        FROM INTERPRO.ETAXI
        """
    )

    taxonomy = {}
    for row in cur:
        taxon_id = row[0]

        taxonomy[taxon_id] = {
            "parent": row[1],
            "sci_name": row[2],
            "full_name": row[3],
            "rank": row[4],
            "children": set(),
            "lineage": [taxon_id]
        }

    cur.close()
    con.close()

    for taxon_id, taxon in taxonomy.items():
        node_id = taxon_id
        parent_id = taxon["parent"]

        # Traverse lineage from child to parent
        while parent_id is not None:
            taxon["lineage"].append(parent_id)
            taxonomy[parent_id]["children"].add(node_id)

            # We move to the parent
            node_id = parent_id
            parent_id = taxonomy[parent_id]["parent"]

    for taxon_id, info in taxonomy.items():
        info["children"] = list(info["children"])
        info["lineage"] = list(map(str, reversed(info["lineage"])))

    datadump(output, taxonomy)


class EntrySet(object):
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


def get_sets(url: str) -> Dict[str, EntrySet]:
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

    sets = {}
    for row in cur:
        set_acc = row[0]
        name = row[1]
        descr = row[2]
        database = row[3]
        member_acc = row[4]
        seq_length = row[5]
        score = row[6]

        try:
            s = sets[set_acc]
        except KeyError:
            s = sets[set_acc] = EntrySet(set_acc, name, descr, database)

        s.members.append((member_acc, score, seq_length))

    cur.close()
    con.close()

    return sets


def iter_set_alignments(url: str):
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


def get_isoforms(url: str):
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    cur.execute(
        """
        SELECT EM.METHOD_AC, EM.ENTRY_AC 
        FROM INTERPRO.ENTRY2METHOD EM
        INNER JOIN INTERPRO.ENTRY E ON EM.ENTRY_AC = E.ENTRY_AC
        WHERE E.CHECKED = 'Y'
        """
    )
    integrated = dict(cur.fetchall())

    cur.execute(
        """
        SELECT V.PROTEIN_AC, V.VARIANT, V.LENGTH, P.SEQ_SHORT, P.SEQ_LONG
        FROM INTERPRO.VARSPLIC_MASTER V
        INNER JOIN UNIPARC.PROTEIN P 
          ON V.CRC64 = P.CRC64
        """
    )

    isoforms = {}
    for row in cur:
        variant_acc = row[0] + '-' + str(row[1])
        isoforms[variant_acc] = {
            "protein_acc": row[0],
            "length": row[2],
            "sequence": row[4].read() if row[4] is not None else row[3],
            "matches": []
        }

    cur.execute(
        """
        SELECT PROTEIN_AC, METHOD_AC, MODEL_AC, POS_FROM, POS_TO, FRAGMENTS
        FROM INTERPRO.VARSPLIC_MATCH M
        """
    )

    for row in cur:
        # PROTEIN_AC is actually PROTEIN-VARIANT (e.g. Q13733-1)
        variant_acc = row[0]
        signature_acc = row[1]
        model_acc = row[2] if signature_acc != row[2] else None
        pos_start = row[3]
        pos_end = row[4]
        fragments_str = row[5]

        try:
            isoform = isoforms[variant_acc]
        except KeyError:
            continue

        if fragments_str is None:
            fragments = [{
                "start": pos_start,
                "end": pos_end,
                "dc-status": DC_STATUSES['S']  # Continuous
            }]
        else:
            fragments = []
            for frag in fragments_str.split(','):
                # Format: START-END-STATUS
                s, e, t = frag.split('-')
                fragments.append({
                    "start": int(s),
                    "end": int(e),
                    "dc-status": DC_STATUSES[t]
                })
            fragments.sort(key=repr_fragment)

        isoform["matches"].append({
            "accession": signature_acc,
            "condense": integrated.get(signature_acc),
            "fragments": fragments,
            "model": model_acc
        })

    cur.close()
    con.close()

    for accession, variant in isoforms.items():
        yield (
            accession,
            variant["protein_acc"],
            variant["length"],
            variant["sequence"],
            _post_matches(variant["matches"])
        )
