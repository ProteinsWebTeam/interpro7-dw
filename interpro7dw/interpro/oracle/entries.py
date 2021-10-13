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
]
FEATURE_DATABASES = [
    'g',    # MobiDB Lite
    'j',    # Phobius
    'n',    # Signal Euk
    'q',    # TMHMM
    's',    # SignalP Gram positive
    'v',    # SignalP Gram negative
    'x',    # COILS
]
SEQUENCE_DATABASES = [
    'S',    # Swiss-Prot
    'T',    # TrEMBL
    'u',    # UniProtKB
]


def dump_databases(url: str, version: str, date: str, file: str,
                   update: bool = False):
    """Exports information on databases/data sources used in InterPro.

    :param url: The Oracle connection string.
    :param version: The version of the upcoming InterPro release.
    :param date: The date of the upcoming InterPro release (YYYY-MM-DD).
    :param file: The output file.
    :param update: If True, update the production table.
    """
    con = cx_Oracle.connect(url)
    cur = con.cursor()

    cur.execute("SELECT COUNT(*) FROM INTERPRO.ENTRY WHERE CHECKED = 'Y'")
    num_interpro_entries, = cur.fetchone()

    cur.execute("SELECT VERSION FROM INTERPRO.DB_VERSION WHERE DBCODE = 'I'")
    prod_version, = cur.fetchone()

    if prod_version == version:
        # DB_VERSION is already up-to-date
        use_db_version = True
    elif update:
        # DB_VERSION is outdated, but will be up-to-date
        use_db_version = True
        cur.execute(
            """
            UPDATE INTERPRO.DB_VERSION
            SET VERSION = :1,
                FILE_DATE = :2,
                ENTRY_COUNT = :3
            WHERE DBCODE = 'I'
            """, (version, datetime.strptime(date, "%Y-%m-%d"),
                  num_interpro_entries)
        )
        con.commit()
    else:
        # DB_VERSION is outdated and will stay outdated
        # This run is a test done on the production database (SRSLY?!!111)
        use_db_version = False

    """
    Using RN=2 to join with the second most recent action in DB_VERSION_AUDIT
    (the most recent is the same record as in DB_VERSION)
    """
    cur.execute(
        """
        SELECT
          DB.DBCODE, LOWER(DB.DBSHORT), DB.DBSHORT, DB.DBNAME, 
          DB.DESCRIPTION, V.VERSION, V.FILE_DATE, V.ENTRY_COUNT, VA.VERSION, 
          VA.FILE_DATE
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
    for rec in cur:
        code = rec[0]
        identifier = rec[1]
        short_name = rec[2]
        name = rec[3]
        description = rec[4]
        release_version = rec[5]
        release_date = rec[6]
        num_entries = rec[7]
        prev_release_version = rec[8]
        prev_release_date = rec[9]

        if code in ENTRY_DATABASES:
            db_type = "entry"
        elif code in FEATURE_DATABASES:
            db_type = "feature"
        elif code in SEQUENCE_DATABASES:
            if code == 'S':
                identifier = "reviewed"
            elif code == 'T':
                identifier = "unreviewed"

            db_type = "protein"
        else:
            db_type = "other"

        if code == 'I' and not use_db_version:
            # DB_VERSION is outdated:
            # it contains info for the live release (soon to be 'previous')
            num_entries = num_interpro_entries
            prev_release_version = release_version
            prev_release_date = release_date
            release_version = version
            release_date = datetime.strptime(date, "%Y-%m-%d")

        databases.append((
            identifier,
            name,
            short_name,
            description,
            db_type,
            num_entries,
            release_version,
            release_date,
            prev_release_version,
            prev_release_date
        ))

    cur.close()
    con.close()

    dumpobj(databases, file)


def dump_domain_organisation(url: str, proteins_src: str, matches_src: str,
                             domorgs_dst: str, **kwargs):
    """Calculates and exports the domain architectures/organisations of
    UniProt entries based on the Pfam matches.

    :param url: The Oracle connection string.
    :param proteins_src: The file containing protein information.
    :param matches_src: The file containing protein matches.
    :param domorgs_dst: The output domain organisation file.
    """
    tempdir = kwargs.get("tempdir")
    workers = kwargs.get("workers", 1)

    # Loads Pfam signatures, and the InterPro entries they are integrated in
    logger.info("loading Pfam signatures")
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    cur = cur.execute(
        """
        SELECT M.METHOD_AC, EM.ENTRY_AC
        FROM INTERPRO.METHOD M
        LEFT OUTER JOIN (
            SELECT EM.METHOD_AC, E.ENTRY_AC
            FROM INTERPRO.ENTRY2METHOD EM
            INNER JOIN INTERPRO.ENTRY E ON EM.ENTRY_AC = E.ENTRY_AC
            WHERE E.CHECKED = 'Y'
        ) EM ON M.METHOD_AC = EM.METHOD_AC
        WHERE M.DBCODE = 'H'  -- Pfam
        """
    )

    pfam_signatures = {}
    for pfam_acc, interpro_acc in cur:
        pfam_signatures[pfam_acc] = interpro_acc

    cur.close()
    con.close()

    logger.info("iterating proteins")
    all_domains = {}
    with SimpleStore(tempdir=tempdir) as tmp:
        with Store(proteins_src, "r") as st1, Store(matches_src, "r") as st2:
            keys = st1.file_keys

            for i, (protein_acc, matches) in enumerate(st2.items()):
                locations = []
                for entry_acc in matches:
                    try:
                        interpro_acc = pfam_signatures[entry_acc]
                    except KeyError:
                        # Not a Pfam match
                        continue

                    for loc in matches[entry_acc]:
                        locations.append({
                            "pfam": entry_acc,
                            "interpro": interpro_acc,
                            # We do not consider fragmented locations
                            "start": loc["fragments"][0]["start"],
                            "end": max(f["end"] for f in loc["fragments"])
                        })

                if (i + 1) % 10000000 == 0:
                    logger.info(f"{i + 1:>15,}")

                if not locations:
                    continue  # No Pfam matches: no domain organisation

                domains = []
                members = set()
                for loc in sorted(locations, key=lambda l: (l["start"],
                                                            l["end"])):
                    pfam_acc = loc["pfam"]
                    interpro_acc = loc["interpro"]

                    if interpro_acc:
                        domains.append(f"{pfam_acc}:{interpro_acc}")
                        members.add(interpro_acc)
                    else:
                        domains.append(pfam_acc)

                    members.add(pfam_acc)

                dom_str = '-'.join(domains)
                dom_id = hashlib.sha1(dom_str.encode("utf-8")).hexdigest()
                tmp.add((protein_acc, dom_str, dom_id, members))

                # string (YYYY-MM-DD) which is enough to compare dates
                date = st1[protein_acc]["date"]

                # Selects the oldest protein to represent
                # this domain organisation.
                try:
                    other_date, _ = all_domains[dom_id]
                except KeyError:
                    all_domains[dom_id] = (date, protein_acc)
                else:
                    if date < other_date:
                        all_domains[dom_id] = (date, protein_acc)

            logger.info(f"{i + 1:>15,}")

        size = tmp.size

        logger.info("exporting domain organisations")
        with Store(domorgs_dst, mode="w", keys=keys, tempdir=tempdir) as st:
            for i, (protein_acc, dom_str, dom_id, members) in enumerate(tmp):
                _, repr_acc = all_domains[dom_id]

                st.add(protein_acc, (dom_str, dom_id, members, repr_acc))

                if (i + 1) % 10000000 == 0:
                    logger.info(f"{i + 1:>15,}")

            logger.info(f"{i + 1:>15,}")

            size += st.size
            st.merge(workers, apply=st.get_first)

        logger.info(f"temporary files: {size / 1024 / 1024:.0f} MB")

    logger.info("done")


def dump_similar_entries(url: str, matches_file: str, overlapping_file: str,
                         min_similarity: float = 0.75):
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    cur.execute(
        """
        SELECT E.ENTRY_AC, LOWER(ET.ABBREV)
        FROM INTERPRO.ENTRY E
        INNER JOIN INTERPRO.CV_ENTRY_TYPE ET ON E.ENTRY_TYPE = ET.CODE
        WHERE E.CHECKED = 'Y'
        """
    )
    entries = dict(cur.fetchall())
    cur.close()
    con.close()

    num_proteins = {}  # number of proteins matched by entry
    num_overlaps = {}  # number of proteins where two entries overlap >= 50%

    logger.info("iterating proteins")
    with Store(matches_file, "r") as store:
        for i, (protein_acc, matches) in enumerate(store.items()):
            entries = {}
            for entry_acc, locations in matches.items():
                if entry_acc in entries:
                    entries[entry_acc] = []

                    for loc in locations:
                        # InterPro locations have one fragment only
                        entries[entry_acc].append((
                            loc["fragments"][0]["start"],
                            loc["fragments"][0]["end"],
                        ))

            # Evaluate how entries overlap
            for entry_acc, locations in entries.items():
                try:
                    num_proteins[entry_acc] += 1
                except KeyError:
                    num_proteins[entry_acc] = 1

                for other_acc, other_locations in entries.items():
                    if other_acc >= entry_acc:
                        continue

                    try:
                        entry_overlaps = num_overlaps[entry_acc]
                    except KeyError:
                        entry_overlaps = num_overlaps[entry_acc] = {}

                    try:
                        overlaps = entry_overlaps[other_acc]
                    except KeyError:
                        overlaps = entry_overlaps[other_acc] = [0, 0]

                    flag = 0
                    for start1, end1 in locations:
                        length1 = end1 - start1 + 1

                        for start2, end2 in other_locations:
                            length2 = end2 - start2 + 1
                            overlap = min(end1, end2) - max(start1, start2) + 1

                            if not flag & 1 and overlap >= length1 * 0.5:
                                # 1st time fragments overlap >= 50% of entry 1
                                flag |= 1
                                overlaps[0] += 1

                            if not flag & 2 and overlap >= length2 * 0.5:
                                # 1st time fragments overlap >= 50% of entry 2
                                flag |= 2
                                overlaps[1] += 1

                        if flag == 3:
                            # Both cases already happened
                            break

            if (i + 1) % 10000000 == 0:
                logger.info(f"{i + 1:>15,}")

        logger.info(f"{i + 1:>15,}")

    supfam = "homologous_superfamily"
    types = (supfam, "domain", "family", "repeat")
    overlapping_entries = []
    for entry_acc, overlaps in num_overlaps.items():
        entry_cnt = num_proteins[entry_acc]

        for other_acc, (cnt1, cnt2) in overlaps.items():
            other_cnt = num_proteins[other_acc]

            # Independent coefficients
            coef1 = cnt1 / (entry_cnt + other_cnt - cnt1)
            coef2 = cnt2 / (entry_cnt + other_cnt - cnt2)

            # Final coefficient (average of independent coefficients)
            coef = (coef1 + coef2) * 0.5

            # Containment indices
            cont1 = cnt1 / entry_cnt
            cont2 = cnt2 / other_cnt

            if all(e < min_similarity for e in (coef, cont1, cont2)):
                continue

            # Entries are deemed similar
            type1 = entries[entry_acc]
            type2 = entries[other_acc]
            if ((type1 == supfam and type2 in types)
                    or (type2 == supfam and type1 in types)):
                overlapping_entries.append((entry_acc, other_acc))

    dumpobj(overlapping_entries, overlapping_file)
    logger.info("done")


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
            "interpro": {
                "id": row[4],
                "name": row[5],
                "type": row[6],
                "parent_id": row[7],
            } if row[4] else None
        }
    return signatures


class Entry:
    def __init__(self, accession: str, short_name: str, name: str,
                 entry_type: str, database: str):
        self.accession = accession
        self.shot_name = short_name
        self.name = name
        self.type = entry_type
        self.database = database
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
        self.integrates = {}        # InterPro only
        self.literature = {}        # all entries
        self.overlaps_with = []     # InterPro only
        self.pathways = {}          # InterPro only
        self.ppi = []               # prot-prot interactions (InterPro only)
        self.xrefs = {}             # InterPro only


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
                                              _type, "interpro")
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
            e = entries[acc] = Entry(acc, short_name, name, _type, "interpro")
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
            WHERE CHECKED = 'Y'
          )
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

        e = entries[acc] = Entry(acc, short_name, name, _type, database)
        e.creation_date = date
        e.descriptions.append(descr_text)
        e.integrated_in = interpro_acc

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

            entry.xrefs["ec"] = sorted(xrefs["enzymes"])

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
                "taxa": len(xrefs["taxa"]),
            }

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
                    cur.executemany(
                        """
                        INSERT INTO INTERPRO.ENTRY2PATHWAY 
                        VALUES (:1, :2, :3, :4)
                        """, (entry.accession, dbcode, pathway_id, name)
                    )

        con.commit()
        cur.close()
        con.close()

    logger.info("done")
