import bisect

import cx_Oracle

from interpro7dw import intact, uniprot
from interpro7dw.utils import logger, oracle


def update_pathways(uri: str, entry2pathways: dict[str, list[tuple]]):
    con = cx_Oracle.connect(uri)
    cur = con.cursor()
    cur.execute("TRUNCATE TABLE INTERPRO.ENTRY2PATHWAY")
    cur.execute(
        """
        SELECT LOWER(DBSHORT), DBCODE
        FROM INTERPRO.CV_DATABASE
        """
    )
    id2dbcode = dict(cur.fetchall())

    params = []
    for entry_acc, databases in entry2pathways.values():
        for database, pathways in databases.items():
            dbcode = id2dbcode[database.lower()]

            for pathway_id, pathway_name in pathways:
                params.append((entry_acc, dbcode, pathway_id, pathway_name))

    for i in range(0, len(params), 1000):
        cur.executemany(
            """
            INSERT INTO INTERPRO.ENTRY2PATHWAY (ENTRY_AC, DBCODE, AC, NAME)
            VALUES (:1, :2, :3, :4)
            """,
            params[i:i+1000]
        )

    con.commit()
    cur.close()
    con.close()


def _get_active_interpro_entries(cur: cx_Oracle.Cursor) -> dict:
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

    for accession, _type, name, short_name, date, descr_id, text in cur:
        try:
            entry = entries[accession]
        except KeyError:
            entry = entries[accession] = {
                "accession": accession,
                "creation_date": date,
                "descriptions": [],
                "database": "InterPro",
                "name": name,
                "short_name": short_name,
                "type": _type,
            }

        if text:
            entry["descriptions"].append((descr_id, text))

    # Sorts descriptions
    for accession, entry in entries.items():
        if not entry["descriptions"]:
            raise ValueError(f"{accession}: no descriptions")

        descriptions = []
        for descr_id, text in sorted(entry["descriptions"]):
            descriptions.append(text)

        entry["descriptions"] = descriptions

    return entries


def _add_go_terms(cur: cx_Oracle.Cursor, goa_url: str, entries: dict):
    interpro2go = {}
    cur.execute("SELECT ENTRY_AC, GO_ID FROM INTERPRO.INTERPRO2GO")
    for accession, go_id in cur:
        if accession not in entries:
            continue

        try:
            interpro2go[accession].append(go_id)
        except KeyError:
            interpro2go[accession] = [go_id]

    # Gets GO terms from GOA.
    go_terms = uniprot.goa.get_terms(goa_url)

    while interpro2go:
        accession, term_ids = interpro2go.popitem()
        terms = []

        for go_id in term_ids:
            try:
                name, aspect, aspect_full, order = go_terms[go_id]
            except KeyError:
                logger.error(f"{accession}: term {go_id} not found")
                continue

            terms.append((order, go_id, name, aspect, aspect_full))

        entry = entries[accession]
        entry["go_terms"] = []

        # Sort terms
        for _, go_id, name, aspect, aspect_full in sorted(terms):
            entry["go_terms"].append({
                "identifier": go_id,
                "name": name,
                "category": {
                    "code": aspect,  # e.g. F
                    "name": aspect_full  # e.g. molecular_function
                }
            })


def _add_hierarchies(cur: cx_Oracle.Cursor, entries: dict):
    child2parent = {}
    parent2children = {}
    cur.execute(
        """
        SELECT ENTRY_AC, PARENT_AC
        FROM INTERPRO.ENTRY2ENTRY
        WHERE RELATION = 'TY'
        """
    )
    for entry_acc, parent_acc in cur:
        if entry_acc not in entries:
            raise KeyError(f"{entry_acc}: unchecked entry in hierarchy")
        elif parent_acc not in entries:
            raise KeyError(f"{parent_acc}: unchecked entry in hierarchy")

        entries[entry_acc]["parent"] = parent_acc


def _add_xrefs(cur: cx_Oracle.Cursor, entries: dict):
    cur.execute(
        """
        SELECT X.ENTRY_AC, X.AC, D.DBSHORT
        FROM INTERPRO.ENTRY_XREF X
        INNER JOIN INTERPRO.CV_DATABASE D ON X.DBCODE = D.DBCODE
        """
    )

    xrefs = {}
    for accession, xref_id, xref_db in cur:
        if accession not in xrefs:
            xrefs[accession] = {xref_db: [xref_id]}
        elif xref_db in xrefs[accession]:
            xrefs[accession][xref_db].append(xref_id)
        else:
            xrefs[accession][xref_db] = [xref_id]

    for accession in xrefs:
        if accession in entries:
            entries[accession]["cross_references"] = xrefs[accession]


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


def _get_past_names(cur: cx_Oracle.Cursor) -> dict[str, list[str]]:
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


def _get_past_integrations(cur: cx_Oracle.Cursor) -> dict:
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
        SELECT M.METHOD_AC, DB.DBSHORT, E.ENTRY_AC
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


def _get_retired_interpro_entries(cur: cx_Oracle.Cursor) -> dict[str, dict]:
    """Returns a list of InterPro entries that are not public anymore
    (i.e. deleted of checked=N, including in the upcoming release).

    Only entries that were public at least once are returned.
    We use production freeze times to evaluate if an entry was still existing
    and checked=Y for a release.

    :param cur: Oracle cursor.
    :return: A dict of entries.
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

        if acc in entries:
            # Update entry
            entries[acc].update({
                "deletion_date": timestamp,
                "name": name,
                "short_name": short_name,
                "type": _type,
            })
        else:
            entries[acc] = {
                "accession": acc,
                "creation_date": timestamp,
                "database": "InterPro",
                "deletion_date": timestamp,
                "name": name,
                "short_name": short_name,
                "type": _type,
            }

    results = {}
    for acc, entry in entries.items():
        if any(entries_per_release[acc].values()):
            # Entry public for at least once release
            results[acc] = entry

    return results


def _get_signatures(cur: cx_Oracle.Cursor) -> dict[str, dict]:
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

    signatures = {}
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

        signatures[acc] = {
            "accession": acc,
            "creation_date": date,
            "descriptions": [descr_text] if descr_text else [],
            "database": database,
            "evidence": evidence,
            "integrated_in": interpro_acc,
            "name": name,
            "short_name": short_name,
            "type": _type,
        }

    return signatures


def _add_citations(cur: cx_Oracle.Cursor, entries: dict[str, dict],
                   signatures: dict[str, dict]):
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
        """
    )

    entry2pub = {}
    for accession, pub_id in cur:
        try:
            entry2pub[accession][pub_id] = citations[pub_id]
        except KeyError:
            entry2pub[accession] = {pub_id: citations[pub_id]}

    for acc in entry2pub:
        if acc in entries:
            entries[acc]["literature"] = entry2pub[acc]

    cur.execute(
        """
        SELECT METHOD_AC, PUB_ID
        FROM INTERPRO.METHOD2PUB        
        """
    )

    entry2pub.clear()
    for accession, pub_id in cur:
        try:
            entry2pub[accession][pub_id] = citations[pub_id]
        except KeyError:
            entry2pub[accession] = {pub_id: citations[pub_id]}

    for acc in entry2pub:
        if acc in signatures:
            signatures[acc]["literature"] = entry2pub[acc]


def export_entries(interpro_uri: str, goa_uri: str, intact_uri: str):
    logger.info("loading from Oracle databases")
    con = cx_Oracle.connect(interpro_uri)
    cur = con.cursor()
    # fetch CLOB object as strings
    cur.outputtypehandler = oracle.lob_as_str

    # Starts with InterPro entries
    entries = _get_active_interpro_entries(cur)

    # Adds entry hierarchies
    _add_hierarchies(cur, entries)

    # Adds GO terms
    _add_go_terms(cur, goa_uri, entries)

    # Adds cross-references
    _add_xrefs(cur, entries)

    # Adds protein-protein interactions from IntAct
    for acc, ppi in intact.get_interactions(interpro_uri).items():
        if acc in entries:
            entries[acc]["ppi"] = ppi

    # Add past names
    for acc, old_names in _get_past_names(cur).items():
        if acc in entries:
            entries[acc]["old_names"] = old_names

    # Add past integrations
    for acc, mem_db in _get_past_integrations(cur).items():
        if acc in entries:
            entries[acc]["old_integrations"] = mem_db

    # Adds retired entries (that were at least public in one release)
    for acc, entry in _get_retired_interpro_entries(cur).items():
        pass

    signatures = _get_signatures(cur)

    # Adds literature references
    _add_citations(cur, entries, signatures)

    cur.close()
    con.close()

