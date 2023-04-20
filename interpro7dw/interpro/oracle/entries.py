import bisect
import json
import os
import pickle
import re
from dataclasses import dataclass, field
from datetime import datetime

import cx_Oracle

from interpro7dw import intact, uniprot
from interpro7dw.utils import logger, oracle


@dataclass
class Entry:
    accession: str
    database: str
    name: str
    short_name: str
    type: str
    creation_date: datetime
    descriptions: list = field(default_factory=list, init=False)
    go_terms: list = field(default_factory=list, init=False)
    literature: dict = field(default_factory=dict, init=False)
    public: str = True

    # InterPro only
    cross_references: dict = field(default_factory=dict, init=False)
    deletion_date: datetime = field(default=None, init=False)
    old_names: list = field(default_factory=list, init=False)
    old_integrations: dict = field(default_factory=dict, init=False)
    parent: str = field(default=None, init=False)
    ppi: list = field(default_factory=list, init=False)

    # Member database only
    evidence: str = field(default=None, init=False)
    integrated_in: str = field(default=None, init=False)


DoE = dict[str, Entry]


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
    for entry_acc, databases in entry2pathways.items():
        for database, pathways in databases:
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


def _get_active_interpro_entries(cur: cx_Oracle.Cursor) -> DoE:
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
            entry = entries[accession] = Entry(accession, "INTERPRO", name,
                                               short_name, _type, date)

        if text:
            entry.descriptions.append((descr_id, text))

    # Sorts descriptions
    for accession, entry in entries.items():
        if not entry.descriptions:
            raise ValueError(f"{accession}: no descriptions")

        descriptions = []
        for descr_id, text in sorted(entry.descriptions):
            descriptions.append(text)

        entry.descriptions = descriptions

    return entries


def _add_go_terms(cur: cx_Oracle.Cursor, goa_url: str, entries: DoE):
    entry2go = {}
    cur.execute(
        """
        SELECT ENTRY_AC, GO_ID FROM INTERPRO.INTERPRO2GO
        UNION ALL
        SELECT DISTINCT METHOD_AC, GO_ID FROM INTERPRO.PANTHER2GO
        """
    )
    for accession, go_id in cur:
        if accession not in entries:
            continue

        try:
            entry2go[accession].append(go_id)
        except KeyError:
            entry2go[accession] = [go_id]

    # Gets GO terms from GOA.
    go_terms = uniprot.goa.get_terms(goa_url)

    while entry2go:
        accession, term_ids = entry2go.popitem()
        terms = []

        for go_id in term_ids:
            try:
                name, aspect, aspect_full, order = go_terms[go_id]
            except KeyError:
                logger.error(f"{accession}: term {go_id} not found")
                continue

            terms.append((order, go_id, name, aspect, aspect_full))

        entry = entries[accession]

        # Sort terms
        for _, go_id, name, aspect, aspect_full in sorted(terms):
            entry.go_terms.append({
                "identifier": go_id,
                "name": name,
                "category": {
                    "code": aspect,  # e.g. F
                    "name": aspect_full  # e.g. molecular_function
                }
            })


def _add_hierarchies(cur: cx_Oracle.Cursor, entries: dict[str, Entry]):
    cur.execute(
        """
        SELECT ENTRY_AC, PARENT_AC
        FROM INTERPRO.ENTRY2ENTRY
        WHERE RELATION = 'TY'
        """
    )
    for entry_acc, parent_acc in cur:
        # TODO: raise again
        if entry_acc not in entries:
            continue
            # raise KeyError(f"{entry_acc}: unchecked entry in hierarchy")
        elif parent_acc not in entries:
            continue
            # raise KeyError(f"{parent_acc}: unchecked entry in hierarchy")

        entries[entry_acc].parent = parent_acc


def _add_xrefs(cur: cx_Oracle.Cursor, entries: DoE):
    cur.execute(
        """
        SELECT X.ENTRY_AC, D.DBSHORT, X.AC
        FROM INTERPRO.ENTRY_XREF X
        INNER JOIN INTERPRO.CV_DATABASE D ON X.DBCODE = D.DBCODE
        """
    )

    for entry_acc, xref_db, xref_id in cur:
        if entry_acc not in entries:
            continue

        entry = entries[entry_acc]
        if xref_db in entry.cross_references:
            entry.cross_references[xref_db].append(xref_id)
        else:
            entry.cross_references[xref_db] = [xref_id]


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


def _get_retired_interpro_entries(cur: cx_Oracle.Cursor) -> DoE:
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
            entry = entries[acc]
        else:
            entry = entries[acc] = Entry(acc, "INTERPRO", name, short_name,
                                         _type, timestamp)

        entry.public = False
        entry.deletion_date = timestamp
        entry.name = name
        entry.short_name = short_name
        entry.type = _type

    results = {}
    for acc, entry in entries.items():
        if any(entries_per_release[acc].values()):
            # Entry public for at least once release
            results[acc] = entry

    return results


def _get_signatures(cur: cx_Oracle.Cursor) -> DoE:
    cur.execute(
        """
        SELECT
          M.METHOD_AC, M.NAME, M.DESCRIPTION, M.ABSTRACT, M.ABSTRACT_LONG,
          M.METHOD_DATE, ET.ABBREV, DB.DBSHORT, E2M.ENTRY_AC, EVI.ABBREV
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
        UNION ALL
        SELECT  -- FunFams
          F.METHOD_AC, F.NAME, F.DESCRIPTION, NULL, NULL,
          F.METHOD_DATE, 'Region', DB.DBSHORT, NULL, EVI.ABBREV
        FROM INTERPRO.FEATURE_METHOD F
        INNER JOIN INTERPRO.CV_DATABASE DB ON F.DBCODE = DB.DBCODE
        LEFT OUTER JOIN INTERPRO.IPRSCAN2DBCODE I2D
          ON F.DBCODE = I2D.DBCODE
        LEFT OUTER JOIN INTERPRO.CV_EVIDENCE EVI
          ON I2D.EVIDENCE = EVI.CODE
        WHERE F.DBCODE = 'f'
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

        signature = Entry(acc, database, name, short_name, _type, date)
        signature.evidence = evidence
        signature.integrated_in = interpro_acc
        if descr_text:
            signature.descriptions.append(descr_text)

        signatures[acc] = signature

    # Update PANTHER subfamilies and CATH FunFams
    panther_subfamily = re.compile(r"(PTHR\d+):SF\d+")
    cath_funfams = re.compile(r"(G3DSA:\d+\.\d+\.\d+\.\d+):FF:\d+")
    for acc, signature in signatures.items():
        m = panther_subfamily.fullmatch(acc)
        if m:
            family_acc = m.group(1)

            if family_acc in signatures:
                signatures[acc].integrated_in = family_acc
                signatures[acc].parent = family_acc
                signatures[acc].public = False
            else:
                raise KeyError(f"PANTHER family {family_acc} not found "
                               f"for subfamily {acc}")

            continue

        m = cath_funfams.fullmatch(acc)
        if m:
            supfam_acc = m.group(1)

            if supfam_acc in signatures:
                signatures[acc].integrated_in = supfam_acc
                signatures[acc].parent = supfam_acc
                signatures[acc].public = False
            else:
                raise KeyError(f"CATH-Gene3D superfamily {supfam_acc} "
                               f"not found for family {acc}")

    return signatures


def _add_citations(cur: cx_Oracle.Cursor, entries: DoE, signatures: DoE):
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
            entries[acc].literature = entry2pub[acc]

    cur.execute(
        """
        SELECT METHOD_AC, PUB_ID
        FROM INTERPRO.METHOD2PUB        
        """
    )

    entry2pub = {}
    for accession, pub_id in cur:
        try:
            entry2pub[accession][pub_id] = citations[pub_id]
        except KeyError:
            entry2pub[accession] = {pub_id: citations[pub_id]}

    for acc in entry2pub:
        if acc in signatures:
            signatures[acc].literature = entry2pub[acc]


def export_entries(interpro_uri: str, goa_uri: str, intact_uri: str,
                   output: str):
    """Export InterPro entries and member database signatures.

    :param interpro_uri: InterPro Oracle connection string.
    :param goa_uri: GOA Oracle connection string.
    :param intact_uri:  IntAct Oracle connection string.
    :param output: Output file.
    """
    logger.info("loading from Oracle databases")
    con = cx_Oracle.connect(interpro_uri)
    cur = con.cursor()
    # fetch CLOB object as strings
    cur.outputtypehandler = oracle.lob_as_str

    # Starts with InterPro entries
    entries = _get_active_interpro_entries(cur)

    # Adds entry hierarchies
    _add_hierarchies(cur, entries)

    # Adds cross-references
    _add_xrefs(cur, entries)

    # Adds protein-protein interactions from IntAct
    for acc, ppi in intact.get_interactions(intact_uri).items():
        if acc in entries:
            entries[acc].ppi = ppi

    # Adds retired entries (that were at least public in one release)
    for acc, entry in _get_retired_interpro_entries(cur).items():
        entries[acc] = entry

    # Add past names
    for acc, old_names in _get_past_names(cur).items():
        if acc in entries:
            entries[acc].old_names = old_names

    # Add past integrations
    for acc, mem_dbs in _get_past_integrations(cur).items():
        if acc in entries:
            entries[acc].old_integrations = mem_dbs

    signatures = _get_signatures(cur)

    # Adds literature references
    _add_citations(cur, entries, signatures)

    while signatures:
        k, v = signatures.popitem()
        entries[k] = v

    # Adds GO terms (InterPro + PANTHER)
    _add_go_terms(cur, goa_uri, entries)

    cur.close()
    con.close()

    with open(output, "wb") as fh:
        pickle.dump(entries, fh)

    logger.info("done")


def update_ipr_data(ipr_data, ipr, pid):
    if ipr in ipr_data:
        pdata1 = ipr_data[ipr]
    else:
        pdata1 = []
    pdata1.append(pid)
    ipr_data[ipr] = pdata1

    return ipr_data


def _export_pathways(cur: cx_Oracle.Cursor, output_path: str):
    cur.execute("""
        SELECT ENTRY_AC, DBCODE, AC, NAME
        FROM INTERPRO.ENTRY2PATHWAY
        """)
    pathway_data = cur.fetchall()

    pdata = {}
    ipr_data = {}
    for el in pathway_data:
        ipr = el[0]
        pid = el[2]
        ipr_data = update_ipr_data(ipr_data, ipr, pid)
        pdata[pid] = el

    with open(os.path.join(output_path, "pathways.json"), "wt") as fh:
        json.dump(pdata, fh)

    with open(os.path.join(output_path, "pathways.ipr.json"), "wt") as fh:
        json.dump(ipr_data, fh)

    logger.info("Done.")


def _export_go_terms(cur: cx_Oracle.Cursor, output_path: str):
    cur.execute("""
        SELECT i2g.entry_ac, g.go_id, g.name, g.category
        FROM INTERPRO.INTERPRO2GO i2g
        INNER JOIN INTERPRO.ENTRY e ON e.entry_ac = i2g.entry_ac
        JOIN go.terms@GOAPRO g ON i2g.go_id = g.go_id
        WHERE e.checked='Y'
        """)
    goterms_data = cur.fetchall()

    godata = {}
    ipr_data = {}

    for el in goterms_data:
        ipr = el[0]
        goid = el[1]
        ipr_data = update_ipr_data(ipr_data, ipr, goid)
        godata[goid] = el

    with open(os.path.join(output_path, "goterms.json"), "wt") as fh:
        json.dump(godata, fh)

    with open(os.path.join(output_path, "goterms.ipr.json"), "wt") as fh:
        json.dump(ipr_data, fh)

    logger.info("Done.")


def export_for_interproscan(uri: str, outdir: str):
    con = cx_Oracle.connect(uri)
    cur = con.cursor()

    _export_pathways(cur, outdir)
    _export_go_terms(cur, outdir)

    cur.close()
    con.close()
