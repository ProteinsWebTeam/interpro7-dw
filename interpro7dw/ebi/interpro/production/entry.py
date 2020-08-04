# -*- coding: utf-8 -*-

import bisect
from typing import List, Optional

import cx_Oracle

from interpro7dw import kegg, logger, metacyc
from interpro7dw.ebi import intact, uniprot
from interpro7dw.utils import Store, dumpobj, loadobj


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
                "type": e.type,
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
        SELECT E.ENTRY_AC, LOWER(T.ABBREV), E.NAME, E.SHORT_NAME, 
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
        ORDER BY GC.SORT_ORDER, I2G.GO_ID
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
            e.cross_references[ref_db].append(ref_id)
        except KeyError:
            e.cross_references[ref_db] = [ref_id]

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


def export_entries(url: str, p_metacyc: str, p_clans: str,
                   p_uniprot2matches: str, p_entries: str,
                   p_uniprot2entries: str, processes: int=1,
                   dir: Optional[str]=None):
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

    logger.info("loading Reactome pathways")
    u2reactome = uniprot.get_swissprot2reactome(cur)

    logger.info("loading ENZYME")
    u2enzyme = uniprot.get_swissprot2enzyme(cur)
    cur.close()
    con.close()

    logger.info("loading KEGG pathways")
    ec2kegg = kegg.get_ec2pathways()

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

    logger.info("exporting UniProt-entries mapping")
    u2matches = Store(p_uniprot2matches)
    u2entries = Store(p_uniprot2entries, u2matches.get_keys(), dir)
    interpro2enzyme = {}
    interpro2reactome = {}
    i = 0
    for uniprot_acc, matches in u2matches.items():
        protein_entries = []
        for entry_acc in matches:
            entry = entries[entry_acc]

            if entry.database == "interpro":
                terms = entry.go_terms

                for ecno in u2enzyme.get(uniprot_acc, []):
                    try:
                        interpro2enzyme[entry_acc].add(ecno)
                    except KeyError:
                        interpro2enzyme[entry_acc] = {ecno}

                for pathway in u2reactome.get(uniprot_acc, []):
                    try:
                        interpro2reactome[entry_acc].add(pathway)
                    except KeyError:
                        interpro2reactome[entry_acc] = {pathway}
            else:
                terms = []

            protein_entries.append((
                entry.accession,
                entry.database,
                entry.clan["accession"] if entry.clan else None,
                terms
            ))

        u2entries[uniprot_acc] = protein_entries

        i += 1
        if not i % 1000000:
            u2entries.sync()

            if not i % 10000000:
                logger.info(f"{i:>12,}")

    logger.info(f"{i:>12,}")
    u2matches.close()
    size = u2entries.merge(processes=processes)
    u2entries.close()
    logger.info(f"temporary files: {size/1024/1024:.0f} MB")

    logger.info("updating ENZYME cross-references and pathways")
    for entry_acc, ecnos in interpro2enzyme.items():
        entry = entries[entry_acc]
        entry.cross_references["ec"] = list(ecnos)

        for ecno in ecnos:
            for pathway in sorted(ec2kegg.get(ecno, [])):
                try:
                    pathways = entry.pathways["kegg"]
                except KeyError:
                    pathways = entry.pathways["kegg"] = []

                pathways.append(dict(zip(("id", "name"), pathway)))

            for pathway in sorted(ec2metacyc.get(ecno, [])):
                try:
                    pathways = entry.pathways["metacyc"]
                except KeyError:
                    pathways = entry.pathways["metacyc"] = []

                pathways.append(dict(zip(("id", "name"), pathway)))

    for entry_acc in interpro2reactome:
        entry = entries[entry_acc]
        try:
            pathways = entry.pathways["reactome"]
        except KeyError:
            pathways = entry.pathways["reactome"] = []

        for pathway in sorted(interpro2reactome[entry_acc]):
            pathways.append(dict(zip(("id", "name"), pathway)))

    dumpobj(p_entries, entries)
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
