import json
import sys
from typing import Generator, Optional, Tuple

from .. import dbms, goa, io
from . import DC_STATUSES, repr_frag


def get_databases(uri: str) -> list:
    # todo: do not hardcode this value!
    member_dbs = {
        'B', 'F', 'H', 'I', 'J', 'M', 'N', 'P', 'Q', 'R', 'U', 'V',
        'X', 'Y', 'g'
    }

    con, cur = dbms.connect(uri)

    """
    Using RN=2 to join with the second most recent action
    (the most recent is the current record)
    """
    cur.execute(
        """
        SELECT
          LOWER(DB.DBSHORT), DB.DBCODE, DB.DBNAME, DB.DESCRIPTION,
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
        name = row[0]
        code = row[1]

        if code in member_dbs:
            db_type = "entry"
        elif code in ('S', 'T', 'u'):
            if code == 'S':
                name = "reviewed"
            elif code == 'T':
                name = "unreviewed"

            db_type = "protein"
        else:
            db_type = "other"

        databases.append({
            "code": code,
            "name": name,
            "name_long": row[2],
            "description": row[3],
            "type": db_type,
            "version": {
                "code": row[4],
                "date": row[5]
            },
            "previous_version": {
                "code": row[6],
                "date": row[7]
            }
        })

    cur.close()
    con.close()

    return databases


class EntryHierarchyTree(object):
    def __init__(self, relationships):
        parent_of = {}
        children_of = {}

        for child_ac, parent_ac in relationships:
            parent_of[child_ac] = parent_ac

            if parent_ac not in children_of:
                children_of[parent_ac] = []

            children_of[parent_ac].append(child_ac)

        self._parent_of = parent_of
        self._children_of = children_of

    def children_of(self, entry_ac):
        return self._children_of.get(entry_ac, [])

    def parent_of(self, entry_ac):
        return self._parent_of.get(entry_ac)

    def get_root(self, entry_ac):
        parent_ac = self.parent_of(entry_ac)

        while parent_ac is not None:
            entry_ac = parent_ac
            parent_ac = self.parent_of(entry_ac)

        return entry_ac


def format_node(hierarchy, entries, accession):
    entry = entries[accession]
    children = hierarchy.children_of(accession)
    return {
        'accession': accession,
        'name': entry['name'],
        'type': entry['type'],
        'children': [
            format_node(hierarchy, entries, child_ac)
            for child_ac in children
        ]
    }


def get_deleted_entries(uri: str) -> list:
    con, cur = dbms.connect(uri)

    cur.execute(
        """
        SELECT
            E.ENTRY_AC, LOWER(T.ABBREV), E.NAME,
            E.SHORT_NAME, E.TIMESTAMP, E.CHECKED
        FROM INTERPRO.ENTRY_AUDIT E
        LEFT OUTER JOIN INTERPRO.CV_ENTRY_TYPE T
          ON E.ENTRY_TYPE = T.CODE
        WHERE E.ENTRY_AC NOT IN (
            SELECT ENTRY_AC
            FROM INTERPRO.ENTRY
            WHERE CHECKED='Y'
        )
        ORDER BY E.ENTRY_AC, E.TIMESTAMP
        """
    )

    entries = {}
    for row in cur:
        acc = row[0]
        _type = "unknown" if row[1] is None else row[1]

        if acc in entries:
            entries[acc].update({
                "type": _type,
                "name": row[2],
                "short_name": row[3],
                "deletion_date": row[4]
            })
            if row[5] == 'Y':
                entries[acc]["was_public"] = True
        else:
            entries[acc] = {
                "accession": acc,
                "type": _type,
                "name": row[2],
                "short_name": row[3],
                "database": "interpro",
                "creation_date": row[4],
                "deletion_date": None,
                "was_public": row[5] == 'Y'
            }

    cur.close()
    con.close()

    public_entries = []
    for e in entries.values():
        if e["was_public"]:
            # the entry was once public (otherwise we don't expose it)

            if e["deletion_date"] is None:
                """
                Some entries have only one record in the audit table:
                they miss the record for their deletion so we fallback to
                the creation date for the deletion date
                """
                e["deletion_date"] = e["creation_date"]

            public_entries.append(e)

    return public_entries


def get_relationships(cur):
    cur.execute(
        """
        SELECT ENTRY_AC, PARENT_AC
        FROM INTERPRO.ENTRY2ENTRY
        WHERE RELATION = 'TY'
        """
    )

    return EntryHierarchyTree(cur.fetchall())


def get_entries(uri: str) -> list:
    con, cur = dbms.connect(uri)

    # InterPro entries (+ description)
    cur.execute(
        """
        SELECT
          E.ENTRY_AC, LOWER(ET.ABBREV), E.NAME, E.SHORT_NAME,
          E.CREATED, E.CHECKED, CA.TEXT
        FROM INTERPRO.ENTRY E
        INNER JOIN INTERPRO.CV_ENTRY_TYPE ET
          ON E.ENTRY_TYPE = ET.CODE
        LEFT OUTER JOIN INTERPRO.ENTRY2COMMON E2C
          ON E.ENTRY_AC = E2C.ENTRY_AC
        LEFT OUTER JOIN INTERPRO.COMMON_ANNOTATION CA
          ON E2C.ANN_ID = CA.ANN_ID
        """
    )

    entries = {}
    for row in cur:
        entry_ac = row[0]

        if entry_ac in entries:
            e = entries[entry_ac]
        else:
            e = entries[entry_ac] = {
                "accession": entry_ac,
                "type": row[1],
                "name": row[2],
                "short_name": row[3],
                "database": "interpro",
                "date": row[4],
                "is_checked": row[5] == 'Y',
                "descriptions": [],
                "integrated": None,
                "member_databases": {},
                "go_terms": [],
                "hierarchy": {},
                "citations": {},
                "cross_references": {}
            }

        if row[6] is not None:
            # todo: formatting descriptions
            """
            Some annotations contain multiple blocks of text,
            but some blocks might not be surrounded by <p> and </p> tags.

            Other blocks miss the opening <p> tag
            but have the closing </p> tag (or reverse).
            """
            e["descriptions"].append(row[6])

    # InterPro entry contributing signatures
    cur.execute(
        """
        SELECT
          EM.ENTRY_AC, M.METHOD_AC,
          LOWER(DB.DBSHORT), M.NAME, M.DESCRIPTION
        FROM INTERPRO.ENTRY2METHOD EM
        INNER JOIN INTERPRO.METHOD M
          ON EM.METHOD_AC = M.METHOD_AC
        INNER JOIN INTERPRO.CV_DATABASE DB
          ON M.DBCODE = DB.DBCODE
        """
    )

    for row in cur:
        entry_ac = row[0]
        method_ac = row[1]
        method_db = row[2]
        method_name = row[3]
        method_descr = row[4]

        if entry_ac not in entries:
            continue

        databases = entries[entry_ac]["member_databases"]

        if method_db not in databases:
            databases[method_db] = {}

        if method_descr:
            databases[method_db][method_ac] = method_descr
        else:
            databases[method_db][method_ac] = method_name

    # Only keep InterPro entries with contributing signatures
    entries = {
        entry_ac: entry
        for entry_ac, entry in entries.items()
        if entry["member_databases"]
    }

    # GO terms (InterPro entries)
    for entry_ac, go_id, name, cat_code, cat_name in goa.get_terms(cur):
        if entry_ac in entries:
            entries[entry_ac]["go_terms"].append({
                "identifier": go_id,
                "name": name,
                "category": {
                    "code": cat_code,
                    "name": cat_name
                }
            })

    # Hierarchy (InterPro entries)
    hierarchy = get_relationships(cur)

    for entry_ac in entries:
        accession = hierarchy.get_root(entry_ac)
        entries[entry_ac]["hierarchy"] = format_node(hierarchy, entries,
                                                     accession)

    # Member database entries (with literature references, and integration)
    methods = {}
    cur.execute(
        """
        SELECT
          M.METHOD_AC, M.NAME, LOWER(ET.ABBREV), M.DESCRIPTION,
          LOWER(DB.DBSHORT), M.ABSTRACT, M.ABSTRACT_LONG, M.METHOD_DATE,
          E2M.ENTRY_AC
        FROM INTERPRO.METHOD M
        INNER JOIN INTERPRO.CV_ENTRY_TYPE ET
          ON M.SIG_TYPE = ET.CODE
        INNER JOIN INTERPRO.CV_DATABASE DB
          ON M.DBCODE = DB.DBCODE
        LEFT OUTER JOIN INTERPRO.ENTRY2METHOD E2M
          ON M.METHOD_AC = E2M.METHOD_AC
        WHERE M.DBCODE != 'D'
        """
    )

    for row in cur:
        method_ac = row[0]

        abstr = row[5]
        abstr_long = row[6]

        if abstr is not None:
            descr = [abstr.lstrip("<p>").rstrip("</p>")]
        elif abstr_long is not None:
            descr = [abstr_long.read().lstrip("<p>").rstrip("</p>")]
        else:
            descr = []

        entry_ac = row[8]
        if entry_ac and not entries[entry_ac]["is_checked"]:
            entry_ac = None

        methods[method_ac] = {
            "accession": method_ac,
            "type": row[2],
            "name": row[3],
            "short_name": row[1],
            "database": row[4],
            "date": row[7],
            "descriptions": descr,
            "integrated": entry_ac,
            "member_databases": {},
            "go_terms": [],
            "hierarchy": {},
            "citations": {},
            "cross_references": {}
        }

    # Merging Interpro entries and Member DB entries
    entries.update(methods)

    # References
    citations = {}
    entries2citations = {}
    cur.execute(
        """
        SELECT
          E.ENTRY_AC, C.PUB_ID, C.PUBMED_ID, C.ISBN, C.VOLUME, C.ISSUE,
          C.YEAR, C.TITLE, C.URL, C.RAWPAGES, C.MEDLINE_JOURNAL,
          C.ISO_JOURNAL, C.AUTHORS, C.DOI_URL
        FROM (
          SELECT ENTRY_AC, PUB_ID
          FROM INTERPRO.ENTRY2PUB
          UNION ALL
          SELECT ENTRY_AC, PUB_ID
          FROM INTERPRO.PDB_PUB_ADDITIONAL
          UNION ALL
          SELECT ENTRY_AC, PUB_ID
          FROM INTERPRO.SUPPLEMENTARY_REF
          UNION ALL
          SELECT METHOD_AC AS ENTRY_AC, PUB_ID
          FROM INTERPRO.METHOD2PUB
        ) E
        INNER JOIN INTERPRO.CITATION C ON E.PUB_ID = C.PUB_ID
        """
    )

    for row in cur:
        entry_ac = row[0]
        if entry_ac not in entries:
            continue
        elif entry_ac not in entries2citations:
            entries2citations[entry_ac] = set()

        pub_id = row[1]
        entries2citations[entry_ac].add(pub_id)

        if pub_id not in citations:
            if row[12] is None:
                authors = []
            else:
                authors = [name.strip() for name in row[12].split(",")]

            citations[pub_id] = {
                "authors": authors,
                "DOI_URL": row[13],
                "ISBN": row[3],
                "issue": row[5],
                "ISO_journal": row[11],
                "medline_journal": row[10],
                "raw_pages": row[9],
                "PMID": row[2],
                "title": row[7],
                "URL": row[8],
                "volume": row[4],
                "year": row[6]
            }

    for entry_ac in entries2citations:
        for pub_id in entries2citations[entry_ac]:
            entries[entry_ac]['citations'][pub_id] = citations[pub_id]

    """
    Cross-references (InterPro entries only)
    Exclude the following databases:
        * C: PANDIT (outdated)
        * E: MSDsite (incorporated in PDB)
        * b: PDB (structures accessible from the "Structures" tab)
        * L: Blocks (outdated)
    """
    cur.execute(
        """
        SELECT X.ENTRY_AC, LOWER(D.DBSHORT), X.AC
        FROM INTERPRO.ENTRY_XREF X
        INNER JOIN INTERPRO.CV_DATABASE D ON X.DBCODE = D.DBCODE
        WHERE X.DBCODE NOT IN ('C', 'E', 'b', 'L')
        """
    )

    for entry_ac, xref_db, xref_ac in cur:
        if entry_ac not in entries:
            continue

        cross_references = entries[entry_ac]["cross_references"]

        if xref_db not in cross_references:
            cross_references[xref_db] = []

        cross_references[xref_db].append(xref_ac)

    cur.close()
    con.close()

    # Remove entries with a "is_checked" property (InterPro) set to False
    return [e for e in entries.values() if e.get("is_checked", True)]


def get_profile_alignments2(uri: str, threshold: float=1e-2,
                            chunk_size: int=1000, processes: int=1,
                            sync_frequency: int=1000000,
                            tmpdir: Optional[str]=None):
    con, cur = dbms.connect(uri)
    cur.execute(
        """
        SELECT S.SET_AC, LENGTH(S.SEQUENCE), DB.DBSHORT, S.METHOD_AC
        FROM INTERPRO.METHOD_SET S
        INNER JOIN INTERPRO.CV_DATABASE DB
            ON S.DBCODE = DB.DBCODE
        WHERE S.SET_AC IS NOT NULL
        """
    )

    sets = {}
    set_extra = {}
    entry2set = {}
    for set_ac, length, database, entry_ac in cur:
        entry2set[entry_ac] = set_ac
        if set_ac in sets:
            sets[set_ac].append({
                "accession": entry_ac,
                "type": "entry",
                "score": 1
            })
        else:
            sets[set_ac] = [{
                "accession": entry_ac,
                "type": "entry",
                "score": 1
            }]
            set_extra[set_ac] = (length, database)

    keys = sorted(entry2set.keys())
    keys = [keys[i] for i in range(0, len(keys), chunk_size)]
    with io.Store(None, keys, tmpdir) as store:
        cur.execute(
            """
            SELECT QUERY_AC, TARGET_AC, EVALUE, EVALUE_STR, DOMAINS
            FROM INTERPRO.METHOD_SCAN
            """
        )

        i = 0
        for row in cur:
            query_ac = row[0]
            if query_ac in entry2set:
                query_set_ac = entry2set[query_ac]
            else:
                continue

            target_ac = row[1]
            if row[2]:
                evalue = row[2]
            else:
                # evalue (BINARY_DOUBLE) == 0: check the evalue stored as string
                if row[3] == "0":
                    # Zero as well: take the smallest possible value
                    evalue = sys.float_info.min
                else:
                    # Somehow the BINARY_DOUBLE was rounded to 0...
                    evalue = float(row[3])

            if evalue > threshold:
                continue

            domains = []
            # DOMAINS is of type CLOB, hence the .read()
            for dom in json.loads(row[4].read()):
                domains.append({
                    "start": dom["start"],
                    "end": dom["end"],
                })
            domains.sort(key=repr_frag)

            length, database = set_extra[query_set_ac]
            target_set_ac = entry2set.get(target_ac)

            store.update(
                query_ac,
                {
                    "database": database,
                    "length": length,
                    "targets": [{
                        "accession": target_ac,
                        "set": target_set_ac,
                        "score": evalue,
                        "domains": domains
                    }]
                }
            )

            i += 1
            if sync_frequency and not i % sync_frequency:
                store.sync()

        cur.close()
        con.close()
        store.merge(processes=processes)


def get_profile_alignments(uri: str, database: str,
                           threshold: float=1e-2) -> Generator[Tuple, None, None]:
    con, cur = dbms.connect(uri)

    # Get sets and their members
    cur.execute(
        """
        SELECT S.SET_AC, S.METHOD_AC
        FROM INTERPRO.METHOD_SET S
        INNER JOIN INTERPRO.CV_DATABASE DB
            ON S.DBCODE = DB.DBCODE
        WHERE S.SET_AC IS NOT NULL
        AND DB.DBSHORT = UPPER(:1)
        """,
        (database,)
    )

    sets = {}
    entry2set = {}
    for set_ac, entry_ac in cur:
        entry2set[entry_ac] = set_ac
        if set_ac in sets:
            sets[set_ac].append({
                "accession": entry_ac,
                "type": "entry",
                "score": 1
            })
        else:
            sets[set_ac] = [{
                "accession": entry_ac,
                "type": "entry",
                "score": 1
            }]

    """
    Then get alignments between signatures
    Ordering by SET_AC allow us not to map everything in memory
    """
    cur.execute(
        """
        SELECT
          S.SET_AC, LENGTH(S.SEQUENCE), A.QUERY_AC,
          A.TARGET_AC, A.EVALUE, A.EVALUE_STR, A.DOMAINS
        FROM INTERPRO.METHOD_SET S
        INNER JOIN INTERPRO.CV_DATABASE DB
            ON S.DBCODE = DB.DBCODE
        INNER JOIN INTERPRO.METHOD_SCAN A
            ON S.METHOD_AC = A.QUERY_AC
        WHERE S.SET_AC IS NOT NULL
        AND DB.DBSHORT = :1
        AND A.EVALUE <= :2
        ORDER BY S.SET_AC
        """,
        (database.upper(), threshold)
    )

    _set_acc = None
    links = {}
    alignments = {}
    for row in cur:
        set_acc = row[0]
        if set_acc != _set_acc:
            if _set_acc:
                relationships = {
                    "nodes": sets[_set_acc],
                    "links": []
                }

                for query_acc, targets in links.items():
                    for target_acc, evalue in targets.items():
                        relationships["links"].append({
                            "source": query_acc,
                            "target": target_acc,
                            "score": evalue
                        })

                yield _set_acc, relationships, alignments

            _set_acc = set_acc
            links = {}
            alignments = {}

        length = row[1]
        query_acc = row[2]
        target_acc = row[3]
        evalue = row[4]
        if not evalue:
            # evalue (BINARY_DOUBLE) == 0: check the evalue stored as string
            if row[5] == "0":
                # Zero as well: take the smallest possible value
                evalue = sys.float_info.min
            else:
                # Somehow the BINARY_DOUBLE was rounded to 0...
                evalue = float(row[5])

        target_set_acc = entry2set.get(target_acc)
        domains = []
        # DOMAINS is of type CLOB, hence the .read()
        for dom in json.loads(row[6].read()):
            domains.append({
                "start": dom["start"],
                "end": dom["end"],
            })
        domains.sort(key=lambda x: x["start"])

        if query_acc in alignments:
            targets = alignments[query_acc]
        else:
            targets = alignments[query_acc] = {}

        targets[target_acc] = {
            "set_acc": target_set_acc,
            "score": evalue,
            "length": length,
            "domains": domains
        }

        if set_acc == target_set_acc:
            # Query and target belong to the same set
            # Keep only one edge, and the smallest e-value
            if query_acc > target_acc:
                query_acc, target_acc = target_acc, query_acc

            if query_acc not in links:
                links[query_acc] = {target_acc: evalue}
            elif target_acc not in links[query_acc]:
                links[query_acc][target_acc] = evalue
            elif evalue < links[query_acc][target_acc]:
                links[query_acc][target_acc] = evalue

    cur.close()
    con.close()

    if _set_acc:
        relationships = {
            "nodes": sets[_set_acc],
            "links": []
        }

        for query_acc, targets in links.items():
            for target_acc, evalue in targets.items():
                relationships["links"].append({
                    "source": query_acc,
                    "target": target_acc,
                    "score": evalue
                })

        yield _set_acc, relationships, alignments


def get_taxa(uri):
    con, cur = dbms.connect(uri)
    cur.execute(
        """
        SELECT
            TO_CHAR(TAX_ID), TO_CHAR(PARENT_ID), SCIENTIFIC_NAME,
            FULL_NAME, RANK
        FROM INTERPRO.ETAXI
        """
    )

    taxa = {}
    for row in cur:
        tax_id = row[0]

        taxa[tax_id] = {
            'id': tax_id,
            'parent_id': row[1],
            'sci_name': row[2],
            'full_name': row[3],
            'rank': row[4],
            'lineage': [tax_id],
            'children': set()
        }

    cur.close()
    con.close()

    for tax_id, taxon in taxa.items():
        child_id = tax_id
        parent_id = taxon['parent_id']

        while parent_id is not None:
            parent = taxa[parent_id]
            taxon['lineage'].append(parent['id'])
            parent['children'].add(child_id)

            child_id = parent_id
            parent_id = parent['parent_id']

    # taxa with short lineage first
    taxa = sorted(taxa.values(), key=lambda t: len(t['lineage']))

    for taxon in taxa:
        taxon['lineage'] = list(map(str, reversed(taxon['lineage'])))
        taxon['children'] = list(taxon['children'])

    return taxa


def get_structural_predictions(uri: str) -> dict:
    con, cur = dbms.connect(uri)

    # Only ModBase and SWISS-MODEL matches
    cur.execute(
        """
        SELECT
          M.PROTEIN_AC,
          LOWER(D.DBSHORT),
          M.DOMAIN_ID,
          S.FAM_ID,
          M.POS_FROM,
          M.POS_TO
        FROM INTERPRO.MATCH_STRUCT M
        INNER JOIN INTERPRO.STRUCT_CLASS S ON M.DOMAIN_ID = S.DOMAIN_ID
        INNER JOIN INTERPRO.CV_DATABASE D ON M.DBCODE = D.DBCODE
        WHERE M.DBCODE IN ('A', 'W')
        """
    )

    proteins = {}
    for acc, database, domain_id, fam_id, start, end in cur:
        if acc in proteins:
            p = proteins[acc]
        else:
            p = proteins[acc] = {}

        if database in p:
            db = p[database]
        else:
            db = p[database] = {}

        if domain_id in db:
            dom = db[domain_id]
        else:
            dom = db[domain_id] = {
                "class_id": domain_id,
                "domain_id": fam_id,
                "coordinates": []
            }

        dom["coordinates"].append({"start": start, "end": end})

    cur.close()
    con.close()

    for p in proteins.values():
        for db in p.values():
            for dom in db.values():
                dom["coordinates"].sort(key=lambda x: (x["start"], x["end"]))

    return proteins


def get_isoforms(uri: str) -> list:
    con, cur = dbms.connect(uri)
    cur.execute(
        """
        SELECT
          XV.PROTEIN_AC,
          XV.VARIANT,
          P.LEN,
          P.SEQ_SHORT,
          P.SEQ_LONG,
          MA.METHOD_AC,
          MA.SEQ_START,
          MA.SEQ_END,
          MA.FRAGMENTS,
          MA.MODEL_AC
        FROM (
          SELECT
            SUBSTR(AC, 1, INSTR(AC, '-') - 1) AS PROTEIN_AC,
            SUBSTR(AC, INSTR(AC, '-') + 1) AS VARIANT,
            UPI
          FROM UNIPARC.XREF
          WHERE DBID IN (24, 25)
            AND DELETED = 'N'
        ) XV
        INNER JOIN UNIPARC.PROTEIN P
          ON XV.UPI = P.UPI
        INNER JOIN UNIPARC.XREF XP
          ON XV.PROTEIN_AC = XP.AC
            AND XP.DBID IN (2, 3)
            AND XP.DELETED = 'N'
        INNER JOIN IPRSCAN.MV_IPRSCAN MA
          ON XV.UPI = MA.UPI
        INNER JOIN INTERPRO.IPRSCAN2DBCODE I2D
          ON MA.ANALYSIS_ID = I2D.IPRSCAN_SIG_LIB_REL_ID
        INNER JOIN INTERPRO.METHOD ME
          ON MA.METHOD_AC = ME.METHOD_AC
        WHERE XV.UPI != XP.UPI
        """
    )

    isoforms = {}
    for row in cur:
        accession = row[0] + '-' + row[1]
        if accession in isoforms:
            isoform = isoforms[accession]
        else:
            isoform = isoforms[accession] = {
                "accession": accession,
                "protein_acc": row[0],
                "length": row[2],
                "sequence": row[3] if row[3] else row[4].read(),
                "features": {}
            }

        if row[8] is None:
            fragments = [{
                "start": row[6],
                "end": row[7],
                "dc-status": "CONTINUOUS"
            }]
        else:
            fragments = []
            for frag in row[8].split(','):
                # Format: START-END-STATUS
                s, e, t = frag.split('-')
                fragments.append({
                    "start": int(s),
                    "end": int(e),
                    "dc-status": DC_STATUSES[t]
                })
            fragments.sort(key=repr_frag)

        signature_acc = row[5]
        if signature_acc in isoform["features"]:
            isoform["features"][signature_acc].append({
                "fragments": fragments,
                "model_acc": row[9]
            })
        else:
            isoform["features"][signature_acc] = [{
                "fragments": fragments,
                "model_acc": row[9] if row[9] != signature_acc else None
            }]

    cur.close()
    con.close()

    for isoform in isoforms.values():
        for locations in isoform["features"].values():
            locations.sort(key=lambda x: repr_frag(x["fragments"][0]))

    return list(isoforms.values())
