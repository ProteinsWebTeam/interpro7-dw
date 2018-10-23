#!/usr/bin/env python
# -*- coding: utf-8 -*-

import base64
import gzip
import json
import logging
import os
import re
import tempfile
import urllib.parse
import urllib.error
import urllib.request

from i7dw import dbms
from i7dw.disk import Store
from i7dw.ebi import hmmer

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s: %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)


def chunk_proteins(uri, dst, chunk_size=200000):
    con, cur = dbms.connect(uri)
    cur.execute(
        """
        SELECT PROTEIN_AC 
        FROM INTERPRO.PROTEIN
        ORDER BY PROTEIN_AC
        """
    )
    
    cnt = 0
    chunks = []
    for row in cur:
        cnt += 1
        if cnt % chunk_size == 1:
            chunks.append(row[0])
        
    cur.close()
    con.close()
        
    with open(dst, "wt") as fh:
        json.dump(chunks, fh)


def export_protein2structures(uri, src, dst, tmpdir=None, flush=1000000):
    logging.info("starting")

    with open(src, "rt") as fh:
        keys = json.load(fh)

    s = Store(dst, keys, tmpdir)
    con, cur = dbms.connect(uri)
    cur.execute(
        """
        SELECT
          M.PROTEIN_AC,
          LOWER(D.DBNAME),
          LOWER(M.DOMAIN_ID),
          LOWER(S.FAM_ID),
          M.POS_FROM,
          M.POS_TO
        FROM INTERPRO.MATCH_STRUCT M
        INNER JOIN INTERPRO.STRUCT_CLASS S ON M.DOMAIN_ID = S.DOMAIN_ID
        INNER JOIN INTERPRO.CV_DATABASE D ON M.DBCODE = D.DBCODE
        """
    )

    i = 0
    for acc, database, domain_id, fam_id, start, end in cur:
        if database in ('pdb', 'scop', 'cath'):
            k = 'feature'
        elif database in ('modbase', 'swiss-model'):
            k = 'prediction'
        else:
            continue

        s.update(
            acc,
            {
                k: {
                    database: {
                        domain_id: {
                            "class_id": domain_id,
                            "domain_id": fam_id,
                            "coordinates": [{"start": start, "end": end}]
                        }
                    }
                }
            }
        )

        i += 1
        if not i % flush:
            s.flush()

        if not i % 1000000:
            logging.info("{:>12}".format(i))

    cur.close()
    con.close()
    logging.info("{:>12}".format(i))
    size = s.merge(func=sort_struct_coordinates)
    logging.info("temporary files: {} bytes".format(size))


def sort_struct_coordinates(item: dict) -> dict:
    for databases in item.values():
        for db in databases.values():
            for domain in db.values():
                domain["coordinates"].sort(key=lambda x: (x["start"],
                                                          x["end"]))
    return item


def export_protein2matches(uri, src, dst, tmpdir=None, flush=1000000):
    logging.info("starting")

    with open(src, "rt") as fh:
        keys = json.load(fh)

    s = Store(dst, keys, tmpdir)
    con, cur = dbms.connect(uri)
    cur.execute(
        """
        SELECT 
          M.PROTEIN_AC PROTEIN_AC, LOWER(M.METHOD_AC), M.MODEL_AC, NULL,
          LOWER(E2M.ENTRY_AC), E.CHECKED, M.POS_FROM, M.POS_TO, M.FRAGMENTS
        FROM INTERPRO.MATCH M
        LEFT OUTER JOIN INTERPRO.ENTRY2METHOD E2M 
          ON M.METHOD_AC = E2M.METHOD_AC
        LEFT OUTER JOIN INTERPRO.ENTRY E 
          ON E2M.ENTRY_AC = E.ENTRY_AC
        AND M.STATUS = 'T'
        AND M.POS_FROM IS NOT NULL
        AND M.POS_TO IS NOT NULL   
        UNION ALL
        SELECT 
          FM.PROTEIN_AC PROTEIN_AC, LOWER(FM.METHOD_AC), NULL, FM.SEQ_FEATURE,
          NULL, NULL, FM.POS_FROM, FM.POS_TO, NULL
        FROM INTERPRO.FEATURE_MATCH FM
        WHERE FM.DBCODE = 'g'
        """
    )

    i = 0
    for row in cur:
        protein_acc = row[0]
        method_acc = row[1]
        model_acc = row[2]
        seq_feature = row[3]
        entry_acc = row[4]
        is_checked = row[5] == 'Y'
        pos_start = row[6]
        pos_end = row[7]
        fragments_str = row[8]

        if fragments_str is None:
            fragments = [{"start": pos_start, "end": pos_end}]
        else:
            fragments = []
            for frag in fragments_str.split(','):
                """
                Format: START-END-TYPE 
                Types:
                    * S: Continuous single chain domain
                    * N: N-terminal discontinuous
                    * C: C-terminal discontinuous
                    * NC: N and C -terminal discontinuous
                """
                pos_start, pos_end, t = frag.split('-')
                fragments.append({
                    "start": int(pos_start),
                    "end": int(pos_end)
                })

        s.append(protein_acc, {
            "method_ac": method_acc,
            "model_ac": model_acc,
            "seq_feature": seq_feature,
            "entry_ac": entry_acc if is_checked else None,
            "fragments": fragments
        })

        i += 1
        if not i % flush:
            s.flush()

        if not i % 1000000:
            logging.info("{:>12}".format(i))

    cur.close()
    con.close()
    logging.info("{:>12}".format(i))
    size = s.merge(func=sort_matches)
    logging.info("temporary files: {} bytes".format(size))


def sort_fragments(fragments: list) -> tuple:
    start = end = None
    for f in fragments:
        if start is None or f["start"] < start:
            start = f["start"]

        if end is None or f["end"] < end:
            end = f["end"]

    return start, end


def sort_matches(matches: list) -> list:
    return sorted(matches, key=lambda m: sort_fragments(m["fragments"]))


def export_protein2features(uri, src, dst, tmpdir=None, flush=1000000):
    logging.info("starting")

    with open(src, "rt") as fh:
        keys = json.load(fh)

    s = Store(dst, keys, tmpdir)
    con, cur = dbms.connect(uri)
    cur.execute(
        """
        SELECT 
          FM.PROTEIN_AC, LOWER(FM.METHOD_AC), LOWER(DB.DBSHORT), 
          FM.POS_FROM, FM.POS_TO
        FROM INTERPRO.FEATURE_MATCH FM
        INNER JOIN INTERPRO.CV_DATABASE DB ON FM.DBCODE = DB.DBCODE
        WHERE FM.DBCODE != 'g'
        """
    )

    i = 0
    for protein_acc, method_acc, database, start, end in cur:
        s.update(
            protein_acc,
            {
                method_acc: {
                    "accession": method_acc,
                    "source_database": database,
                    "locations": [{"start": start, "end": end}]
                }
            }
        )

        i += 1
        if not i % flush:
            s.flush()

        if not i % 1000000:
            logging.info("{:>12}".format(i))

    cur.close()
    con.close()
    logging.info("{:>12}".format(i))
    size = s.merge(func=sort_feature_locations)
    logging.info("temporary files: {} bytes".format(size))


def sort_feature_locations(item: dict) -> dict:
    for method in item.values():
        method["locations"].sort(key=lambda x: (x["start"], x["end"]))
    return item


def export_protein2residues(uri, src, dst, tmpdir=None, flush=1000000):
    logging.info("starting")

    with open(src, "rt") as fh:
        keys = json.load(fh)

    s = Store(dst, keys, tmpdir)
    con, cur = dbms.connect(uri)
    cur.execute(
        """
        SELECT 
          S.PROTEIN_AC, LOWER(S.METHOD_AC), M.NAME, LOWER(D.DBSHORT), 
          S.DESCRIPTION, S.RESIDUE, S.RESIDUE_START, S.RESIDUE_END
        FROM INTERPRO.SITE_MATCH S
        INNER JOIN INTERPRO.METHOD M ON S.METHOD_AC = M.METHOD_AC
        INNER JOIN INTERPRO.CV_DATABASE D ON M.DBCODE = D.DBCODE
        """
    )

    i = 0
    for row in cur:
        protein_acc = row[0]
        method_acc = row[1]
        method_name = row[2]
        database = row[3]
        description = row[4]
        residue = row[5]
        start = row[6]
        end = row[7]

        s.update(
            protein_acc,
            {
                method_acc: {
                    "accession": method_acc,
                    "name": method_name,
                    "source_database": database,
                    "locations": {
                        description: {
                            "description": description,
                            "fragments": [{
                                "residues": residue,
                                "start": start,
                                "end": end
                            }]
                        }
                    }
                }
            }
        )

        i += 1
        if not i % flush:
            s.flush()

        if not i % 1000000:
            logging.info("{:>12}".format(i))

    cur.close()
    con.close()
    logging.info("{:>12}".format(i))
    size = s.merge(func=sort_residues)
    logging.info("temporary files: {} bytes".format(size))


def sort_residues(item: dict) -> dict:
    for method in item.values():
        locations = []
        for loc in method["locations"].values():
            locations.append(sorted(loc["fragments"],
                                    key=lambda x: (x["start"], x["end"])))

        method["locations"] = locations

    return item


def export_proteins(uri, src, dst_proteins, dst_sequences,
                    tmpdir=None, flush=1000000):
    logging.info("starting")

    with open(src, "rt") as fh:
        keys = json.load(fh)

    proteins = Store(dst_proteins, keys, tmpdir)
    sequences = Store(dst_sequences, keys, tmpdir)
    con, cur = dbms.connect(uri)
    cur.execute(
        """
        SELECT
          IP.PROTEIN_AC,
          TO_CHAR(IP.TAX_ID),
          IP.NAME,
          IP.DBCODE,
          IP.FRAGMENT,
          IP.LEN,
          UP.SEQ_SHORT,
          UP.SEQ_LONG
        FROM INTERPRO.PROTEIN IP
        INNER JOIN UNIPARC.XREF UX 
          ON IP.PROTEIN_AC = UX.AC AND UX.DELETED = 'N'
        INNER JOIN UNIPARC.PROTEIN UP 
          ON UX.UPI = UP.UPI
        """
    )

    i = 0
    for acc, tax_id, name, dbcode, frag, length, seq_short, seq_long in cur:
        proteins[acc] = {
            "taxon": tax_id,
            "identifier": name,
            "isReviewed": dbcode == 'S',
            "isFrag": frag == 'Y',
            "length": length
        }

        if seq_long is not None:
            sequences[acc] = seq_long.read()
        else:
            sequences[acc] = seq_short

        i += 1
        if not i % flush:
            proteins.flush()
            sequences.flush()

        if not i % 1000000:
            logging.info("{:>12}".format(i))

    cur.close()
    con.close()
    logging.info("{:>12}".format(i))
    size1 = proteins.merge()
    size2 = sequences.merge()
    logging.info("temporary files: "
                 "{} bytes (proteins), "
                 "{} bytes (sequences)".format(size1, size2))


def get_taxa(uri):
    con, cur = dbms.connect(uri)
    cur.execute(
        """
        SELECT TO_CHAR(TAX_ID), TO_CHAR(PARENT_ID), SCIENTIFIC_NAME, FULL_NAME, RANK, LEFT_NUMBER, RIGHT_NUMBER
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
            'left_number': row[5],
            'right_number': row[6],
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


def get_databases(uri):
    # todo: do not hardcode this value!
    member_dbs = {
        'B', 'D', 'F', 'H', 'I', 'J', 'M', 'N', 'P', 'Q', 'R', 'U', 'V',
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
        name_long = row[2]
        description = row[3]
        version = row[4]
        release_date = row[5]
        prev_version = row[6]
        prev_release_date = row[7]

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

        databases.append((
            name,
            name_long,
            description,
            version,
            release_date,
            db_type,
            prev_version,
            prev_release_date
        ))

    cur.close()
    con.close()

    return databases


class _EntryHierarchyTree(object):
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


def _format_node(hierarchy, entries, accession):
    entry = entries[accession]
    children = hierarchy.children_of(accession)
    return {
        'accession': accession,
        'name': entry['name'],
        'type': entry['type'],
        'children': [_format_node(hierarchy, entries, child_ac) for child_ac in children]
    }


def _get_relationships(cur):
    cur.execute(
        """
        SELECT LOWER(ENTRY_AC), LOWER(PARENT_AC)
        FROM INTERPRO.ENTRY2ENTRY
        WHERE RELATION = 'TY'
        """
    )

    return _EntryHierarchyTree(cur.fetchall())


def get_entries(uri: str) -> list:
    con, cur = dbms.connect(uri)

    # InterPro entries (+ description)
    cur.execute(
        """
        SELECT 
          LOWER(E.ENTRY_AC), LOWER(ET.ABBREV), E.NAME, E.SHORT_NAME, 
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
                'accession': entry_ac,
                'type': row[1],
                'name': row[2],
                'short_name': row[3],
                'database': 'interpro',
                'date': row[4],
                'is_checked': row[5] == 'Y',
                'descriptions': [],
                'integrated': None,
                'member_databases': {},
                'go_terms': [],
                'hierarchy': {},
                'citations': {},
                'cross_references': {}
            }

        if row[6] is not None:
            # todo: formatting descriptions
            '''
            Some annotations contain multiple blocks of text, 
            but some blocks might not be surrounded by <p> and </p> tags.
            
            Other blocks miss the opening <p> tag 
            but have the closing </p> tag (or reverse).
            
            Shit's broken, yo.
            '''
            e['descriptions'].append(row[6])

    # InterPro entry contributing signatures
    cur.execute(
        """
        SELECT 
          LOWER(EM.ENTRY_AC), LOWER(M.METHOD_AC), 
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

        databases = entries[entry_ac]['member_databases']

        if method_db not in databases:
            databases[method_db] = {}

        if method_descr:
            databases[method_db][method_ac] = method_descr
        else:
            databases[method_db][method_ac] = method_name

    # Remove InterPro entries without contributing signatures
    entries = {
        entry_ac: entry
        for entry_ac, entry in entries.items()
        if not entry["member_databases"]
    }

    # GO terms (InterPro entries)
    cur.execute(
        """
        SELECT 
          LOWER(I2G.ENTRY_AC), GT.GO_ID, GT.NAME, GT.CATEGORY, GC.TERM_NAME
        FROM INTERPRO.INTERPRO2GO I2G
        INNER JOIN GO.TERMS@GOAPRO GT 
          ON I2G.GO_ID = GT.GO_ID
        INNER JOIN GO.CV_CATEGORIES@GOAPRO GC 
          ON GT.CATEGORY = GC.CODE
        """
    )

    for entry_ac, term_id, term_name, cat_code, cat_name in cur:
        if entry_ac in entries:
            entries[entry_ac]['go_terms'].append({
                'identifier': term_id,
                'name': term_name,
                'category': {
                    'code': cat_code,
                    'name': cat_name
                }
            })

    # Hierarchy (InterPro entries)
    hierarchy = _get_relationships(cur)

    for entry_ac in entries:
        accession = hierarchy.get_root(entry_ac)
        entries[entry_ac]['hierarchy'] = _format_node(hierarchy, entries,
                                                      accession)

    # Member database entries (with literature references, and integration)
    methods = {}
    cur.execute(
        """
        SELECT 
          LOWER(M.METHOD_AC), M.NAME, LOWER(ET.ABBREV), M.DESCRIPTION, 
          LOWER(DB.DBSHORT), M.ABSTRACT, M.ABSTRACT_LONG, M.METHOD_DATE, LOWER(E2M.ENTRY_AC)
        FROM INTERPRO.METHOD M
        INNER JOIN INTERPRO.CV_ENTRY_TYPE ET 
          ON M.SIG_TYPE = ET.CODE
        INNER JOIN INTERPRO.CV_DATABASE DB 
          ON M.DBCODE = DB.DBCODE
        LEFT OUTER JOIN INTERPRO.ENTRY2METHOD E2M 
          ON M.METHOD_AC = E2M.METHOD_AC
        """
    )

    for row in cur:
        method_ac = row[0]

        abstr = row[5]
        abstr_long = row[6]

        if abstr is not None:
            descr = [abstr.lstrip('<p>').rstrip('</p>')]
        elif abstr_long is not None:
            descr = [abstr_long.read().lstrip('<p>').rstrip('</p>')]
        else:
            descr = []

        entry_ac = row[8]
        if entry_ac and not entries[entry_ac]["is_checked"]:
            entry_ac = None

        methods[method_ac] = {
            'accession': method_ac,
            'type': row[2],
            'name': row[3],
            'short_name': row[1],
            'database': row[4],
            'date': row[7],
            'descriptions': descr,
            'integrated': entry_ac,
            'member_databases': {},
            'go_terms': [],
            'hierarchy': {},
            'citations': {},
            'cross_references': {}
        }

    # Merging Interpro entries and Member DB entries
    entries.update(methods)

    # References
    citations = {}
    entries2citations = {}
    cur.execute(
        """
        SELECT 
          LOWER(E.ENTRY_AC), C.PUB_ID, C.PUBMED_ID, C.ISBN, C.VOLUME, C.ISSUE, 
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
                authors = [name.strip() for name in row[12].split(',')]

            citations[pub_id] = {
                'authors': authors,
                'DOI_URL': row[13],
                'ISBN': row[3],
                'issue': row[5],
                'ISO_journal': row[11],
                'medline_journal': row[10],
                'raw_pages': row[9],
                'PMID': row[2],
                'title': row[7],
                'URL': row[8],
                'volume': row[4],
                'year': row[6]
            }

    for entry_ac in entries2citations:
        for pub_id in entries2citations[entry_ac]:
            entries[entry_ac]['citations'][pub_id] = citations[pub_id]

    """
    Cross-references (InterPro entries only)
    Exclude the following databases:
        * C: PANDIT (not updated)
        * E: MSDsite (incorporated in PDB)
        * b: PDB (structures accessible from the "Structures" tab)
        * L: Blocks (not updated)
    """
    cur.execute(
        """
        SELECT LOWER(X.ENTRY_AC), LOWER(D.DBSHORT), LOWER(X.AC)
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


def get_deleted_entries(uri: str) -> list:
    con, cur = dbms.connect(uri)

    cur.execute(
        """
        SELECT 
            LOWER(ENTRY_AC), ENTRY_TYPE, NAME, 
            SHORT_NAME, TIMESTAMP, CHECKED
        FROM INTERPRO.ENTRY_AUDIT
        WHERE ENTRY_AC NOT IN (
            SELECT ENTRY_AC 
            FROM INTERPRO.ENTRY 
            WHERE CHECKED='Y'
        )
        ORDER BY ENTRY_AC, TIMESTAMP
        """
    )

    entries = {}
    for row in cur:
        acc = row[0]
        if acc in entries:
            entries[acc].update({
                "type": row[1],
                "name": row[2],
                "short_name": row[3],
                "deletion_date": row[4]
            })
            if row[5] == 'Y':
                entries[acc]["was_public"] = True
        else:
            entries[acc] = {
                "accession": acc,
                "type": row[1],
                "name": row[2],
                "short_name": row[3],
                "database": "interpro",
                "creation_date": row[4],
                "deletion_date": None,
                "was_public": row[5] == 'Y'
            }

    cur.close()
    con.close()
    return [e for e in entries.values() if e["was_public"]]


def get_pfam_wiki(uri):
    base_url = 'https://en.wikipedia.org/api/rest_v1/page/summary/'

    # Pfam DB in LATIN1, with special charachters in Wikipedia title
    con, cur = dbms.connect(uri, encoding='latin1')
    cur.execute(
        """
        SELECT LOWER(p.pfamA_acc), w.title
        FROM pfamA_wiki p
        INNER JOIN wikipedia w ON p.auto_wiki = w.auto_wiki
        """
    )

    rows = cur.fetchall()
    cur.close()
    con.close()

    entries = {}
    for acc, title in rows:
        # cursor returns bytes instead of string due to latin1
        acc = acc.decode()
        title = title.decode()

        try:
            url = base_url + urllib.parse.quote(title)
            res = urllib.request.urlopen(url)
        except urllib.error.HTTPError as e:
            # Content can be retrieved with e.fp.read()
            continue
        else:
            obj = json.loads(res.read().decode('utf-8'))
            thumbnail = obj.get('thumbnail')

            if thumbnail:
                try:
                    filename, headers = urllib.request.urlretrieve(thumbnail['source'])
                except urllib.error.ContentTooShortError:
                    b64str = None
                else:
                    with open(filename, 'rb') as fh:
                        b64str = base64.b64encode(fh.read()).decode('utf-8')

                    os.unlink(filename)
            else:
                b64str = None

            entries[acc] = {
                'title': title,
                'extract': obj['extract'],
                'extract_html': obj['extract_html'],
                'thumbnail': b64str
            }

    return entries


def get_pfam_annotations(cur):
    cur.execute(
        """
        SELECT LOWER(a.pfamA_acc), a.alignment, h.hmm
        FROM alignment_and_tree a
        INNER JOIN pfamA_HMM h ON a.pfamA_acc = h.pfamA_acc AND a.type = 'seed'
        """
    )

    entries = {}
    for pfam_ac, aln, hmm in cur:
        if pfam_ac not in entries:
            entries[pfam_ac] = []

        if aln is not None:
            entries[pfam_ac].append({
                'type': 'alignment',
                'value': aln,
                'mime_type': 'application/octet-stream'
            })

        if hmm is not None:
            entries[pfam_ac].append({
                'type': 'hmm',
                'value': hmm,
                'mime_type': 'application/octet-stream'
            })

            fd, filename = tempfile.mkstemp()
            os.close(fd)

            with open(filename, 'wb') as fh:
                fh.write(hmm)

            # hmm_logo = call_skylign(filename)
            hmm_logo = hmmer.hmm_to_logo(filename, method='info_content_all', processing='hmm')
            os.unlink(filename)

            entries[pfam_ac].append({
                'type': 'logo',
                'value': json.dumps(hmm_logo),
                'mime_type': 'application/json'
            })

    annotations = []
    for pfam_ac in entries:
        for anno in entries[pfam_ac]:
            annotations.append((
                '{}--{}'.format(pfam_ac, anno['type']),
                pfam_ac,
                anno['type'],
                anno['value'],
                anno['mime_type']
            ))

    return annotations


def get_pfam_clans(cur):
    cur.execute(
        """
        SELECT LOWER(pfamA_acc), num_full
        FROM pfamA
        """
    )

    entries = dict(cur.fetchall())

    cur.execute(
        """
        SELECT LOWER(c.clan_acc), c.clan_id, c.clan_description, c.number_sequences, LOWER(m.pfamA_acc), l.pfamA_acc_2, l.evalue
        FROM clan c
        INNER JOIN clan_membership m ON c.clan_acc = m.clan_acc
        LEFT OUTER JOIN pfamA2pfamA_hhsearch l ON m.pfamA_acc = l.pfamA_acc_1
        INNER JOIN clan_membership m2 ON l.pfamA_acc_2 = m2.pfamA_acc and c.clan_acc = m2.clan_acc
        """
    )

    clans = {}
    sizes = {}
    nodes = {}
    links = {}
    for row in cur:
        clan_ac = row[0]

        if clan_ac not in clans:
            clans[clan_ac] = {
                'accession': clan_ac,
                'name': row[1],
                'description': row[2]
            }

            sizes[clan_ac] = row[3]
            nodes[clan_ac] = {}
            links[clan_ac] = []

        pfam_ac_1 = row[4]

        if pfam_ac_1 not in nodes[clan_ac]:
            nodes[clan_ac][pfam_ac_1] = {
                'accession': pfam_ac_1,
                'type': 'entry',
                'score': entries[pfam_ac_1] / sizes[clan_ac]
            }

        pfam_ac_2 = row[5]
        if pfam_ac_2:
            pfam_ac_2 = pfam_ac_2.lower()
            evalue = row[6]

            if pfam_ac_2 not in nodes[clan_ac]:
                nodes[clan_ac][pfam_ac_2] = {
                    'accession': pfam_ac_2,
                    'type': 'entry',
                    'score': entries[pfam_ac_2] / sizes[clan_ac]
                }

            links[clan_ac].append({
                'source': pfam_ac_1,
                'target': pfam_ac_2,
                'score': evalue
            })

    for clan_ac in clans:
        clans[clan_ac]['relationships'] = {
            'nodes': list(nodes[clan_ac].values()),
            'links': links[clan_ac]
        }

    return list(clans.values())


def get_cdd_superfamilies():
    filename, headers = urllib.request.urlretrieve('ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/cddid_all.tbl.gz')

    sets = {}
    nodes = {}
    parent_of = {}
    domains = {}
    with gzip.open(filename, 'rt') as fh:
        for line in fh:
            pssm_id, accession, short_name, descr, pssm_length = line.rstrip().split('\t')
            accession = accession.lower()
            descr = descr.lstrip('N/A. ')

            if re.match(r'cl\d+', accession):
                sets[accession] = {
                    'accession': accession,
                    'name': short_name,
                    'description': descr
                }
                nodes[accession] = set()
            elif re.match(r'cd\d+', accession):
                domains[accession] = {
                    'name': short_name,
                    'description': descr
                }

    os.unlink(filename)

    filename, headers = urllib.request.urlretrieve('ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/family_superfamily_links')
    with open(filename, 'rt') as fh:
        for line in fh:
            cd_id, cd_pssm_id, cl_id, cl_pssm_id = line.rstrip().split('\t')
            cd_id = cd_id.lower()
            cl_id = cl_id.lower()

            if re.match(r'cl\d+', cl_id) and re.match(r'cd\d+', cd_id):
                nodes[cl_id].add(cd_id)

    os.unlink(filename)

    filename, headers = urllib.request.urlretrieve('ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/cdtrack.txt')
    with open(filename, 'rt') as fh:
        for line in fh:
            line = line.rstrip()

            if not line or line[0] == '#':
                continue

            fields = line.split()
            ac = fields[0].lower()
            parent_ac = fields[3].lower()

            if parent_ac in domains:
                parent_of[ac] = parent_ac

    os.unlink(filename)

    final_sets = []
    for cl_id in sets:
        _nodes = []
        _links = []
        for cd_id in nodes[cl_id]:
            cd_id = cd_id.lower()
            _nodes.append({
                'accession': cd_id,
                'type': 'entry',
                'score': 1
            })

            if cd_id in parent_of:
                _links.append({
                    'source': parent_of[cd_id],
                    'target': cd_id,
                    'score': 1
                })

        if _nodes:
            sets[cl_id]['relationships'] = {
                'nodes': _nodes,
                'links': _links
            }

            final_sets.append(sets[cl_id])

    return final_sets


class Supermatch(object):
    def __init__(self, entry_ac, entry_root, start, end):
        self.entries = {(entry_ac, entry_root)}
        self.start = start
        self.end = end

    def in_same_hierarchy(self, other):
        for acc_1, root_1 in self.entries:
            for acc_2, root_2 in other.entries:
                if root_1 != root_2:
                    return False

        return True

    def merge_if_overlap(self, other, min_overlap):
        overlap = min(self.end, other.end) - max(self.start, other.start) + 1
        shortest = min(self.end - self.start, other.end - other.start) + 1

        if (overlap / shortest * 100) >= min_overlap:
            self.entries = self.entries | other.entries
            self.start = min(self.start, other.start)
            self.end = max(self.end, other.end)
            return self
        else:
            return None

    def format(self):
        return '{}:{:.0f}-{:.0f}'.format(
            self.format_entries(),
            self.start,
            self.end
        )

    def format_entries(self):
        return '&'.join(sorted([re.sub(r'IPR0*', '', accession) for accession in self.get_entries()]))

    def get_entries(self):
        return [entry_ac for entry_ac, entry_root in self.entries]

    def __eq__(self, other):
        return (
                isinstance(other, Supermatch) and
                self.start == other.start and
                self.end == other.end and
                self.entries == other.entries
        )


class SupermatchSet(object):
    def __init__(self, supermatch):
        self.supermatches = [supermatch]

    def add(self, candidate, min_overlap):
        if not self.supermatches[0].in_same_hierarchy(candidate):
            return False

        merged = None
        for sm in self.supermatches:
            merged = sm.merge_if_overlap(candidate, min_overlap)

            if merged is not None:
                break

        if merged is None:
            self.supermatches.append(candidate)
        else:
            # Merged supermatch: we now need to remove overlaps between the newly merged supermatch and others
            indexes_ok = set()

            while True:
                index = None

                for i, sm in enumerate(self.supermatches):
                    if sm == merged or i in indexes_ok:
                        continue

                    if merged.merge_if_overlap(sm, min_overlap):
                        # Overlap so merged, we now have to remove the merged supermatch (sm)
                        index = i
                        break
                    else:
                        # No overlap, might be skipped during next iteration
                        indexes_ok.add(i)

                if index is None:
                    # No move overlaps
                    break
                else:
                    self.supermatches.pop(index)

        return True


def merge_supermatches(supermatches, min_overlap=20):
    sets = []

    for sm in supermatches:
        in_set = False

        for s in sets:
            in_set = s.add(sm, min_overlap)

            if in_set:
                break

        if not in_set:
            sets.append(SupermatchSet(sm))

    return sets
