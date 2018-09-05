#!/usr/bin/env python
# -*- coding: utf-8 -*-

import base64
import gzip
import json
import logging
import os
import re
import tempfile
import time
import urllib.error
import urllib.request

from i7dw import dbms, disk
from i7dw.disk import Store
from i7dw.ebi import hmmer

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s: %(levelname)s: %(message)s',
    datefmt='%y-%m-%d %H:%M:%S'
)


def _sort_struct_matches(proteins):
    for acc in proteins:
        for k in proteins[acc]:  # `k` being `feature` or `prediction`
            for dbname in proteins[acc][k]:
                db = proteins[acc][k][dbname]
                for domain_id in db:
                    db[domain_id]['coordinates'].sort(key=lambda m: (m['start'], m['end']))


def export_struct_matches(uri, dst, chunk_size=1000000):
    logging.info('starting')
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
        ORDER BY M.PROTEIN_AC
        """
    )

    store = Store(dst)
    proteins = {}
    cnt = 0
    for acc, dbname, domain_id, fam_id, start, end in cur:
        if dbname in ('pdb', 'scop', 'cath'):
            k = 'feature'
        elif dbname in ('modbase', 'swiss-model'):
            k = 'prediction'
        else:
            continue

        if acc not in proteins:
            if len(proteins) == chunk_size:
                cnt += len(proteins)
                _sort_struct_matches(proteins)
                store.add(proteins)
                proteins = {}
                logging.info('{:>12}'.format(cnt))

            proteins[acc] = {
                'feature': {},
                'prediction': {}
            }

        p = proteins[acc][k]

        if dbname in p:
            db = p[dbname]
        else:
            db = p[dbname] = {}

        if domain_id not in db:
            db[domain_id] = {
                'class_id': domain_id,
                'domain_id': fam_id,
                'coordinates': []
            }

        db[domain_id]['coordinates'].append({'start': start, 'end': end})

    cur.close()
    con.close()

    cnt += len(proteins)
    _sort_struct_matches(proteins)
    store.add(proteins)
    store.close()

    logging.info('{:>12}'.format(cnt))


def export_prot_matches_extra(uri, dst, chunk_size=1000000):
    logging.info('starting')
    con, cur = dbms.connect(uri)
    cur.execute(
        """
        SELECT 
          FM.PROTEIN_AC, FM.METHOD_AC, DB.DBSHORT, FM.POS_FROM, FM.POS_TO
        FROM INTERPRO.FEATURE_MATCH FM
        INNER JOIN INTERPRO.CV_DATABASE DB ON FM.DBCODE = DB.DBCODE
        WHERE FM.DBCODE != 'g'
        ORDER BY PROTEIN_AC, FM.POS_FROM, FM.POS_TO
        """
    )

    store = Store(dst)
    proteins = {}
    cnt = 0
    for acc, method_ac, source_db, start, end in cur:
        if acc in proteins:
            p = proteins[acc]
        else:
            if len(proteins) == chunk_size:
                store.add(proteins)
                proteins = {}
                logging.info('{:>12}'.format(cnt))

            cnt += 1
            p = proteins[acc] = {}

        if method_ac in p:
            method = p[method_ac]
        else:
            method = p[method_ac] = {
                'accession': method_ac,
                'source_database': source_db,
                'locations': []
            }

        method['locations'].append({
            'start': start,
            'end': end
        })

    cur.close()
    con.close()

    store.add(proteins)
    store.close()

    logging.info('{:>12}'.format(cnt))


def _sort_prot_matches(proteins):
    return {acc: sorted(
        proteins[acc],
        key=lambda m: (
            min([f['start'] for f in m['fragments']]),
            min([f['end'] for f in m['fragments']])
        )
    ) for acc in proteins}


def export_prot_matches(uri, dst, chunk_size=1000000):
    logging.info('starting')
    con, cur = dbms.connect(uri)
    cur.execute(
        """
        SELECT 
          M.PROTEIN_AC PROTEIN_AC, LOWER(M.METHOD_AC), M.MODEL_AC, NULL,
          LOWER(E2M.ENTRY_AC), M.POS_FROM, M.POS_TO, M.FRAGMENTS
        FROM INTERPRO.MATCH M
        LEFT OUTER JOIN INTERPRO.ENTRY2METHOD E2M ON M.METHOD_AC = E2M.METHOD_AC
        WHERE M.STATUS = 'T'
        AND M.POS_FROM IS NOT NULL
        AND M.POS_TO IS NOT NULL   
        UNION ALL
        SELECT 
          FM.PROTEIN_AC PROTEIN_AC, LOWER(FM.METHOD_AC), NULL, FM.SEQ_FEATURE,
          NULL, FM.POS_FROM, FM.POS_TO, NULL
        FROM INTERPRO.FEATURE_MATCH FM
        WHERE FM.DBCODE = 'g'
        ORDER BY PROTEIN_AC   
        """
    )

    store = Store(dst)
    proteins = {}
    cnt = 0
    for row in cur:
        acc = row[0]

        if acc in proteins:
            p = proteins[acc]
        else:
            if len(proteins) == chunk_size:
                cnt += len(proteins)
                store.add(_sort_prot_matches(proteins))
                proteins = {}
                logging.info('{:>12}'.format(cnt))

            p = proteins[acc] = []

        if row[7] is not None:
            fragments = []
            for frag in row[7].split(','):
                """
                Format: START-END-TYPE 
                Types:
                    * S: Continuous single chain domain
                    * N: N-terminal discontinuous
                    * C: C-terminal discontinuous
                    * NC: N and C -terminal discontinuous
                """
                s, e, t = frag.split('-')
                fragments.append({'start': int(s), 'end': int(e)})
        else:
            fragments = [{'start': row[5], 'end': row[6]}]

        p.append({
            'method_ac': row[1],
            'model_ac': row[2],
            'seq_feature': row[3],
            'entry_ac': row[4],
            'fragments': fragments
        })

    cur.close()
    con.close()

    cnt += len(proteins)
    store.add(_sort_prot_matches(proteins))
    store.close()

    logging.info('{:>12}'.format(cnt))


def _sort_residues(proteins):
    for acc in proteins:
        for method_ac in proteins[acc]:
            locations = []

            for loc in proteins[acc][method_ac]['locations'].values():
                # Sort residues by position
                loc['fragments'].sort(key=lambda f: (f['start'], f['end']))
                locations.append(loc)

            # Sort locations by description
            proteins[acc][method_ac]['locations'] = sorted(locations, key=lambda l: l['description'])


def export_residues(uri, dst, chunk_size=100000):
    logging.info('starting')
    con, cur = dbms.connect(uri)
    cur.execute(
        """
        SELECT 
          S.PROTEIN_AC, LOWER(S.METHOD_AC), M.NAME, D.DBSHORT, 
          S.DESCRIPTION, S.RESIDUE, S.RESIDUE_START, S.RESIDUE_END
        FROM INTERPRO.SITE_MATCH S
        INNER JOIN INTERPRO.METHOD M ON S.METHOD_AC = M.METHOD_AC
        INNER JOIN INTERPRO.CV_DATABASE D ON M.DBCODE = D.DBCODE
        ORDER BY S.PROTEIN_AC    
        """
    )

    store = Store(dst)
    proteins = {}
    cnt = 0
    for row in cur:
        acc = row[0]

        if acc in proteins:
            p = proteins[acc]
        else:
            if len(proteins) == chunk_size:
                _sort_residues(proteins)
                store.add(proteins)
                proteins = {}

            p = proteins[acc] = {}
            cnt += 1

            if not cnt % 1000000:
                logging.info('{:>12}'.format(cnt))

        method_ac = row[1]
        if method_ac in p:
            method = p[method_ac]
        else:
            method = p[method_ac] = {
                'entry_accession': method_ac,
                'accession': row[2],
                'source_database': row[3],
                'locations': {}
            }

        descr = row[4]
        if descr in method['locations']:
            loc = method['locations'][descr]
        else:
            loc = method['locations'][descr] = {
                'description': descr,
                'fragments': []
            }

        loc['fragments'].append({
            'residues': row[5],
            'start': row[6],
            'end': row[7]
        })

    cur.close()
    con.close()

    _sort_residues(proteins)
    store.add(proteins)
    store.close()

    logging.info('{:>12}'.format(cnt))


def export_proteins(uri, dst, chunk_size=1000000):
    logging.info('starting')
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
        INNER JOIN UNIPARC.XREF UX ON IP.PROTEIN_AC = UX.AC AND UX.DELETED = 'N'
        INNER JOIN UNIPARC.PROTEIN UP ON UX.UPI = UP.UPI
        ORDER BY IP.PROTEIN_AC
        """
    )

    store = Store(dst)
    cnt = 0
    proteins = {}
    for acc, tax_id, name, dbcode, frag, length, seq_short, seq_long in cur:
        if seq_long is not None:
            seq = seq_long.read()
        else:
            seq = seq_short

        proteins[acc] = {
            'taxon': tax_id,
            'identifier': name,
            'isReviewed': dbcode == 'S',
            'isFrag': frag == 'Y',
            'length': length,
            'sequence': seq
        }

        if len(proteins) == chunk_size:
            cnt += chunk_size
            store.add(proteins)
            proteins = {}
            logging.info('{:>12}'.format(cnt))

    cur.close()
    con.close()

    cnt += len(proteins)
    store.add(proteins)
    logging.info('{:>12}'.format(cnt))

    store.close()


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
    member_dbs = {'B', 'D', 'F', 'H', 'I', 'J', 'M', 'N', 'P', 'Q', 'R', 'U', 'V', 'X', 'Y', 'g'}

    con, cur = dbms.connect(uri)

    # Using RN=2 to join with the second most recent action (the most recent is the current record)
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
                              ROW_NUMBER() OVER (PARTITION BY DBCODE ORDER BY TIMESTAMP DESC) RN
                            FROM INTERPRO.DB_VERSION_AUDIT
                            WHERE ACTION = 'U'
                          ) VA ON DB.DBCODE = VA.DBCODE AND VA.RN = 2
        """
    )

    databases = []
    for name, code, name_long, descr, version, rel_date, prev_version, prev_rel_date in cur:
        if code in member_dbs:
            db_type = 'entry'
        elif code in ('S', 'T', 'u'):
            if code == 'S':
                name = 'reviewed'
            elif code == 'T':
                name = 'unreviewed'

            db_type = 'protein'
        else:
            db_type = 'other'

        databases.append((
            name,
            name_long,
            descr,
            version,
            rel_date,
            db_type,
            prev_version,
            prev_rel_date
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


def get_entries(uri):
    con, cur = dbms.connect(uri)

    # InterPro entries (+ description)
    cur.execute(
        """
        SELECT LOWER(E.ENTRY_AC), LOWER(ET.ABBREV), E.NAME, E.SHORT_NAME, E.CREATED, CA.TEXT
        FROM INTERPRO.ENTRY E
        INNER JOIN INTERPRO.CV_ENTRY_TYPE ET ON E.ENTRY_TYPE = ET.CODE
        LEFT OUTER JOIN INTERPRO.ENTRY2COMMON E2C ON E.ENTRY_AC = E2C.ENTRY_AC
        LEFT OUTER JOIN INTERPRO.COMMON_ANNOTATION CA ON E2C.ANN_ID = CA.ANN_ID
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
                'descriptions': [],
                'integrated': None,
                'member_databases': {},
                'go_terms': [],
                'hierarchy': {},
                'citations': {},
                'cross_references': {}
            }

        if row[5] is not None:
            # todo: formatting descriptions
            '''
            Some annotations contain multiple blocks of text, 
            but some blocks might not be surrounded by <p> and </p> tags.
            
            Other blocks miss the opening <p> tag but have the closing </p> tag (or reverse).
            
            Shit's broken, yo.
            '''
            e['descriptions'].append(row[5])

    # InterPro entry contributing signatures
    cur.execute(
        """
        SELECT LOWER(EM.ENTRY_AC), LOWER(M.METHOD_AC), LOWER(DB.DBSHORT), M.NAME, M.DESCRIPTION
        FROM INTERPRO.ENTRY2METHOD EM
        INNER JOIN INTERPRO.METHOD M ON EM.METHOD_AC = M.METHOD_AC
        INNER JOIN INTERPRO.CV_DATABASE DB ON M.DBCODE = DB.DBCODE
        """
    )

    for row in cur:
        entry_ac = row[0]
        method_ac = row[1]
        method_db = row[2]
        method_name = row[3]
        method_descr = row[4]

        databases = entries[entry_ac]['member_databases']

        if method_db not in databases:
            databases[method_db] = {}

        databases[method_db][method_ac] = method_descr if method_descr else method_name

    # Remove InterPro entries without contributing signatures
    entries = {entry_ac: entry for entry_ac, entry in entries.items() if entry['member_databases']}

    # GO terms (InterPro entries)
    cur.execute(
        """
        SELECT LOWER(I2G.ENTRY_AC), GT.GO_ID, GT.NAME, GT.CATEGORY, GC.TERM_NAME
        FROM INTERPRO.INTERPRO2GO I2G
        INNER JOIN GO.TERMS@GOAPRO GT ON I2G.GO_ID = GT.GO_ID
        INNER JOIN GO.CV_CATEGORIES@GOAPRO GC ON GT.CATEGORY = GC.CODE
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
        entries[entry_ac]['hierarchy'] = _format_node(hierarchy, entries, accession)

    # Member database entries (with literature references, and integration)
    methods = {}
    cur.execute(
        """
        SELECT LOWER(M.METHOD_AC), M.NAME, LOWER(ET.ABBREV), M.DESCRIPTION, LOWER(DB.DBSHORT), M.ABSTRACT, M.ABSTRACT_LONG, M.METHOD_DATE, LOWER(E2M.ENTRY_AC)
        FROM INTERPRO.METHOD M
        INNER JOIN INTERPRO.CV_ENTRY_TYPE ET ON M.SIG_TYPE = ET.CODE
        INNER JOIN INTERPRO.CV_DATABASE DB ON M.DBCODE = DB.DBCODE
        LEFT OUTER JOIN INTERPRO.ENTRY2METHOD E2M ON M.METHOD_AC = E2M.METHOD_AC
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

        methods[method_ac] = {
            'accession': method_ac,
            'type': row[2],
            'name': row[3],
            'short_name': row[1],
            'database': row[4],
            'date': row[7],
            'descriptions': descr,
            'integrated': row[8],
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
        SELECT LOWER(E.ENTRY_AC), C.PUB_ID, C.PUBMED_ID, C.ISBN, C.VOLUME, C.ISSUE, C.YEAR, C.TITLE, C.URL, 
               C.RAWPAGES, C.MEDLINE_JOURNAL, C.ISO_JOURNAL, C.AUTHORS, C.DOI_URL
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

        cross_references = entries[entry_ac]['cross_references']

        if xref_db not in cross_references:
            cross_references[xref_db] = []

        cross_references[xref_db].append(xref_ac)

    cur.close()
    con.close()

    return list(entries.values())


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
            res = urllib.request.urlopen(base_url + title)
        except urllib.error.HTTPError as e:
            # Content can be retrieved with e.fp.read()
            continue
        except UnicodeEncodeError:
            print(acc, title)
        else:
            obj = json.loads(res.read().decode('utf-8'))
            thumbnail = obj.get('thumbnail')

            if thumbnail:
                filename, headers = urllib.request.urlretrieve(thumbnail['source'])

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


# def build_ida(upi, length, supermatches, min_overlap=20):
#     # Merge domain supermatches
#     sets = merge_supermatches(supermatches, min_overlap)
#
#     # Sort supermatches by position
#     supermatches = sorted([sm for sms in sets for sm in sms.supermatches], key=lambda e: (e.start, e.end))
#
#     if supermatches:
#         ida = '{}/{}#{}#'.format(length, upi, '~'.join([sm.format() for sm in supermatches]))
#         domain = '~'.join([sm.format_entries() for sm in supermatches])
#         # ida_key = int(hashlib.sha256(domain.encode('utf-8')).hexdigest(), 16) % 10 ** 8
#         ida_key = int(hashlib.md5(domain.encode('utf-8')).hexdigest(), 16)
#     else:
#         ida = ida_key = None
#
#     return ida, ida_key


def insert_supermatches(uri, proteins_f, prot_matches_f, **kwargs):
    chunk_size = kwargs.get('chunk_size', 100000)
    min_overlap = kwargs.get('min_overlap', 20)

    # Opening stores
    proteins = disk.Store(proteins_f)
    prot_matches = disk.Store(prot_matches_f)

    # Oracle connection
    con, cur = dbms.connect(uri)

    hierarchy = _get_relationships(cur)

    logging.info('dropping and creating table')
    try:
        cur.execute('DROP TABLE INTERPRO.SUPERMATCH2 CASCADE CONSTRAINTS')
    except :
        pass
    finally:
        cur.execute(
            """
            CREATE TABLE INTERPRO.SUPERMATCH2
            (
                PROTEIN_AC VARCHAR2(15) NOT NULL,
                ENTRY_AC VARCHAR2(9) NOT NULL,
                POS_FROM NUMBER(5) NOT NULL,
                POS_TO NUMBER(5) NOT NULL,
                DBCODE CHAR(1) NOT NULL 
            ) NOLOGGING
            """
        )

    cnt = 0
    data = []
    n_records = 0
    ts = time.time()
    logging.info('starting')
    for acc, protein in proteins.iter():
        matches = prot_matches.get(acc, [])

        # Merge matches into supermatches
        supermatches = []
        for m in matches:
            entry_ac = m['entry_ac']

            if entry_ac:
                entry_ac = entry_ac.upper()

                supermatches.append(
                    Supermatch(entry_ac, hierarchy.get_root(entry_ac), m['start'], m['end'])
                )

        # Merge overlapping supermatches
        sets = merge_supermatches(supermatches, min_overlap=min_overlap)
        for s in sets:
            for sm in s.supermatches:
                for entry_ac in sm.get_entries():
                    data.append((acc, entry_ac, sm.start, sm.end, 'S' if protein['isReviewed'] else 'T'))

                    if len(data) == chunk_size:
                        cur.executemany(
                            """
                            INSERT /*+APPEND*/ INTO INTERPRO.SUPERMATCH2 (
                                PROTEIN_AC, ENTRY_AC, POS_FROM, POS_TO, DBCODE
                            )
                            VALUES (:1, :2, :3, :4, :5)
                            """,
                            data
                        )
                        con.commit()
                        n_records += len(data)
                        data = []

        cnt += 1
        if not cnt % 1000000:
            logging.info('{:>12} ({:.0f} proteins/sec)'.format(cnt, cnt // (time.time() - ts)))

    if data:
        cur.executemany(
            """
            INSERT /*+APPEND*/ INTO INTERPRO.SUPERMATCH2 (
                PROTEIN_AC, ENTRY_AC, POS_FROM, POS_TO, DBCODE
            )
            VALUES (:1, :2, :3, :4, :5)
            """,
            data
        )
        con.commit()
        n_records += len(data)
        data = []

    logging.info('{:>12} ({:.0f} proteins/sec)'.format(cnt, cnt // (time.time() - ts)))
    logging.info('{} supermatches inserted'.format(n_records))

    logging.info('adding constraints')
    cur.execute(
        """
        ALTER TABLE INTERPRO.SUPERMATCH2
        ADD CONSTRAINT PK_SUPERMATCH2
        PRIMARY KEY (PROTEIN_AC, ENTRY_AC, POS_FROM, POS_TO, DBCODE)
        """
    )

    try:
        cur.execute(
            """
            ALTER TABLE INTERPRO.SUPERMATCH2
            ADD CONSTRAINT FK_SUPERMATCH2$PROTEIN_AC
            FOREIGN KEY (PROTEIN_AC) REFERENCES INTERPRO.PROTEIN (PROTEIN_AC)
            ON DELETE CASCADE
            """
        )
    except:
        logging.warning('could not create PROTEIN_AC constraint on INTERPRO.SUPERMATCH2')

    try:
        cur.execute(
            """
            ALTER TABLE INTERPRO.SUPERMATCH2
            ADD CONSTRAINT FK_SUPERMATCH2$ENTRY_AC
            FOREIGN KEY (ENTRY_AC) REFERENCES INTERPRO.ENTRY (ENTRY_AC)
            ON DELETE CASCADE
            """
        )
    except:
        logging.warning('could not create ENTRY_AC constraint on INTERPRO.SUPERMATCH2')

    logging.info('creating indexes')
    cur.execute('CREATE INDEX I_SUPERMATCH2$PROTEIN ON INTERPRO.SUPERMATCH2 (PROTEIN_AC) NOLOGGING')
    cur.execute('CREATE INDEX I_SUPERMATCH2$ENTRY ON INTERPRO.SUPERMATCH2 (ENTRY_AC) NOLOGGING')
    cur.execute('CREATE INDEX I_SUPERMATCH2$DBCODE$ENTRY ON INTERPRO.SUPERMATCH2 (DBCODE, ENTRY_AC) NOLOGGING')

    logging.info('gathering statistics')
    cur.execute(
        """
            BEGIN
                DBMS_STATS.GATHER_TABLE_STATS(:1, :2, cascade => TRUE);
            END;
        """,
        ('INTERPRO', 'SUPERMATCH2')
    )

    # Privileges
    cur.execute('GRANT SELECT ON INTERPRO.SUPERMATCH2 TO INTERPRO_SELECT')

    cur.close()
    con.close()

    logging.info('complete')


def intersect_matches(matches, sets, intersections):
    for acc1 in matches:
        if acc1 in sets:
            sets[acc1] += 1
        else:
            sets[acc1] = 1

        for acc2 in matches:
            if acc1 >= acc2:
                continue
            elif acc1 not in intersections:
                intersections[acc1] = {acc2: [0, 0]}
            elif acc2 not in intersections[acc1]:
                intersections[acc1][acc2] = [0, 0]

            m1 = matches[acc1][0]
            m2 = matches[acc2][0]
            o = min(m1[1], m2[1]) - max(m1[0], m2[0]) + 1

            l1 = m1[1] - m1[0] + 1
            l2 = m2[1] - m2[0] + 1

            if o > l1 * 0.5:
                # acc1 is in acc2 (because it overlaps acc2 at least 50%)
                intersections[acc1][acc2][0] += 1

            if o > l2 * 0.5:
                # acc2 is in acc1
                intersections[acc1][acc2][1] += 1


def jaccard(uri):
    con, cur = dbms.connect(uri)

    cur.execute(
        """
        SELECT PROTEIN_AC, ENTRY_AC, POS_FROM, POS_TO 
        FROM INTERPRO.SUPERMATCH2
        ORDER BY PROTEIN_AC, ENTRY_AC, POS_FROM, POS_TO        
        """
    )

    sets = {}
    intersections = {}
    matches = {}
    protein = None
    for protein_ac, entry_ac, start, end in cur:
        if protein_ac != protein:
            if protein:
                intersect_matches(matches, sets, intersections)

            matches = {}
            protein = protein_ac

        entry_ac = entry_ac.lower()

        # Consider only leftmost match
        if entry_ac not in matches:
            matches[entry_ac] = [(start, end)]

    cur.close()
    con.close()

    intersect_matches(matches, sets, intersections)

    entries = {}
    for acc1 in intersections:
        s1 = sets[acc1]

        for acc2 in intersections[acc1]:
            s2 = sets[acc2]
            i1, i2 = intersections[acc1][acc2]
            coef1 = i1 / (s1 + s2 - i1)
            coef2 = i1 / (s1 + s2 - i1)
            coef = (coef1 + coef2) * 0.5
            c1 = i1 / s1
            c2 = i2 / s2

            if coef >= 0.75 or c1 >= 0.51 or c2 >= 0.51:
                if acc1 in entries:
                    entries[acc1].append(acc2)
                else:
                    entries[acc1] = [acc2]

                if acc2 in entries:
                    entries[acc2].append(acc1)
                else:
                    entries[acc2] = [acc1]

    return entries
