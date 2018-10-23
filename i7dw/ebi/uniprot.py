#!/usr/bin/env python
# -*- coding: utf-8 -*-

import json
import logging

from i7dw import dbms
from i7dw.disk import Store

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s: %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)


def export_protein2comments(uri, src, dst, tmpdir=None, flush=1000000):
    logging.info("starting")

    with open(src, "rt") as fh:
        keys = json.load(fh)

    s = Store(dst, keys, tmpdir)
    con, cur = dbms.connect(uri)

    # Topic #2 is "FUNCTION"
    cur.execute(
        """
        SELECT E.ACCESSION, CSS.TEXT
        FROM SPTR.DBENTRY@SWPREAD E
        INNER JOIN SPTR.COMMENT_BLOCK@SWPREAD CB 
          ON E.DBENTRY_ID = CB.DBENTRY_ID
        INNER JOIN SPTR.COMMENT_STRUCTURE@SWPREAD CS 
          ON CB.COMMENT_BLOCK_ID = CS.COMMENT_BLOCK_ID
        INNER JOIN SPTR.COMMENT_SUBSTRUCTURE@SWPREAD CSS 
          ON CS.COMMENT_STRUCTURE_ID = CSS.COMMENT_STRUCTURE_ID
        WHERE E.ENTRY_TYPE IN (0, 1)
        AND E.MERGE_STATUS != 'R'
        AND E.DELETED = 'N'
        AND E.FIRST_PUBLIC IS NOT NULL
        AND CB.COMMENT_TOPICS_ID = 2
        """
    )

    i = 0
    for acc, text in cur:
        s.append(acc, text)

        i += 1
        if not i % flush:
            s.flush()

        if not i % 1000000:
            logging.info("{:>12}".format(i))

    cur.close()
    con.close()
    logging.info("{:>12}".format(i))
    size = s.merge()
    logging.info("temporary files: {} bytes".format(size))


def parse_descriptions(item: list) -> tuple:
    """
    item is a list of UniProt descriptions as tuples
    (name, category, subcategory, order)
    """
    rec_name = {
        'fullName': None,
        'shortNames': []
    }
    alt_names = []
    sub_names = []

    for name, catg, subcatg, order in sorted(item, key=lambda x: x[3]):
        if catg == 'RecName':
            if subcatg == 'Full':
                rec_name['fullName'] = name
            else:
                rec_name['shortNames'].append(name)
        elif catg == 'AltName':
            if subcatg == 'Full':
                alt_names.append({
                    'fullName': name,
                    'shortNames': []
                })
            else:
                try:
                    alt_name = alt_names[-1]
                except IndexError:
                    pass
                else:
                    alt_name['shortNames'].append(name)
        elif catg == 'SubName' and subcatg == 'Full':
            sub_names.append({'fullName': name})

    other_names = {}
    if rec_name['fullName']:
        other_names['recommendedName'] = rec_name

    if alt_names:
        other_names['alternativeNames'] = alt_names

    if sub_names:
        other_names['submittedNames'] = sub_names

    if other_names.get('recommendedName'):
        name = other_names['recommendedName']['fullName']
    elif other_names.get('submittedNames'):
        name = other_names['submittedNames'][0]['fullName']
    else:
        name = None

    return name, other_names


def export_protein2names(uri, src, dst, tmpdir=None, flush=1000000):
    logging.info("starting")

    with open(src, "rt") as fh:
        keys = json.load(fh)

    s = Store(dst, keys, tmpdir)
    con, cur = dbms.connect(uri)

    # Topic #2 is "FUNCTION"
    cur.execute(
        """
        SELECT 
          E.ACCESSION, E2D.DESCR, CV.CATG_TYPE, CV.SUBCATG_TYPE, CV.ORDER_IN
        FROM SPTR.DBENTRY@SWPREAD E
        INNER JOIN SPTR.DBENTRY_2_DESC@SWPREAD E2D 
          ON E.DBENTRY_ID = E2D.DBENTRY_ID
        INNER JOIN SPTR.CV_DESC@SWPREAD CV 
          ON E2D.DESC_ID = CV.DESC_ID
        WHERE E.ENTRY_TYPE IN (0, 1)
        AND E.MERGE_STATUS != 'R'
        AND E.DELETED = 'N'
        AND E.FIRST_PUBLIC IS NOT NULL
        AND CV.SECTION_TYPE = 'Main'
        AND CV.SUBCATG_TYPE IN ('Full', 'Short')
        """
    )

    i = 0
    for row in cur:
        s.append(row[0], (row[1], row[2], row[3], row[4]))

        i += 1
        if not i % flush:
            s.flush()

        if not i % 1000000:
            logging.info("{:>12}".format(i))

    cur.close()
    con.close()
    logging.info("{:>12}".format(i))
    size = s.merge(func=parse_descriptions)
    logging.info("temporary files: {} bytes".format(size))


def export_protein2supplementary(uri, src, dst, tmpdir=None, flush=1000000):
    logging.info("starting")

    with open(src, "rt") as fh:
        keys = json.load(fh)

    s = Store(dst, keys, tmpdir)
    con, cur = dbms.connect(uri)
    cur.execute(
        """
        SELECT ACCESSION, PROTEIN_EXISTENCE_ID, NAME
        FROM (
          SELECT
            E.ACCESSION,
            E.PROTEIN_EXISTENCE_ID,
            GN.NAME,
            ROW_NUMBER() OVER (
              PARTITION BY E.ACCESSION 
              ORDER BY GN.GENE_NAME_TYPE_ID
            ) RN
          FROM SPTR.DBENTRY@SWPREAD E
          LEFT OUTER JOIN SPTR.GENE@SWPREAD G 
            ON E.DBENTRY_ID = G.DBENTRY_ID
          LEFT OUTER JOIN SPTR.GENE_NAME@SWPREAD GN 
            ON G.GENE_ID = GN.GENE_ID
          WHERE E.ENTRY_TYPE IN (0, 1)
          AND E.MERGE_STATUS != 'R'
          AND E.DELETED = 'N'
          AND E.FIRST_PUBLIC IS NOT NULL
        )
        WHERE RN = 1
        """
    )

    i = 0
    for acc, evi, gene in cur:
        s[acc] = (evi, gene)

        i += 1
        if not i % flush:
            s.flush()

        if not i % 1000000:
            logging.info("{:>12}".format(i))

    cur.close()
    con.close()
    logging.info("{:>12}".format(i))
    size = s.merge()
    logging.info("temporary files: {} bytes".format(size))


def export_protein2proteome(uri, src, dst, tmpdir=None, flush=1000000):
    logging.info("starting")

    with open(src, "rt") as fh:
        keys = json.load(fh)

    s = Store(dst, keys, tmpdir)
    con, cur = dbms.connect(uri)
    cur.execute(
        """
        SELECT
          DISTINCT E.ACCESSION, LOWER(P.UPID)
        FROM SPTR.DBENTRY@SWPREAD E
        INNER JOIN SPTR.PROTEOME2UNIPROT@SWPREAD P2U 
          ON E.ACCESSION = P2U.ACCESSION
        INNER JOIN SPTR.PROTEOME@SWPREAD P 
          ON P2U.PROTEOME_ID = P.PROTEOME_ID
        WHERE E.ENTRY_TYPE IN (0, 1)
        AND E.MERGE_STATUS != 'R'
        AND E.DELETED = 'N'
        AND E.FIRST_PUBLIC IS NOT NULL
        """
    )

    i = 0
    for acc, upid in cur:
        s[acc] = upid

        i += 1
        if not i % flush:
            s.flush()

        if not i % 1000000:
            logging.info("{:>12}".format(i))

    cur.close()
    con.close()
    logging.info("{:>12}".format(i))
    size = s.merge()
    logging.info("temporary files: {} bytes".format(size))


def get_proteomes(uri: str) -> dict:
    con, cur = dbms.connect(uri)

    cur.execute(
        """
        SELECT
          LOWER(P.UPID),
          P.PROTEOME_NAME,
          P.IS_REFERENCE,
          P.GC_SET_ACC,
          TO_CHAR(P.PROTEOME_TAXID),
          SN.NAME
        FROM SPTR.PROTEOME@SWPREAD P
        LEFT OUTER JOIN TAXONOMY.SPTR_STRAIN@SWPREAD S ON P.PROTEOME_TAXID = S.TAX_ID
        LEFT OUTER JOIN TAXONOMY.SPTR_STRAIN_NAME@SWPREAD SN ON S.STRAIN_ID = SN.STRAIN_ID
        ORDER BY P.UPID
        """
    )

    proteomes = {}
    for row in cur:
        upid = row[0]

        if upid not in proteomes:
            proteomes[upid] = {
                "accession": upid,
                "name": row[1],
                "is_reference": bool(row[2]),
                "assembly": row[3],
                "tax_id": row[4],
                "strain": row[5]
            }

    cur.close()
    con.close()

    return proteomes
