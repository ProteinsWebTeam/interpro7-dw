#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging

from i7dw import dbms
from i7dw.disk import Store

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s: %(levelname)s: %(message)s',
    datefmt='%y-%m-%d %H:%M:%S'
)


def export_annotations(uri, dst, chunk_size=1000000):
    logging.info('starting')
    con, cur = dbms.connect(uri)
    cur.execute(
        """
        SELECT A.ENTITY_ID, T.GO_ID, T.NAME, D.DEFINITION, T.CATEGORY, C.TERM_NAME
        FROM GO.ANNOTATIONS@GOAPRO A
        INNER JOIN GO.ECO2EVIDENCE@GOAPRO EE ON A.ECO_ID = EE.ECO_ID
        INNER JOIN GO.CV_SOURCES@GOAPRO S ON A.SOURCE = S.CODE
        INNER JOIN GO.TERMS@GOAPRO T ON A.GO_ID = T.GO_ID
        INNER JOIN GO.CV_CATEGORIES@GOAPRO C ON T.CATEGORY = C.CODE
        INNER JOIN GO.DEFINITIONS@GOAPRO D ON T.GO_ID = D.GO_ID
        WHERE A.ENTITY_TYPE = 'protein'
        AND EE.GO_EVIDENCE != 'IEA'
        AND S.IS_PUBLIC = 'Y'
        ORDER BY A.ENTITY_ID
        """
    )

    store = Store(dst)
    proteins = {}
    cnt = 0
    for row in cur:
        acc = row[0]
        term_id = row[1]
        term_name = row[2]
        term_def = row[3]
        cat_code = row[4]
        cat_name = row[5]

        if acc in proteins:
            p = proteins[acc]
        else:
            if len(proteins) == chunk_size:
                cnt += len(proteins)
                store.add(proteins)
                proteins = {}
                logging.info('{:>12}'.format(cnt))

            p = proteins[acc] = []

        p.append({
            'identifier': term_id,
            'name': term_name,
            'definition': term_def,
            'category': {
                'code': cat_code,
                'name': cat_name
            }
        })

    cur.close()
    con.close()

    cnt += len(proteins)
    store.add(proteins)
    store.close()

    logging.info('{:>12}'.format(cnt))
