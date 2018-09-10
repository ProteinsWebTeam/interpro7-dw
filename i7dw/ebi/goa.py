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

    """
    If `EE.GO_EVIDENCE` == 'IEA', it's an electronically assigned GO term 
        see:
            * https://www.uniprot.org/help/gene_ontology
            * http://geneontology.org/page/guide-go-evidence-codes#IEA 
    """
    cur.execute(
        """
        SELECT A.ENTITY_ID, T.GO_ID, T.NAME, D.DEFINITION, T.CATEGORY, C.TERM_NAME, A.SOURCE
        FROM GO.ANNOTATIONS@GOAPRO A
        INNER JOIN GO.ECO2EVIDENCE@GOAPRO EE ON A.ECO_ID = EE.ECO_ID
        INNER JOIN GO.CV_SOURCES@GOAPRO S ON A.SOURCE = S.CODE
        INNER JOIN GO.TERMS@GOAPRO T ON A.GO_ID = T.GO_ID
        INNER JOIN GO.CV_CATEGORIES@GOAPRO C ON T.CATEGORY = C.CODE
        INNER JOIN GO.DEFINITIONS@GOAPRO D ON T.GO_ID = D.GO_ID
        WHERE A.ENTITY_TYPE = 'protein'
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
        from_interpro = row[6] == 'IPRO'

        if acc in proteins:
            terms = proteins[acc]
        else:
            if len(proteins) == chunk_size:
                for _acc in proteins:
                    terms = proteins[_acc]
                    proteins[_acc] = sorted(terms.values(), key=lambda x: x['identifier'])

                cnt += len(proteins)
                store.add(proteins)
                proteins = {}

                if not cnt % 1000000:
                    logging.info('{:>12}'.format(cnt))

            terms = proteins[acc] = {}

        if term_id not in terms:
            terms[term_id] = {
                'identifier': term_id,
                'name': term_name,
                # 'definition': term_def,
                'category': {
                    'code': cat_code,
                    'name': cat_name
                },
                'interpro': from_interpro
            }
        elif from_interpro:
            terms[term_id]['interpro'] = from_interpro

    cur.close()
    con.close()

    for _acc in proteins:
        terms = proteins[_acc]
        proteins[_acc] = sorted(terms.values(), key=lambda x: x['identifier'])

    cnt += len(proteins)
    store.add(proteins)
    store.close()

    logging.info('{:>12}'.format(cnt))
