# -*- coding: utf-8 -*-

import json
from typing import Optional

import cx_Oracle

from i7dw import logger
from i7dw.io import Store


def export_ppi(url: str, proteins_file: str, dst: str, processes: int=1,
               tmpdir: Optional[str]=None, sync_frequency: int=1000000):

    with open(proteins_file, "rt") as fh:
        keys = json.load(fh)

    with Store(dst, keys, tmpdir) as store:
        con = cx_Oracle.connect(url)
        cur = con.cursor()
        cur.execute(
            """
            SELECT 
              C1.INTERACTION_AC, I1.AC, IX1.PRIMARYID, CVT1.SHORTLABEL, 
              I2.AC, IX2.PRIMARYID, CVT2.SHORTLABEL, PX.PRIMARYID
            FROM IA_INTERACTOR_XREF IX1
            INNER JOIN IA_INTERACTOR I1 
              ON IX1.PARENT_AC = I1.AC
            INNER JOIN IA_COMPONENT C1 
              ON I1.AC = C1.INTERACTOR_AC
            INNER JOIN IA_COMPONENT C2 
              ON C1.INTERACTION_AC = C2.INTERACTION_AC 
              AND C1.INTERACTOR_AC != C2.INTERACTOR_AC
            INNER JOIN IA_INTERACTOR I2 
              ON C2.INTERACTOR_AC = I2.AC
            INNER JOIN IA_INTERACTOR_XREF IX2 
              ON I2.AC = IX2.PARENT_AC
            INNER JOIN IA_CONTROLLEDVOCAB CVT1 
              ON I1.INTERACTORTYPE_AC = CVT1.AC
            INNER JOIN IA_CONTROLLEDVOCAB CVT2 
              ON I2.INTERACTORTYPE_AC = CVT2.AC
            INNER JOIN IA_INT2EXP I2E 
              ON C1.INTERACTION_AC = I2E.INTERACTION_AC
            INNER JOIN IA_EXPERIMENT E 
              ON I2E.EXPERIMENT_AC = E.AC
            INNER JOIN IA_PUBLICATION P 
              ON E.PUBLICATION_AC = P.AC
            INNER JOIN IA_PUBLICATION_XREF PX 
              ON P.AC = PX.PARENT_AC
            WHERE IX1.DATABASE_AC = 'EBI-31'      -- uniprotkb
            AND IX1.QUALIFIER_AC = 'EBI-28'       -- identity
            AND IX2.DATABASE_AC = 'EBI-31'        -- uniprotkb
            AND IX2.QUALIFIER_AC = 'EBI-28'       -- identity
            AND PX.QUALIFIER_AC = 'EBI-49940'     -- primary-reference
            AND PX.DATABASE_AC = 'EBI-34'         -- pubmed
            """
        )

        i = 0
        for row in cur:
            intact_id = row[0]
            # interactor_1_id = row[1]
            interactor_1_ac = row[2]
            interactor_1_type = row[3]
            # interactor_2_id = row[4]
            interactor_2_ac = row[5]
            interactor_2_type = row[6]

            try:
                pubmed_id = int(row[7])
            except (ValueError, TypeError):
                continue

            store.update(interactor_1_ac, {
                "type": interactor_1_type,
                "interactions": {
                    interactor_2_ac: {
                        "type": interactor_2_type,
                        "pubmed_id": pubmed_id
                    }
                }
            })

            i += 1
            if sync_frequency and not i % sync_frequency:
                store.sync()

            if not i % 10000000:
                logger.info(f"{i:>15,}")

        cur.close()
        con.close()

        logger.info(f"{i:>15,}")

        size = store.merge(processes=processes)
        logger.info(f"temporary files: {size/1024/1024:.0f} MB")
