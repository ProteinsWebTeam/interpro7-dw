# -*- coding: utf-8 -*-

from typing import List, Optional, Sequence

import cx_Oracle

from interpro7dw import logger
from interpro7dw.utils import Store, datadump


def export_comments(url: str, input: str, output: str,
                    dir: Optional[str]=None, processes: int=1):
    logger.info("starting")
    with Store(output, input, dir) as store:
        con = cx_Oracle.connect(url)
        cur = con.cursor()
        """
        Note on the TEXT structure: 
        Some comments have a title (e.g. Q01299) which is not retrieved 
        when joining on CC_STRUCTURE_TYPE_ID = 1
        """
        cur.execute(
            """
            SELECT 
              E.ACCESSION, B.ORDER_IN, B.TEXT, SS.TEXT
            FROM SPTR.DBENTRY@SWPREAD E
            INNER JOIN SPTR.COMMENT_BLOCK@SWPREAD B
              ON E.DBENTRY_ID = B.DBENTRY_ID
              AND B.COMMENT_TOPICS_ID = 2        -- FUNCTION comments
            LEFT OUTER JOIN SPTR.COMMENT_STRUCTURE@SWPREAD S
              ON B.COMMENT_BLOCK_ID = S.COMMENT_BLOCK_ID
              AND S.CC_STRUCTURE_TYPE_ID = 1      -- TEXT structure
            LEFT OUTER JOIN SPTR.COMMENT_SUBSTRUCTURE@SWPREAD SS
              ON S.COMMENT_STRUCTURE_ID = SS.COMMENT_STRUCTURE_ID
            WHERE E.ENTRY_TYPE IN (0, 1)          -- Swiss-Prot and TrEMBL
              AND E.MERGE_STATUS != 'R'           -- not 'Redundant'
              AND E.DELETED = 'N'                 -- not deleted
              AND E.FIRST_PUBLIC IS NOT NULL      -- published
            """
        )

        i = 0
        for accession, block_number, text, text2 in cur:
            store.append(accession, (block_number, text or text2))

            i += 1
            if not i % 1000000:
                store.sync()

                if not i % 10000000:
                    logger.info(f"{i:>12,}")

        cur.close()
        con.close()

        logger.info(f"{i:>12,}")
        size = store.merge(fn=_post_comments, processes=processes)
        logger.info(f"temporary files: {size/1024/1024:.0f} MB")


def _post_comments(blocks: Sequence[tuple]) -> List[str]:
    return [text for order, text in sorted(blocks)]


def export_descriptions(url: str, input: str, output: str,
                        dir: Optional[str]=None, processes: int=1):
    logger.info("starting")
    with Store(output, input, dir) as store:
        con = cx_Oracle.connect(url)
        cur = con.cursor()
        cur.execute(
            """
            SELECT ACCESSION, DESCR
            FROM (
                SELECT
                  E.ACCESSION, 
                  D.DESCR, 
                  ROW_NUMBER() OVER (
                    PARTITION BY E.ACCESSION 
                    ORDER BY CV.DESC_ID,    -- 1=RecName, 2=AltName, 3=SubName
                             CV.ORDER_ID,   -- Swiss-Prot manual order
                             D.DESCR        -- TrEMBL alphabetic order
                  ) RN
                FROM SPTR.DBENTRY@SWPREAD E
                INNER JOIN SPTR.DBENTRY_2_DESC@SWPREAD D
                  ON E.DBENTRY_ID = D.DBENTRY_ID
                  AND D.DESC_ID IN (1,4,11,13,16,23,25,28,35)  --Full description section
                INNER JOIN SPTR.CV_DESC@SWPREAD CV
                  ON D.DESC_ID = CV.DESC_ID
                WHERE E.ENTRY_TYPE IN (0, 1)
                  AND E.MERGE_STATUS != 'R'
                  AND E.DELETED = 'N'
                  AND E.FIRST_PUBLIC IS NOT NULL
            )
            WHERE RN = 1
            """
        )

        i = 0
        for accession, description in cur:
            store[accession] = description

            i += 1
            if not i % 1000000:
                store.sync()

                if not i % 10000000:
                    logger.info(f"{i:>12,}")

        cur.close()
        con.close()

        logger.info(f"{i:>12,}")
        size = store.merge(processes=processes)
        logger.info(f"temporary files: {size/1024/1024:.0f} MB")


def export_evidences(url: str, input: str, output: str,
                     dir: Optional[str]=None, processes: int=1):
    logger.info("starting")
    with Store(output, input, dir) as store:
        con = cx_Oracle.connect(url)
        cur = con.cursor()
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
        for accession, evidence, gene in cur:
            store[accession] = (evidence, gene)

            i += 1
            if not i % 1000000:
                store.sync()

                if not i % 10000000:
                    logger.info(f"{i:>12,}")

        cur.close()
        con.close()

        logger.info(f"{i:>12,}")
        size = store.merge(processes=processes)
        logger.info(f"temporary files: {size/1024/1024:.0f} MB")


def export_proteomes(url: str, input: str, output: str,
                     dir: Optional[str]=None, processes: int=1):
    logger.info("starting")
    with Store(output, input, dir) as store:
        con = cx_Oracle.connect(url)
        cur = con.cursor()
        """
        Without the DISTINCT, there would be duplicated rows, e.g.
        A0A059MHQ6  UP000024941
        A0A059MHQ6  UP000024941
        
        Even for duplicated rows, a given UniProt accession is associated
        to one unique UPID.
        
        It's just easier to remove the duplicates at the database level.
        """
        cur.execute(
            """
            SELECT DISTINCT E.ACCESSION, P.UPID
            FROM SPTR.DBENTRY@SWPREAD E
            INNER JOIN SPTR.PROTEOME2UNIPROT@SWPREAD P2U
              ON E.ACCESSION = P2U.ACCESSION AND E.TAX_ID = P2U.TAX_ID
            INNER JOIN SPTR.PROTEOME@SWPREAD P
              ON P2U.PROTEOME_ID = P.PROTEOME_ID
              AND P.IS_REFERENCE = 1
            WHERE E.ENTRY_TYPE IN (0, 1)
            AND E.MERGE_STATUS != 'R'
            AND E.DELETED = 'N'
            AND E.FIRST_PUBLIC IS NOT NULL
            """
        )

        i = 0
        for accession, upid in cur:
            store[accession] = upid

            i += 1
            if not i % 1000000:
                store.sync()

                if not i % 10000000:
                    logger.info(f"{i:>12,}")

        cur.close()
        con.close()

        logger.info(f"{i:>12,}")
        size = store.merge(processes=processes)
        logger.info(f"temporary files: {size/1024/1024:.0f} MB")


def export_ref_proteomes(url: str, output: str):
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    cur.execute(
        """
        SELECT P.UPID, P.PROTEOME_NAME, P.IS_REFERENCE, P.GC_SET_ACC, 
          TO_CHAR(P.PROTEOME_TAXID), SN.NAME
        FROM SPTR.PROTEOME@SWPREAD P
        LEFT OUTER JOIN TAXONOMY.SPTR_STRAIN@SWPREAD S
          ON P.PROTEOME_TAXID = S.TAX_ID
        LEFT OUTER JOIN TAXONOMY.SPTR_STRAIN_NAME@SWPREAD SN
          ON S.STRAIN_ID = SN.STRAIN_ID
        WHERE P.IS_REFERENCE = 1
        """
    )

    proteomes = {}
    for row in cur:
        upid = row[0]

        if upid in proteomes:
            continue

        proteomes[upid] = {
            "name": row[1],
            "is_reference": row[2] != 0,
            "assembly": row[3],
            "taxon_id": row[4],
            "strain": row[5]
        }

    cur.close()
    con.close()

    datadump(output, proteomes)
