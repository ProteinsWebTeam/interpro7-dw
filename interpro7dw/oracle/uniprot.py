import cx_Oracle

from interpro7dw import logger
from interpro7dw.utils import NewStore


def export_functions(url: str, output: str):
    logger.info("starting")

    with NewStore(output, mode="w") as store:
        con = cx_Oracle.connect(url)
        cur = con.cursor()
        """
        Note on the TEXT structure: 
        Some comments have a title (e.g. Q01299) which is not retrieved 
        when joining on CC_STRUCTURE_TYPE_ID = 1
        """
        cur.execute(
            """
            SELECT E.ACCESSION, B.ORDER_IN, NVL(B.TEXT, SS.TEXT)
            FROM SPTR.DBENTRY E
            INNER JOIN SPTR.COMMENT_BLOCK B
              ON E.DBENTRY_ID = B.DBENTRY_ID
              AND B.COMMENT_TOPICS_ID = 2        -- FUNCTION comments
            LEFT OUTER JOIN SPTR.COMMENT_STRUCTURE S
              ON B.COMMENT_BLOCK_ID = S.COMMENT_BLOCK_ID
              AND S.CC_STRUCTURE_TYPE_ID = 1      -- TEXT structure
            LEFT OUTER JOIN SPTR.COMMENT_SUBSTRUCTURE SS
              ON S.COMMENT_STRUCTURE_ID = SS.COMMENT_STRUCTURE_ID
            WHERE E.ENTRY_TYPE IN (0, 1)          -- Swiss-Prot and TrEMBL
              AND E.MERGE_STATUS != 'R'           -- not 'Redundant'
              AND E.DELETED = 'N'                 -- not deleted
              AND E.FIRST_PUBLIC IS NOT NULL      -- published
            """
        )

        for i, (accession, block_number, text) in enumerate(cur):
            store.add(accession, (block_number, text))

            if (i + 1) % 10000000 == 0:
                logger.info(f"{i + 1:>15,}")

        cur.close()
        con.close()

        logger.info(f"{i + 1:>15,}")
        size = store.digest(apply=_sort_functions)
        logger.info(f"temporary files: {size / 1024 / 1024:.0f} MB")


def _sort_functions(values):
    return [text for block_number, text in sorted(values)]


def export_names(url: str, output: str):
    logger.info("starting")

    with NewStore(output, mode="w") as store:
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
                             CV.ORDER_IN,   -- Swiss-Prot manual order
                             D.DESCR        -- TrEMBL alphabetic order
                  ) RN
                FROM SPTR.DBENTRY E
                INNER JOIN SPTR.DBENTRY_2_DESC D
                  ON E.DBENTRY_ID = D.DBENTRY_ID
                  AND D.DESC_ID IN (1,4,11,13,16,23,25,28,35)  --Full description section
                INNER JOIN SPTR.CV_DESC CV
                  ON D.DESC_ID = CV.DESC_ID
                WHERE E.ENTRY_TYPE IN (0, 1)
                  AND E.MERGE_STATUS != 'R'
                  AND E.DELETED = 'N'
                  AND E.FIRST_PUBLIC IS NOT NULL
            )
            WHERE RN = 1
            """
        )

        for i, (accession, name) in enumerate(cur):
            store.add(accession, name)

            if (i + 1) % 10000000 == 0:
                logger.info(f"{i+1:>15,}")

        cur.close()
        con.close()

        logger.info(f"{i + 1:>15,}")
        size = store.digest(apply=lambda values: values[0])
        logger.info(f"temporary files: {size / 1024 / 1024:.0f} MB")


def export_proteomes(url: str, output: str):
    logger.info("starting")

    with NewStore(output, mode="w") as store:
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
            FROM SPTR.DBENTRY E
            INNER JOIN SPTR.PROTEOME2UNIPROT P2U
              ON E.ACCESSION = P2U.ACCESSION AND E.TAX_ID = P2U.TAX_ID
            INNER JOIN SPTR.PROTEOME P
              ON P2U.PROTEOME_ID = P.PROTEOME_ID
              AND P.IS_REFERENCE = 1
            WHERE E.ENTRY_TYPE IN (0, 1)
            AND E.MERGE_STATUS != 'R'
            AND E.DELETED = 'N'
            AND E.FIRST_PUBLIC IS NOT NULL
            """
        )

        for i, (accession, upid) in enumerate(cur):
            store.add(accession, upid)

            if (i + 1) % 10000000 == 0:
                logger.info(f"{i + 1:>15,}")

        cur.close()
        con.close()

        logger.info(f"{i + 1:>15,}")
        size = store.digest(apply=lambda values: values[0])
        logger.info(f"temporary files: {size / 1024 / 1024:.0f} MB")