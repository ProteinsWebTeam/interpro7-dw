import re

import oracledb

from interpro7dw.utils import logger
from interpro7dw.utils.store import KVStore, KVStoreBuilder


def export_entry2functions(uri: str, proteins_file: str, output: str,
                           tempdir: str | None = None):
    logger.info("starting")
    with KVStore(proteins_file) as store:
        keys = store.get_keys()

    with KVStoreBuilder(output, keys=keys, tempdir=tempdir) as store:
        con = oracledb.connect(uri)
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

            if (i + 1) % 1e7 == 0:
                logger.info(f"{i + 1:>15,}")

        cur.close()
        con.close()

        logger.info(f"{i + 1:>15,}")

        store.build(apply=_sort_blocks)
        logger.info(f"temporary files: {store.get_size() / 1024 ** 2:.0f} MB")

    logger.info("done")


def _sort_blocks(blocks: list[tuple[int, str]]) -> list[str]:
    return [text for order, text in sorted(blocks)]


def export_entry2name(uri: str, proteins_file: str, output: str,
                      tempdir: str | None = None):
    logger.info("starting")
    with KVStore(proteins_file) as store:
        keys = store.get_keys()

    with KVStoreBuilder(output, keys=keys, tempdir=tempdir) as store:
        con = oracledb.connect(uri)
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
                    ORDER BY CV.DESC_ID,    -- 1=RecName,2=AltName,3=SubName
                             CV.ORDER_IN,   -- Swiss-Prot manual order
                             D.DESCR        -- TrEMBL alphabetic order
                  ) RN
                FROM SPTR.DBENTRY E
                INNER JOIN SPTR.DBENTRY_2_DESC D
                  ON E.DBENTRY_ID = D.DBENTRY_ID
                                    --Full description section
                  AND D.DESC_ID IN (1,4,11,13,16,23,25,28,35)
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

            if (i + 1) % 1e7 == 0:
                logger.info(f"{i + 1:>15,}")

        cur.close()
        con.close()

        logger.info(f"{i + 1:>15,}")

        store.build(apply=store.get_first)
        logger.info(f"temporary files: {store.get_size() / 1024 ** 2:.0f} MB")

    logger.info("done")


def export_entry2evidence(uri: str, proteins_file: str, output: str,
                          tempdir: str | None = None):
    logger.info("starting")
    with KVStore(proteins_file) as store:
        keys = store.get_keys()

    with KVStoreBuilder(output, keys=keys, tempdir=tempdir) as store:
        con = oracledb.connect(uri)
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
              FROM SPTR.DBENTRY E
              LEFT OUTER JOIN SPTR.GENE G
                ON E.DBENTRY_ID = G.DBENTRY_ID
              LEFT OUTER JOIN SPTR.GENE_NAME GN
                ON G.GENE_ID = GN.GENE_ID
              WHERE E.ENTRY_TYPE IN (0, 1)
              AND E.MERGE_STATUS != 'R'
              AND E.DELETED = 'N'
              AND E.FIRST_PUBLIC IS NOT NULL
            )
            WHERE RN = 1
            """
        )

        for i, (accession, evidence, gene) in enumerate(cur):
            store.add(accession, (evidence, gene))

            if (i + 1) % 1e7 == 0:
                logger.info(f"{i + 1:>15,}")

        cur.close()
        con.close()

        logger.info(f"{i + 1:>15,}")

        store.build(apply=store.get_first)
        logger.info(f"temporary files: {store.get_size() / 1024 ** 2:.0f} MB")

    logger.info("done")


def export_entry2proteome(uri: str, proteins_file: str, output: str,
                          tempdir: str | None = None):
    logger.info("starting")
    with KVStore(proteins_file) as store:
        keys = store.get_keys()

    with KVStoreBuilder(output, keys=keys, tempdir=tempdir) as store:
        con = oracledb.connect(uri)
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

            if (i + 1) % 1e7 == 0:
                logger.info(f"{i + 1:>15,}")

        cur.close()
        con.close()

        logger.info(f"{i + 1:>15,}")

        store.build(apply=store.get_first)
        logger.info(f"temporary files: {store.get_size() / 1024 ** 2:.0f} MB")

    logger.info("done")


def get_swissprot2enzyme(url: str) -> dict[str, list[str]]:
    """Returns EC numbers and the Swiss-Prot proteins they are associated with.

    :param url: UniProt Oracle connection string.
    :return: A dictionary where the key is an UniProt accession
        and the value a list of EC numbers.
    """
    con = oracledb.connect(url)
    cur = con.cursor()
    cur.execute(
        """
        SELECT DISTINCT E.ACCESSION, D.DESCR
        FROM SPTR.DBENTRY E
        INNER JOIN SPTR.DBENTRY_2_DESC D
          ON E.DBENTRY_ID = D.DBENTRY_ID
        INNER JOIN SPTR.CV_DESC C
          ON D.DESC_ID = C.DESC_ID
        WHERE E.ENTRY_TYPE = 0            -- Swiss-Prot
          AND E.MERGE_STATUS != 'R'       -- not 'Redundant'
          AND E.DELETED = 'N'             -- not deleted
          AND E.FIRST_PUBLIC IS NOT NULL  -- published
          AND C.SUBCATG_TYPE = 'EC'
        """
    )

    proteins = {}
    for acc, ecno in cur:
        # Accepts X.X.X.X or X.X.X.-
        # Does not accept preliminary EC numbers (e.g. X.X.X.nX)
        if re.match(r"(\d+\.){3}(\d+|-)$", ecno):
            try:
                proteins[acc].append(ecno)
            except KeyError:
                proteins[acc] = [ecno]

    cur.execute(
      """
      SELECT DISTINCT E.ACCESSION
      FROM SPTR.DBENTRY E
      JOIN SPTR.DBENTRY_2_DESC D ON E.DBENTRY_ID = D.DBENTRY_ID
      JOIN SPTR.CV_DESC C ON D.DESC_ID = C.DESC_ID
      LEFT JOIN (
          SELECT DISTINCT D2.DBENTRY_ID
          FROM SPTR.DBENTRY_2_DESC D2
          JOIN SPTR.CV_DESC C2 ON D2.DESC_ID = C2.DESC_ID
          WHERE C2.SUBCATG_TYPE = 'EC'
      ) EC_ANNOTATED ON E.DBENTRY_ID = EC_ANNOTATED.DBENTRY_ID
      WHERE E.ENTRY_TYPE = 0
        AND E.MERGE_STATUS != 'R'
        AND E.DELETED = 'N'
        AND E.FIRST_PUBLIC IS NOT NULL
        AND EC_ANNOTATED.DBENTRY_ID IS NULL
      """
    )
    no_ec_proteins = set()
    for row in cur:
      no_ec_proteins.add(row[0])

    cur.close()
    con.close()

    return proteins, no_ec_proteins


def get_swissprot2reactome(url: str) -> dict[str, list[tuple]]:
    """Returns the Reactome pathways and the Swiss-Prot proteins
    they are associated with.

    :param url: UniProt Oracle connection string.
    :return: A dictionary where the key is an UniProt accession
        and the value a list of pathways represented as tuples (ID, name).
    """
    con = oracledb.connect(url)
    cur = con.cursor()
    cur.execute(
        """
        SELECT DISTINCT E.ACCESSION, D.PRIMARY_ID, D.SECONDARY_ID
        FROM SPTR.DBENTRY E
        INNER JOIN SPTR.DBENTRY_2_DATABASE D
          ON E.DBENTRY_ID = D.DBENTRY_ID 
          AND D.DATABASE_ID = 'GK'          -- Reactome
        WHERE E.ENTRY_TYPE = 0              -- Swiss-Prot
          AND E.MERGE_STATUS != 'R'         -- not 'Redundant'
          AND E.DELETED = 'N'               -- not deleted
          AND E.FIRST_PUBLIC IS NOT NULL    -- published
        """
    )

    proteins = {}
    for uniprot_acc, pathway_id, pathway_name in cur:
        try:
            proteins[uniprot_acc].append((pathway_id, pathway_name))
        except KeyError:
            proteins[uniprot_acc] = [(pathway_id, pathway_name)]

    cur.close()
    con.close()

    return proteins
