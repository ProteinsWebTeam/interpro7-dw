# -*- coding: utf-8 -*-

import json
from typing import Optional

import cx_Oracle

from i7dw import io, logger


def export_comments(url: str, src: str, dst: str,
                    tmpdir: Optional[str]=None, processes: int=1,
                    sync_frequency: int=1000000):
    logger.info("starting")

    with open(src, "rt") as fh:
        keys = json.load(fh)

    with io.Store(dst, keys, tmpdir) as store:
        con = cx_Oracle.connect(url)
        cur = con.cursor()

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
            store.append(acc, text)

            i += 1
            if sync_frequency and not i % sync_frequency:
                store.sync()

            if not i % 10000000:
                logger.info(f"{i:>12,}")

        cur.close()
        con.close()

        logger.info(f"{i:>12,}")
        size = store.merge(processes=processes)
        logger.info(f"temporary files: {size/1024/1024:.0f} MB")


def select_name(item: list) -> str:
    """
    item is a list of UniProt descriptions as tuples
    (name, order)
    """
    item.sort(key=lambda x: x[1])
    return item[0]


def export_descriptions(url: str, src: str, dst: str, processes: int=1,
                        tmpdir: Optional[str]=None,
                        sync_frequency: int=1000000):
    logger.info("starting")

    with open(src, "rt") as fh:
        keys = json.load(fh)

    with io.Store(dst, keys, tmpdir) as store:
        con = cx_Oracle.connect(url)
        cur = con.cursor()
        cur.execute(
            """
            SELECT E.ACCESSION, E2D.DESCR, CV.ORDER_IN
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
            store.append(row[0], (row[1], row[2]))

            i += 1
            if sync_frequency and not i % sync_frequency:
                store.sync()

            if not i % 10000000:
                logger.info(f"{i:>12,}")

        cur.close()
        con.close()

        logger.info(f"{i:>12,}")
        size = store.merge(func=select_name, processes=processes)
        logger.info(f"temporary files: {size/1024/1024:.0f} MB")


def export_misc(url: str, src: str, dst: str, processes: int=1,
                tmpdir: Optional[str]=None, sync_frequency: int=1000000):
    logger.info("starting")

    with open(src, "rt") as fh:
        keys = json.load(fh)

    with io.Store(dst, keys, tmpdir) as store:
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
        for acc, evi, gene in cur:
            store[acc] = (evi, gene)

            i += 1
            if sync_frequency and not i % sync_frequency:
                store.sync()

            if not i % 10000000:
                logger.info(f"{i:>12,}")

        cur.close()
        con.close()

        logger.info(f"{i:>12,}")
        size = store.merge(processes=processes)
        logger.info(f"temporary files: {size/1024/1024:.0f} MB")


def export_proteomes(url: str, src: str, dst: str, processes: int=1,
                     tmpdir: Optional[str]=None, sync_frequency: int=1000000):
    logger.info("starting")

    with open(src, "rt") as fh:
        keys = json.load(fh)

    with io.Store(dst, keys, tmpdir) as store:
        con = cx_Oracle.connect(url)
        cur = con.cursor()

        """
        DISTINCT not really necessary (but there *are* duplicated rows),
        but by using it the Store will hold less temporary values
        """
        cur.execute(
            """
            SELECT DISTINCT E.ACCESSION, P.UPID
            FROM SPTR.DBENTRY@SWPREAD E
            INNER JOIN SPTR.PROTEOME2UNIPROT@SWPREAD P2U
              ON E.ACCESSION = P2U.ACCESSION AND E.TAX_ID = P2U.TAX_ID
            INNER JOIN SPTR.PROTEOME@SWPREAD P
              ON P2U.PROTEOME_ID = P.PROTEOME_ID
            WHERE E.ENTRY_TYPE IN (0, 1)
            AND E.MERGE_STATUS != 'R'
            AND E.DELETED = 'N'
            AND E.FIRST_PUBLIC IS NOT NULL
            AND P.IS_REFERENCE = 1
            """
        )

        i = 0
        for acc, upid in cur:
            store[acc] = upid

            i += 1
            if sync_frequency and not i % sync_frequency:
                store.sync()

            if not i % 10000000:
                logger.info(f"{i:>12,}")

        cur.close()
        con.close()
        logger.info(f"{i:>12,}")

        size = store.merge(processes=processes)
        logger.info(f"temporary files: {size/1024/1024:.0f} MB")


def get_proteomes(url: str) -> list:
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    cur.execute(
        """
        SELECT
          P.UPID,
          P.PROTEOME_NAME,
          P.IS_REFERENCE,
          P.GC_SET_ACC,
          TO_CHAR(P.PROTEOME_TAXID),
          SN.NAME
        FROM SPTR.PROTEOME@SWPREAD P
        LEFT OUTER JOIN TAXONOMY.SPTR_STRAIN@SWPREAD S
          ON P.PROTEOME_TAXID = S.TAX_ID
        LEFT OUTER JOIN TAXONOMY.SPTR_STRAIN_NAME@SWPREAD SN
          ON S.STRAIN_ID = SN.STRAIN_ID
        WHERE P.IS_REFERENCE = 1
        """
    )

    upids = set()
    proteomes = []
    for row in cur:
        upid = row[0]

        if upid not in upids:
            upids.add(upid)
            proteomes.append({
                "id": upid,
                "name": row[1],
                "is_reference": bool(row[2]),
                "assembly": row[3],
                "tax_id": row[4],
                "strain": row[5]
            })

    cur.close()
    con.close()

    return proteomes


def get_taxa(url: str) -> list:
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    cur.execute(
        """
        SELECT
          TO_CHAR(TAX_ID), TO_CHAR(PARENT_ID), SPTR_SCIENTIFIC,
          NVL(N.SPTR_COMMON, N.NCBI_COMMON), RANK
        FROM TAXONOMY.V_PUBLIC_NODE@SWPREAD
        """
    )

    taxa = {}
    for row in cur:
        tax_id = row[0]

        taxa[tax_id] = {
            "id": tax_id,
            "parent_id": row[1],
            "sci_name": row[2],
            "full_name": row[3],
            "rank": row[4],
            "lineage": [tax_id],
            "children": set()
        }

    cur.close()
    con.close()

    for tax_id, taxon in taxa.items():
        child_id = tax_id
        parent_id = taxon["parent_id"]

        while parent_id is not None:
            parent = taxa[parent_id]
            taxon["lineage"].append(parent["id"])
            parent["children"].add(child_id)

            child_id = parent_id
            parent_id = parent["parent_id"]

    # taxa with short lineage first
    taxa = sorted(taxa.values(), key=lambda t: len(t["lineage"]))

    for taxon in taxa:
        taxon["lineage"] = list(map(str, reversed(taxon["lineage"])))
        taxon["children"] = list(taxon["children"])

    return taxa
