from typing import Generator

from cx_Oracle import Cursor


def get_terms(cursor: Cursor) -> Generator[tuple, None, None]:
    cursor.execute(
        """
        SELECT
          I2G.ENTRY_AC, GT.GO_ID, GT.NAME, GT.CATEGORY, GC.TERM_NAME
        FROM INTERPRO.INTERPRO2GO I2G
        INNER JOIN GO.TERMS@GOAPRO GT
          ON I2G.GO_ID = GT.GO_ID
        INNER JOIN GO.CV_CATEGORIES@GOAPRO GC
          ON GT.CATEGORY = GC.CODE
        """
    )

    for row in cursor:
        yield row
