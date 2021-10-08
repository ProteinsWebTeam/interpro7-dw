from typing import Dict

import cx_Oracle


def get_terms(url: str) -> Dict[str, tuple]:
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    cur.execute(
        """
        SELECT GT.GO_ID, GT.NAME, GT.CATEGORY, GC.TERM_NAME, GC.SORT_ORDER
        FROM GO.TERMS GT
        INNER JOIN GO.CV_CATEGORIES GC
          ON GT.CATEGORY = GC.CODE
        """
    )
    terms = {row[0]: row[1:] for row in cur}
    cur.close()
    con.close()
    return terms
