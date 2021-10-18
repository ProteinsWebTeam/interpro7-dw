from interpro7dw.utils import logger
from interpro7dw.utils.oracle import lob_as_str
from interpro7dw.utils.store import SimpleStore

import cx_Oracle


def export_structural_models(url: str, output: str):
    logger.info("exporting structural models")
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    cur.outputtypehandler = lob_as_str
    cur.execute(
        """
        SELECT METHOD_AC, LOWER(ALGORITHM), CONTACTS, ERRORS, PLDDT, STRUCTURE
        FROM INTERPRO.STRUCT_MODEL
        """
    )

    with SimpleStore(output) as store:
        for record in store:
            store.add(record)

    cur.close()
    con.close()

    logger.info("done")
