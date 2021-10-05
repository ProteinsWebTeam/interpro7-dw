import re
from typing import Dict, List

import cx_Oracle


def get_swissprot2enzyme(url: str) -> Dict[str, List[str]]:
    """Returns EC numbers and the Swiss-Prot proteins they are associated with.

    :param url: UniProt Oracle connection string.
    :return: A dictionary where the key is an UniProt accession
        and the value a list of EC numbers.
    """
    con = cx_Oracle.connect(url)
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

    cur.close()
    con.close()

    return proteins


def get_swissprot2reactome(url: str) -> Dict[str, List[tuple]]:
    """Returns the Reactome pathways and the Swiss-Prot proteins
    they are associated with.

    :param url: UniProt Oracle connection string.
    :return: A dictionary where the key is an UniProt accession
        and the value a list of pathways represented as tuples (ID, name).
    """
    con = cx_Oracle.connect(url)
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
