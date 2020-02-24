# -*- coding: utf-8 -*-

from typing import List

import cx_Oracle


ENTRY_DATABASES = [
    'B',    # SFLD
    'F',    # PRINTS
    'H',    # Pfam
    'I',    # InterPro
    'J',    # CDD
    'M',    # PROSITE profiles
    'N',    # TIGRFAMs
    'P',    # PROSITE patterns
    'Q',    # HAMAP
    'R',    # SMART
    'U',    # PIRSF
    'V',    # PANTHER
    'X',    # CATH-Gene3D
    'Y',    # SUPERFAMILY
    'g',    # MobiDB Lite
]


def get_databases(url: str) -> List[tuple]:
    con = cx_Oracle.connect(url)
    cur = con.cursor()

    """
    Using RN=2 to join with the second most recent action
    (the most recent is the current record)
    """
    cur.execute(
        """
        SELECT
          DB.DBCODE, LOWER(DB.DBSHORT), DB.DBNAME, DB.DESCRIPTION,
          V.VERSION, V.FILE_DATE, VA.VERSION, VA.FILE_DATE
        FROM INTERPRO.CV_DATABASE DB
        LEFT OUTER JOIN INTERPRO.DB_VERSION V ON DB.DBCODE = V.DBCODE
        LEFT OUTER JOIN (
          SELECT
            DBCODE, VERSION, FILE_DATE,
            ROW_NUMBER() OVER (
              PARTITION BY DBCODE ORDER BY TIMESTAMP DESC
            ) RN
          FROM INTERPRO.DB_VERSION_AUDIT
          WHERE ACTION = 'U'
        ) VA ON DB.DBCODE = VA.DBCODE AND VA.RN = 2
        """
    )

    databases = []
    for row in cur:
        code = row[0]
        name_short = row[1]
        name_long = row[2]
        description = row[3]
        release_version = row[4]
        release_date = row[5]
        previous_version = row[6]
        previous_date = row[7]

        if code in ENTRY_DATABASES:
            db_type = "entry"
        elif code in ('S', 'T', 'u'):
            if code == 'S':
                name_short = "reviewed"
            elif code == 'T':
                name_short = "unreviewed"

            db_type = "protein"
        else:
            db_type = "other"

        databases.append((
            name_short,
            name_long,
            description,
            db_type,
            release_version,
            release_date,
            previous_version,
            previous_date
        ))

    cur.close()
    con.close()

    return databases
