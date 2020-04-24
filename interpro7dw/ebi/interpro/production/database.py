# -*- coding: utf-8 -*-

from datetime import datetime
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
]


def get_databases(url: str, version: str, date: str) -> List[tuple]:
    con = cx_Oracle.connect(url)
    cur = con.cursor()

    cur.execute(
        """
        SELECT VERSION
        FROM INTERPRO.DB_VERSION
        WHERE DBCODE = 'I'
        """
    )
    db_version, = cur.fetchone()

    if db_version != version:
        cur.execute(
            """
            UPDATE INTERPRO.DB_VERSION
            SET VERSION = :1,
                FILE_DATE = :2,
                ENTRY_COUNT = (SELECT COUNT(*) 
                               FROM INTERPRO.ENTRY 
                               WHERE CHECKED='Y')
            WHERE DBCODE = 'I'
            """, (version, datetime.strptime(date, "%Y-%m-%d"))
        )
        con.commit()

    """
    Using RN=2 to join with the second most recent action
    (the most recent is the current record)
    """
    cur.execute(
        """
        SELECT
          DB.DBCODE, LOWER(DB.DBSHORT), DB.DBNAME, DB.DESCRIPTION,
          V.VERSION, V.FILE_DATE, V.ENTRY_COUNT, VA.VERSION, VA.FILE_DATE
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
        num_entries = row[6]
        previous_version = row[7]
        previous_date = row[8]

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
            num_entries,
            release_version,
            release_date,
            previous_version,
            previous_date
        ))

    cur.close()
    con.close()

    return databases
