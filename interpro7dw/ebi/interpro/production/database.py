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
FEATURE_DATABASES = [
    'g',    # MobiDB Lite
    'j',    # Phobius
    'n',    # Signal Euk
    'q',    # TMHMM
    's',    # SignalP Gram positive
    'v',    # SignalP Gram negative
    'x',    # COILS
]


def get_databases(url: str, version: str, date: str,
                  update: bool = False) -> List[tuple]:
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

    """
    use_db:
      * True if the upcoming release is in DB_VERSION and the current one 
        is the most recent in DB_VERSION_AUDIT
      * False if the current release is in DB_VERSION and the most recent
        in DB_VERSION_AUDIT is the previous release
    """
    if db_version != version:
        if update:
            cur.execute(
                """
                UPDATE INTERPRO.DB_VERSION
                SET VERSION = :1,
                    FILE_DATE = :2,
                    ENTRY_COUNT = (SELECT COUNT(*) 
                                   FROM INTERPRO.ENTRY 
                                   WHERE CHECKED = 'Y')
                WHERE DBCODE = 'I'
                """, (version, datetime.strptime(date, "%Y-%m-%d"))
            )
            con.commit()
            use_db = True
        else:
            use_db = False
    else:
        use_db = True

    cur.execute("SELECT COUNT(*) FROM INTERPRO.ENTRY WHERE CHECKED = 'Y'")
    num_interpro_entries, = cur.fetchone()

    """
    Using RN=2 to join with the second most recent action in DB_VERSION_AUDIT
    (the most recent is the same record as in DB_VERSION)
    """
    cur.execute(
        """
        SELECT
          DB.DBCODE, LOWER(DB.DBSHORT), DB.DBSHORT, DB.DBNAME, 
          DB.DESCRIPTION, V.VERSION, V.FILE_DATE, V.ENTRY_COUNT, VA.VERSION, 
          VA.FILE_DATE
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
        name = row[1]
        name_alt = row[2]
        name_long = row[3]
        description = row[4]
        release_version = row[5]
        release_date = row[6]
        num_entries = row[7]
        previous_version = row[8]
        previous_date = row[9]

        if code in ENTRY_DATABASES:
            db_type = "entry"
        elif code in FEATURE_DATABASES:
            db_type = "feature"
        elif code in ('S', 'T', 'u'):
            if code == 'S':
                name = "reviewed"
            elif code == 'T':
                name = "unreviewed"

            db_type = "protein"
        else:
            db_type = "other"

        if code == 'I' and not use_db:
            num_entries = num_interpro_entries
            previous_version = release_version
            previous_date = release_date
            release_date = datetime.strptime(date, "%Y-%m-%d")
            release_version = version

        databases.append((
            name,
            name_alt,
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
