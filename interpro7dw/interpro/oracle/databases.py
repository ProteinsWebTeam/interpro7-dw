import pickle
from datetime import datetime

import cx_Oracle


def export(uri: str, version: str, date: str, file: str, update: bool = False):
    """Exports information on databases/data sources used in InterPro.

    :param uri: The Oracle connection string.
    :param version: The version of the upcoming InterPro release.
    :param date: The date of the upcoming InterPro release (YYYY-MM-DD).
    :param file: The output file.
    :param update: If True, update the production table.
    """
    con = cx_Oracle.connect(uri)
    cur = con.cursor()

    cur.execute("SELECT COUNT(*) FROM INTERPRO.ENTRY WHERE CHECKED = 'Y'")
    num_interpro_entries, = cur.fetchone()

    cur.execute("SELECT VERSION FROM INTERPRO.DB_VERSION WHERE DBCODE = 'I'")
    prod_version, = cur.fetchone()

    cur.execute("SELECT DISTINCT DBCODE FROM INTERPRO.METHOD")
    signature_dbcodes = {dbcode for dbcode, in cur}

    cur.execute("SELECT DISTINCT DBCODE FROM INTERPRO.FEATURE_METHOD")
    feature_dbcodes = {dbcode for dbcode, in cur}

    if prod_version == version:
        # DB_VERSION is already up-to-date
        use_db_version = True
    elif update:
        # DB_VERSION is outdated, but will be up-to-date
        use_db_version = True
        cur.execute(
            """
            UPDATE INTERPRO.DB_VERSION
            SET VERSION = :1,
                FILE_DATE = :2,
                ENTRY_COUNT = :3
            WHERE DBCODE = 'I'
            """, (version, datetime.strptime(date, "%Y-%m-%d"),
                  num_interpro_entries)
        )
        con.commit()
    else:
        # DB_VERSION is outdated and will stay outdated
        # This run is a test done on the production database (SRSLY?!!111)
        use_db_version = False

    # Get all releases
    all_releases = {}
    cur.execute(
        """
        SELECT DBCODE, VERSION, FILE_DATE
        FROM INTERPRO.DB_VERSION_AUDIT
        ORDER BY TIMESTAMP DESC
        """
    )

    for dbcode, dbversion, dbdate in cur:
        try:
            all_releases[dbcode].append((dbversion, dbdate))
        except KeyError:
            all_releases[dbcode] = [(dbversion, dbdate)]

    cur.execute(
        """
        SELECT
          DB.DBCODE, LOWER(DB.DBSHORT), DB.DBSHORT, DB.DBNAME,
          DB.DESCRIPTION, V.VERSION, V.FILE_DATE, V.ENTRY_COUNT
        FROM INTERPRO.CV_DATABASE DB
        LEFT OUTER JOIN INTERPRO.DB_VERSION V ON DB.DBCODE = V.DBCODE
        """
    )

    databases = []
    for rec in cur:
        code = rec[0]
        identifier = rec[1]
        short_name = rec[2]
        name = rec[3]
        description = rec[4]
        release_version = rec[5]
        release_date = rec[6]
        num_entries = rec[7]
        prev_release_version = prev_release_date = None

        try:
            prev_releases = all_releases[code]
        except KeyError:
            pass
        else:
            for dbversion, dbdate in prev_releases:
                if dbversion != release_version:
                    prev_release_version = dbversion
                    prev_release_date = dbdate
                    break

        if code == 'I' or code in signature_dbcodes:
            db_type = "entry"
        elif code in feature_dbcodes:
            db_type = "feature"
        elif code in ['S', 'T', 'u']:
            if code == 'S':
                identifier = "reviewed"
            elif code == 'T':
                identifier = "unreviewed"

            db_type = "protein"
        else:
            db_type = "other"

        if code == 'I' and not use_db_version:
            # DB_VERSION is outdated:
            # it contains info for the live release (soon to be 'previous')
            num_entries = num_interpro_entries
            prev_release_version = release_version
            prev_release_date = release_date
            release_version = version
            release_date = datetime.strptime(date, "%Y-%m-%d")

        databases.append((
            identifier,
            short_name,
            name,
            description,
            db_type,
            num_entries,
            release_version,
            release_date,
            prev_release_version,
            prev_release_date
        ))

    cur.close()
    con.close()

    with open(file, "wb") as fh:
        pickle.dump(databases, fh)
