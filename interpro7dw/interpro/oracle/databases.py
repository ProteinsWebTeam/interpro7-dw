import datetime
import pickle

import oracledb

from interpro7dw.uniprot import goa


def export(ipr_uri: str, goa_uri: str, version: str, date: datetime.date,
           output: str, update: bool = False):
    """Exports information on databases/data sources used in InterPro.

    :param ipr_uri: The InterPro Oracle connection string.
    :param goa_uri: The GOA Oracle connection string.
    :param version: The version of the upcoming InterPro release.
    :param date: The date of the upcoming InterPro release.
    :param output: The output file.
    :param update: If True, update the production table.
    """
    con = oracledb.connect(ipr_uri)
    cur = con.cursor()

    cur.execute("SELECT COUNT(*) FROM INTERPRO.ENTRY WHERE CHECKED = 'Y'")
    num_interpro_entries, = cur.fetchone()

    cur.execute("SELECT VERSION FROM INTERPRO.DB_VERSION WHERE DBCODE = 'I'")
    prod_version, = cur.fetchone()

    signature_dbcodes, feature_dbcodes = get_databases_codes(cur)

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
            """,
            [version, date, num_interpro_entries]
        )
        con.commit()
    else:
        # DB_VERSION is outdated and will stay outdated
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
          DB.DBCODE, DB.DBSHORT, DB.DBNAME, DB.DESCRIPTION, V.VERSION, 
          V.FILE_DATE, V.ENTRY_COUNT
        FROM INTERPRO.CV_DATABASE DB
        LEFT OUTER JOIN INTERPRO.DB_VERSION V ON DB.DBCODE = V.DBCODE
        -- Ignore unused databases: COG, OMIM
        WHERE DB.DBCODE NOT IN ('O', 'o')
        """
    )

    databases = {}
    for rec in cur:
        dbcode = rec[0]
        short_name = rec[1]
        name = rec[2]
        description = rec[3]
        release_version = rec[4]
        release_date = rec[5]
        num_entries = rec[6]
        prev_release_version = prev_release_date = None

        try:
            prev_releases = all_releases[dbcode]
        except KeyError:
            pass
        else:
            for dbversion, dbdate in prev_releases:
                if dbversion != release_version:
                    prev_release_version = dbversion
                    prev_release_date = dbdate
                    break

        if dbcode == 'I' or dbcode in signature_dbcodes:
            db_type = "entry"
        elif dbcode in feature_dbcodes:
            db_type = "feature"
        elif dbcode in ['S', 'T', 'u']:
            if dbcode == 'S':
                short_name = "reviewed"
            elif dbcode == 'T':
                short_name = "unreviewed"

            db_type = "protein"
        else:
            db_type = "other"

        if dbcode == 'I' and not use_db_version:
            # DB_VERSION is outdated:
            # it contains info for the live release (soon to be 'previous')
            num_entries = num_interpro_entries
            prev_release_version = release_version
            prev_release_date = release_date
            release_version = version
            release_date = date

        databases[short_name] = {
            "name": name,
            "description": description,
            "type": db_type,
            "entries": num_entries,
            "release": {
                "version": release_version,
                "date": release_date,
            },
            "previous_release": {
                "version": prev_release_version,
                "date": prev_release_date,
            },
        }

    cur.close()
    con.close()

    go_timestamp = goa.get_timestamp(goa_uri)
    databases["go"] = {
        "name": "Gene Ontology",
        "description": "The Gene Ontology (GO) describes knowledge of the "
                       "biological domain with respect to three aspects: "
                       "Molecular function, Biological process, "
                       "and Cellular component.",
        "type": "other",
        "entries": len(goa.get_terms(goa_uri)),
        "release": {
            "version": go_timestamp.strftime("%Y-%m-%d"),
            "date": go_timestamp,
        },
        "previous_release": {
            "version": None,
            "date": None,
        },
    }

    with open(output, "wb") as fh:
        pickle.dump(databases, fh)


def get_databases_codes(cur: oracledb.Cursor) -> tuple[list[str], list[str]]:
    cur.execute("SELECT DISTINCT DBCODE FROM INTERPRO.METHOD")
    signatures = {dbcode for dbcode, in cur}
    signatures.add("a")  # Flag AntiFam as a member database

    cur.execute("SELECT DISTINCT DBCODE FROM INTERPRO.FEATURE_METHOD")
    features = {dbcode for dbcode, in cur}

    return list(signatures), list(features - signatures)
