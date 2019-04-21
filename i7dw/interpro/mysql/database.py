from datetime import datetime

from .. import oracle
from ... import dbms


def insert_databases(ora_uri: str, my_uri: str, version: str,
                     release_date: str):
    data = []
    for db in oracle.get_databases(ora_uri):
        if db["name"] == "interpro" and db["version"]["code"] != version:
            """
            Happens when Oracle hasn't been updated yet
            (DB_VERSION still on the previous release)

            --> move the `version` to `previous_version`
            """
            db["previous_version"] = db["version"]
            db["version"] = {
                "code": version,
                "date": datetime.strptime(release_date, "%Y-%m-%d"),
            }

        data.append((
            db["name"],
            db["name_long"],
            db["description"],
            db["type"],
            db["version"]["code"],
            db["version"]["date"],
            db["previous_version"]["code"],
            db["previous_version"]["date"]
        ))

    con, cur = dbms.connect(my_uri)
    cur.executemany(
        """
        INSERT INTO webfront_database (
          name, name_long, description, type,
          version, release_date, prev_version, prev_release_date
        )
        VALUES (%s, %s, %s, %s, %s, %s, %s, %s)
        """,
        data
    )
    con.commit()
    cur.close()
    con.close()


def get_databases(uri: str) -> dict:
    con, cur = dbms.connect(uri)

    cur.execute(
        """
        SELECT name, name_long, version, release_date
        FROM webfront_database WHERE type = 'entry'
        """
    )
    databases = {}
    for name, name_long, version, release_date in cur:
        databases[name] = {
            "name_long": name_long,
            "version": version,
            "release_date": release_date
        }

    cur.close()
    con.close()

    return databases
