import pickle

import MySQLdb

from interpro7dw.utils import logger
from interpro7dw.utils.mysql import uri2dict


def populate(uri: str, databases_file: str):
    logger.info("creating webfront_database")
    with open(databases_file, "rb") as fh:
        databases = pickle.load(fh)

    params = []
    for key, info in databases.items():
        params.append((
            key.lower(),
            key,
            info["name"],
            info["description"],
            info["type"],
            info["entries"],
            info["release"]["version"],
            info["release"]["date"],
            info["previous_release"]["version"],
            info["release"]["date"]
        ))

    con = MySQLdb.connect(**uri2dict(uri), charset="utf8mb4")
    cur = con.cursor()
    cur.execute("DROP TABLE IF EXISTS webfront_database")
    cur.execute(
        """
        CREATE TABLE webfront_database
        (
            name VARCHAR(10) NOT NULL PRIMARY KEY,
            name_alt VARCHAR(10) NOT NULL,
            name_long VARCHAR(25) NOT NULL,
            description LONGTEXT,
            type ENUM('protein', 'entry', 'feature', 'other') NOT NULL,
            num_entries INTEGER,
            version VARCHAR(20),
            release_date DATETIME,
            prev_version VARCHAR(20),
            prev_release_date DATETIME
        ) CHARSET=utf8mb4 DEFAULT COLLATE=utf8mb4_unicode_ci
        """
    )

    cur.executemany(
        """
        INSERT INTO webfront_database 
        VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
        """,
        params
    )
    con.commit()
    cur.close()
    con.close()

    logger.info("done")
