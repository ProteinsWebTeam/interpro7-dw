import cx_Oracle


def update_pathways(uri: str, entry2pathways: dict[str, list[tuple]]):
    con = cx_Oracle.connect(uri)
    cur = con.cursor()
    cur.execute("TRUNCATE TABLE INTERPRO.ENTRY2PATHWAY")
    cur.execute(
        """
        SELECT LOWER(DBSHORT), DBCODE
        FROM INTERPRO.CV_DATABASE
        """
    )
    id2dbcode = dict(cur.fetchall())

    params = []
    for entry_acc, databases in entry2pathways.values():
        for database, pathways in databases.items():
            dbcode = id2dbcode[database.lower()]

            for pathway_id, pathway_name in pathways:
                params.append((entry_acc, dbcode, pathway_id, pathway_name))

    for i in range(0, len(params), 1000):
        cur.executemany(
            """
            INSERT INTO INTERPRO.ENTRY2PATHWAY (ENTRY_AC, DBCODE, AC, NAME)
            VALUES (:1, :2, :3, :4)
            """,
            params[i:i+1000]
        )

    con.commit()
    cur.close()
    con.close()


def export_entries(interpro_uri: str, goa_uri: str, intact_uri: str):
    pass
