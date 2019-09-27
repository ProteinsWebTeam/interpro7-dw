# -*- coding: utf-8 -*-

import cx_Oracle


def test_database_links(url: str) -> bool:
    con = cx_Oracle.connect(url)
    cur = con.cursor()

    num_errors = 0
    for link in ("GOAPRO", "PDBE_LIVE", "SWPREAD"):
        try:
            cur.execute(f"SELECT * FROM DUAL@{link}")
        except cx_Oracle.DatabaseError:
            print(f"{link:<15} error")
            num_errors += 1
        else:
            print("f{link:<15} ok")

    cur.close()
    con.close()

    return num_errors == 0


from . import export, tables
