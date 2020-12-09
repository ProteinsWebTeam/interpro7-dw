#!/usr/bin/env python

import cx_Oracle

from .clan import export_clans
from .database import get_databases
from .entry import export_entries, get_features, get_hmms, get_signatures
from .protein import chunk_proteins
from .protein import export_features
from .protein import export_matches
from .protein import export_proteins
from .protein import export_residues, export_residues2
from .protein import export_sequences
from .protein import get_isoforms
from .taxonomy import export_taxonomy


def test_db_links(url: str) -> bool:
    success = True

    for link in ["GOAPRO", "INTACPRO", "PDBE_LIVE", "SWPREAD"]:
        # New connection to prevent ORA-02020 (too many database links in use)
        con = cx_Oracle.connect(url)
        cur = con.cursor()

        try:
            cur.execute(f"SELECT * FROM DUAL@{link}")
        except cx_Oracle.DatabaseError:
            print(f"{link:<20} error")
            success = False
        else:
            print(f"{link:<20} ok")
        finally:
            cur.close()
            con.close()

    return success
