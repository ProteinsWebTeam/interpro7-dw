#!/usr/bin/env python

from email.message import EmailMessage
from smtplib import SMTP

import cx_Oracle

from .clan import get_clans, iter_clan_alignments
from .database import get_databases
from .entry import export_entries, get_features, get_signatures
from .protein import chunk_proteins
from .protein import export_features
from .protein import export_matches
from .protein import export_proteins
from .protein import export_residues
from .protein import export_sequences
from .protein import get_isoforms
from .taxonomy import export_taxonomy


def unfreeze(host: str, port: int, addr: str):
    msg = EmailMessage()
    msg.set_content("""\
Dear curators,

All steps from the InterPro release procedure requiring data from \
the production database successfully completed.

You may resume integrating signatures, and creating/annotating entries. \
Have fun!

The InterPro Production Team
""")
    msg["Sender"] = addr
    msg["To"] = addr
    msg["Subject"] = "Database unfrozen"

    with SMTP(host, port=port) as s:
        s.send_message(msg)


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
