# -*- coding: utf-8 -*-

import json
import os

import MySQLdb

from interpro7dw.utils import url2dict


def export(url: str, outdir: str):
    con = MySQLdb.connect(**url2dict(url))
    cur = con.cursor()
    cur.execute(
        """
        SELECT version, release_date, content
        FROM webfront_release_note
        ORDER BY release_date DESC
        LIMIT 1
        """
    )
    version, date, content = cur.fetchone()

    cur.execute("SELECT COUNT(*) FROM webfront_varsplic")
    num_variants, = cur.fetchone()
    cur.close()
    con.close()

    info = json.loads(content)

    # Hash as a placeholder of day + suffix
    date_str = date.strftime("# %B %Y")
    if date.day in (1, 21, 31):
        suffix = "st"
    elif date.day in (2, 22):
        suffix = "nd"
    elif date.day in (3, 23):
        suffix = "rd"
    else:
        suffix = "th"
    date_str = date_str.replace('#', f"{date.day}{suffix}")

    with open(os.path.join(outdir, "release_notes.txt"), "wt") as fh:
        fh.write("Release Notes\n\n")
        fh.write("======================================\n\n")
        fh.write(f"Release {version}, {date_str}\n\n")

        new_entries = len(info["interpro"]["new_entries"])
        databases = list(info["member_databases"].values())
        databases.sort(key=lambda x: x["name"])
        new_databases = []
        upd_databases = []
        new_integrated = []
        for db in databases:
            if db["is_new"]:
                new_databases.append(f"{db['name']} ({db['version']})")
            elif db["is_updated"]:
                upd_databases.append(f"{db['name']} ({db['version']})")

            if db["recently_integrated"]:
                new_integrated.append((db["name"],
                                       len(db["recently_integrated"])))

        if new_entries or new_databases or upd_databases or new_integrated:
            fh.write("New features include:\n\n")

            if new_entries:
                fh.write(f"* The addition of {new_entries} "
                         f"InterPro entries.\n\n")

            if new_databases:
                fh.write(f"* New member database "
                         f"{', '.join(new_databases)}.\n\n")

            if upd_databases:
                fh.write(f"* An update to "
                         f"{', '.join(upd_databases)}.\n\n")

            if new_integrated:
                l_dbs = []
                total = 0
                for name, cnt in new_integrated:
                    l_dbs.append(f"{name} ({cnt})")
                    total += cnt

                fh.write(f"* Integration of {total} new methods "
                         f"from the {', '.join(l_dbs)} databases.\n\n")

        fh.write(
            f"""\
Contents and coverage of InterPro {version}
InterPro protein matches are now calculated for all UniProtKB and UniParc
proteins. The following statistics are for all UniProtKB proteins.
InterPro release {version} contains {info['interpro']['entries']} entries, \
representing:\n"""
        )

        for entry_type in sorted(info["interpro"]["types"]):
            cnt = info["interpro"]["types"][entry_type]
            entry_type = entry_type.replace('_', ' ').capitalize()
            fh.write(f"{entry_type:>22} {cnt:>6}\n")

        fh.write(
            f"""


Last Entry {info['interpro']['latest_entry']}

InterPro cites {info['citations']} publications in PubMed.

Member database information

{'Signature Database':>18}\
{'Version':>12}\
{'Signatures*':>25}\
{'Integrated Signatures**':>33}\n"""
        )

        for db in databases:
            fh.write(f"{db['name']:>18}"
                     f"{db['version']:>12}"
                     f"{db['signatures']:>25}"
                     f"{db['integrated_signatures']:>33}\n")

        fh.write(
            f"""
            
            
* Some signatures may not have matches to UniProtKB proteins.

** Not all signatures of a member database may be integrated at the time
of an InterPro release.

We use MobiDB-lite, a derivative of the MobiDB database, to provide consensus annotation of long-range intrinsic disorder in protein sequences.
Read more about MobiDB-lite in Bioinformatics, 33(9), 2017, 1402–1404, (doi: 10.1093/bioinformatics/btx015).


{'Sequence Database':>20}{'Version':>12}{'Count':>21}{'':16}{'Count of proteins matching':^42}
{'':69}{'any signature':^17}{'':4}{'integrated signatures':^21}\n"""
        )

        for key in ["UniProtKB", "UniProtKB/TrEMBL", "UniProtKB/Swiss-Prot"]:
            db = info["proteins"][key]
            n_p = db["count"]
            n_s = db["signatures"]
            p_s = n_s / n_p * 100
            n_is = db["integrated_signatures"]
            p_is = n_is / n_p * 100
            fh.write(f"{key:>20}"
                     f"{db['version']:>12}"
                     f"{n_p:>21}"
                     f"{'':16}"
                     f"{n_s:>9} ({p_s:.1f}%)"
                     f"{'':6}"
                     f"{n_is:>9} ({p_is:.1f}%)\n")

        num_proteins = info['proteins']['UniProtKB']['count']
        fh.write(
            f"""

Total number of proteins included in InterPro

Canonical sequences: {num_proteins}
Splice variants: {num_variants}
Total proteins: {num_proteins+num_variants}

InterPro to GO

*         Number of GO terms mapped to InterPro  - {info['interpro']['go_terms']}


Feedback
We need your help and would welcome any feedback. If you find errors or
omissions please let us know. You can contact us at:
http://www.ebi.ac.uk/support/interpro-general-query
Copyright
InterPro - Integrated Resource Of Protein Domains And Functional Sites.
Copyright (C) 2001 The InterPro Consortium. This manual and the
accompanying database may be copied and redistributed freely, without
advance permission, provided that this Copyright statement is reproduced
with each copy.\n"""
        )
