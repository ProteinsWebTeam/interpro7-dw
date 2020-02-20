# -*- coding: utf-8 -*-

import json

import MySQLdb

from interpro7dw.ebi.interpro import production as ippro
from interpro7dw.ebi.interpro.utils import Table
from interpro7dw.utils import dataload, url2dict


def insert_isoforms(src_entries: str, pro_url: str, stg_url: str):
    entries = dataload(src_entries)

    con = MySQLdb.connect(**url2dict(stg_url))
    cur = con.cursor()
    cur.execute("DROP TABLE IF EXISTS webfront_varsplic")
    cur.execute(
        """
        CREATE TABLE webfront_varsplic
        (
            accession VARCHAR(20) PRIMARY KEY NOT NULL,
            protein_acc VARCHAR(15) NOT NULL,
            length INT(11) NOT NULL,
            sequence LONGTEXT NOT NULL,
            features LONGTEXT NOT NULL
        ) CHARSET=utf8 DEFAULT COLLATE=utf8_unicode_ci
        """
    )

    sql = """
            INSERT INTO webfront_varsplic VALUES (%s, %s, %s, %s, %s) 
        """
    with Table(con, sql) as table:
        for obj in ippro.get_isoforms(pro_url):
            accession = obj[0]
            protein_acc = obj[1]
            length = obj[2]
            sequence = obj[3]
            features = obj[4]

            enriched_features = {}
            for entry_acc, locations in features.items():
                entry = entries[entry_acc]

                enriched_features[entry_acc] = {
                    "accession": entry_acc,
                    "integrated": entry.integrated_in,
                    "name": entry.name,
                    "type": entry.type,
                    "source_database": entry.database,
                    "locations": locations
                }

            table.insert((accession, protein_acc, length, sequence,
                          json.dumps(enriched_features)))

    con.commit()

    cur = con.cursor()
    cur.execute("CREATE INDEX i_webfront_varsplic "
                "ON webfront_varsplic (protein_acc)")
    cur.close()
    con.close()
