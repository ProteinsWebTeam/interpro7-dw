# -*- coding: utf-8 -*-

import MySQLdb

from interpro7dw.ebi.interpro import production
from interpro7dw.ebi.interpro.utils import Table
from . import url2dict


def init(pro_url: str, stg_url: str):
    con = MySQLdb.connect(**url2dict(stg_url))
    cur = con.cursor()
    cur.execute("DROP TABLE IF EXISTS webfront_taxonomy")
    cur.execute(
        """
        CREATE TABLE webfront_taxonomy
        (
            accession VARCHAR(20) PRIMARY KEY NOT NULL,
            scientific_name VARCHAR(255) NOT NULL,
            full_name VARCHAR(512) NOT NULL,
            lineage LONGTEXT NOT NULL,
            parent_id VARCHAR(20),
            rank VARCHAR(20) NOT NULL,
            children LONGTEXT NOT NULL,
            counts LONGTEXT DEFAULT NULL
        ) CHARSET=utf8 DEFAULT COLLATE=utf8_unicode_ci
        """
    )
    cur.close()

    sql = """
        INSERT INTO webfront_taxonomy (
          accession, scientific_name, full_name, lineage, parent_id, rank, 
          children
        ) VALUES (%s, %s, %s, %s, %s, %s, %s) 
    """
    with Table(con, sql) as table:
        for record in production.get_taxonomy(pro_url):
            table.insert(record)

    con.commit()
    con.close()
