# -*- coding: utf-8 -*-

import MySQLdb

from interpro7dw.ebi import uniprot
from interpro7dw.ebi.interpro.utils import Table
from interpro7dw.utils import url2dict


def init(pro_url: str, stg_url: str):
    con = MySQLdb.connect(**url2dict(stg_url))
    cur = con.cursor()
    cur.execute("DROP TABLE IF EXISTS webfront_proteome")
    cur.execute(
        """
        CREATE TABLE webfront_proteome
        (
            accession VARCHAR(20) PRIMARY KEY NOT NULL,
            name VARCHAR(215) NOT NULL,
            is_reference TINYINT NOT NULL,
            strain VARCHAR(512),
            assembly VARCHAR(512),
            taxonomy_id VARCHAR(20) NOT NULL,
            counts LONGTEXT DEFAULT NULL
        ) CHARSET=utf8 DEFAULT COLLATE=utf8_unicode_ci
        """
    )
    cur.close()

    sql = """
        INSERT INTO webfront_proteome (
          accession, name, is_reference, strain, assembly, taxonomy_id
        ) VALUES (%s, %s, %s, %s, %s, %s) 
    """
    with Table(con, sql) as table:
        for record in uniprot.get_proteomes(pro_url):
            table.insert(record)

    con.commit()
    con.close()
