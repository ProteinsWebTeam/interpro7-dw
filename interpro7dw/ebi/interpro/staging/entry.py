# -*- coding: utf-8 -*-

import json

import MySQLdb

from interpro7dw.ebi import pfam
from interpro7dw.ebi.interpro import production as ippro
from interpro7dw.ebi.interpro.utils import Table
from interpro7dw.utils import datadump, url2dict


def init_entries(pro_url: str, stg_url: str):
    con = MySQLdb.connect(**url2dict(stg_url))
    cur = con.cursor()
    cur.execute("DROP TABLE IF EXISTS webfront_entry")
    cur.execute(
        """
        CREATE TABLE webfront_entry
        (
            entry_id VARCHAR(10) DEFAULT NULL,
            accession VARCHAR(25) PRIMARY KEY NOT NULL,
            type VARCHAR(50) NOT NULL,
            name LONGTEXT,
            short_name VARCHAR(100),
            source_database VARCHAR(10) NOT NULL,
            member_databases LONGTEXT,
            integrated_id VARCHAR(25),
            go_terms LONGTEXT,
            description LONGTEXT,
            wikipedia LONGTEXT,
            literature LONGTEXT,
            hierarchy LONGTEXT,
            cross_references LONGTEXT,
            interactions LONGTEXT,
            pathways LONGTEXT DEFAULT NULL,
            overlaps_with LONGTEXT DEFAULT NULL,
            is_featured TINYINT NOT NULL DEFAULT 0,
            is_alive TINYINT NOT NULL DEFAULT 1,
            entry_date DATETIME NOT NULL,
            history LONGTEXT,
            deletion_date DATETIME,
            counts LONGTEXT DEFAULT NULL
        ) CHARSET=utf8 DEFAULT COLLATE=utf8_unicode_ci
        """
    )
    cur.close()

    sql = """
    
    """

    con.comit()
    con.close()


def insert_annotations(pfam_url: str, stg_url: str):
    con = MySQLdb.connect(**url2dict(stg_url))
    cur = con.cursor()
    cur.execute("DROP TABLE IF EXISTS webfront_entryannotation")
    cur.execute(
        """
        CREATE TABLE webfront_entryannotation
        (
            annotation_id VARCHAR(255) PRIMARY KEY NOT NULL,
            accession_id VARCHAR(25) NOT NULL,
            type VARCHAR(32) NOT NULL,
            value LONGBLOB NOT NULL,
            mime_type VARCHAR(32) NOT NULL
        ) CHARSET=utf8 DEFAULT COLLATE=utf8_unicode_ci
        """
    )
    cur.close()

    sql = """
        INSERT INTO webfront_entryannotation
        VALUES (%s, %s, %s, %s, %s)
    """

    with Table(con, sql) as table:
        for acc, anno_type, value, mime in pfam.get_annotations(pfam_url):
            anno_id = f"{acc}--{anno_type}"
            table.insert((anno_id, acc, anno_type, value, mime))

    con.commit()
    con.close()


def init_sets(pro_url: str, stg_url: str, output: str, threshold: float=1e-2):
    sets = ippro.get_sets(pro_url)

    con = MySQLdb.connect(**url2dict(stg_url))
    cur = con.cursor()
    cur.execute("DROP TABLE IF EXISTS webfront_alignment")
    cur.execute(
        """
        CREATE TABLE webfront_alignment
        (
            id INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
            set_acc VARCHAR(20) NOT NULL,
            entry_acc VARCHAR(25) NOT NULL,
            target_acc VARCHAR(25) NOT NULL,
            target_set_acc VARCHAR(20),
            score DOUBLE NOT NULL,
            seq_length MEDIUMINT NOT NULL,
            domains TEXT NOT NULL
        ) CHARSET=utf8 DEFAULT COLLATE=utf8_unicode_ci
        """
    )
    cur.close()

    sql = """
        INSERT INTO webfront_alignment (set_acc, entry_acc, target_acc, 
                                        target_set_acc, score, seq_length, 
                                        domains)
        VALUES (%s, %s, %s, %s, %s, %s, %s)
    """
    with Table(con, sql) as table:
        for aln in ippro.get_set_alignments(pro_url):
            set_acc = aln[0]
            entry_acc = aln[1]
            seq_length = aln[2]
            target_acc = aln[3]
            target_set_acc = aln[4]
            score = aln[5]
            domains = aln[6]

            if score > threshold:
                continue

            table.insert((set_acc, entry_acc, target_acc, target_set_acc,
                          score, seq_length, json.dumps(domains)))

            if set_acc == target_set_acc:
                # Query and target from the same set: update the set's links
                sets[set_acc].add_link(entry_acc, target_acc, score)

    con.commit()
    con.close()

    datadump(output, sets)
