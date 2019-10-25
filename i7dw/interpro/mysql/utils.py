# -*- coding: utf-8 -*-

import re

import MySQLdb
import MySQLdb.cursors


INSERT_SIZE = 10000


def drop_index(cur: MySQLdb.cursors.Cursor, table: str, index: str):
    try:
        cur.execute(f"DROP INDEX {index} ON {table}")
    except MySQLdb.OperationalError as exc:
        code = exc.args[0]
        if code == 1091:
            # Can't drop '<index>'; check that column/key exists
            pass
        else:
            raise exc


def drop_database(url: str):
    con = MySQLdb.connect(**parse_url(url))
    cur = con.cursor()

    try:
        cur.execute("DROP DATABASE interpro")
    except MySQLdb.OperationalError as exc:
        code = exc.args[0]
        if code == 1008:
            # Can't drop database '<name>'; database doesn't exist
            pass
        else:
            raise exc
    finally:
        cur.close()
        con.close()


def init_tables(url: str):
    con = MySQLdb.connect(**parse_url(url), charset="utf8")
    cur = con.cursor()

    cur.execute('DROP TABLE IF EXISTS webfront_taxonomyperentrydb')
    cur.execute('DROP TABLE IF EXISTS webfront_taxonomyperentry')
    cur.execute('DROP TABLE IF EXISTS webfront_structure')
    cur.execute('DROP TABLE IF EXISTS webfront_alignment')
    cur.execute('DROP TABLE IF EXISTS webfront_set')
    cur.execute('DROP TABLE IF EXISTS webfront_entryannotation')
    cur.execute('DROP TABLE IF EXISTS webfront_entry')
    cur.execute('DROP TABLE IF EXISTS webfront_varsplic')
    cur.execute('DROP TABLE IF EXISTS webfront_protein')
    cur.execute('DROP TABLE IF EXISTS webfront_database')
    cur.execute('DROP TABLE IF EXISTS webfront_proteome')
    cur.execute('DROP TABLE IF EXISTS webfront_taxonomy')

    cur.execute(
        """
        CREATE TABLE webfront_database
        (
            name VARCHAR(10) NOT NULL PRIMARY KEY,
            name_long VARCHAR(25) NOT NULL,
            description LONGTEXT,
            type ENUM('protein', 'entry', 'other') NOT NULL,
            version VARCHAR(20),
            release_date DATETIME,
            prev_version VARCHAR(20),
            prev_release_date DATETIME
        ) CHARSET=utf8 DEFAULT COLLATE=utf8_unicode_ci
        """
    )

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
            overlaps_with LONGTEXT DEFAULT NULL,
            is_featured TINYINT NOT NULL DEFAULT 0,
            is_alive TINYINT NOT NULL DEFAULT 1,
            entry_date DATETIME NOT NULL,
            deletion_date DATETIME,
            counts LONGTEXT DEFAULT NULL
            # CONSTRAINT fk_entry_entry
            #   FOREIGN KEY (integrated_id)
            #   REFERENCES webfront_entry (accession),
            # CONSTRAINT fk_entry_database
            #   FOREIGN KEY (source_database)
            #   REFERENCES webfront_database (name)
        ) CHARSET=utf8 DEFAULT COLLATE=utf8_unicode_ci
        """
    )

    cur.execute(
        """
        CREATE TABLE webfront_entryannotation
        (
            annotation_id VARCHAR(255) PRIMARY KEY NOT NULL,
            accession_id VARCHAR(25) NOT NULL,
            type VARCHAR(32) NOT NULL,
            value LONGBLOB NOT NULL,
            mime_type VARCHAR(32) NOT NULL
            # CONSTRAINT fk_entryannotation_entry
            #   FOREIGN KEY (accession_id)
            #   REFERENCES webfront_entry (accession)
        ) CHARSET=utf8 DEFAULT COLLATE=utf8_unicode_ci
        """
    )

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
            # CONSTRAINT fk_taxonomy_taxonomy
            #   FOREIGN KEY (parent_id)
            #   REFERENCES webfront_taxonomy (accession)
        ) CHARSET=utf8 DEFAULT COLLATE=utf8_unicode_ci
        """
    )

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
            # CONSTRAINT fk_proteome_taxonomy
            #   FOREIGN KEY (taxonomy_id)
            #   REFERENCES webfront_taxonomy (accession)
        ) CHARSET=utf8 DEFAULT COLLATE=utf8_unicode_ci
        """
    )

    cur.execute(
        """
        CREATE TABLE webfront_protein
        (
            accession VARCHAR(15) PRIMARY KEY NOT NULL,
            identifier VARCHAR(16) NOT NULL,
            organism LONGTEXT NOT NULL,
            name VARCHAR(255) NOT NULL,
            other_names LONGTEXT NOT NULL,
            description LONGTEXT NOT NULL,
            sequence LONGTEXT NOT NULL,
            length INT(11) NOT NULL,
            size ENUM('small', 'medium', 'large') NOT NULL,
            proteome VARCHAR(20),
            gene VARCHAR(70),
            go_terms LONGTEXT NOT NULL,
            evidence_code INT(11) NOT NULL,
            source_database VARCHAR(10) NOT NULL,
            residues LONGTEXT NOT NULL,
            is_fragment TINYINT NOT NULL,
            structure LONGTEXT NOT NULL,
            tax_id VARCHAR(20) NOT NULL,
            extra_features LONGTEXT NOT NULL,
            ida_id VARCHAR(40) DEFAULT NULL,
            ida TEXT DEFAULT NULL,
            counts LONGTEXT DEFAULT NULL
            # CONSTRAINT fk_protein_taxonomy
            #   FOREIGN KEY (tax_id)
            #   REFERENCES webfront_taxonomy (accession),
            # CONSTRAINT fk_protein_database
            #   FOREIGN KEY (source_database)
            #   REFERENCES webfront_database (name),
            # CONSTRAINT fk_protein_proteome
            #   FOREIGN KEY (proteome)
            #   REFERENCES webfront_proteome (accession)
        ) CHARSET=utf8 DEFAULT COLLATE=utf8_unicode_ci
        """
    )

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

    cur.execute(
        """
        CREATE TABLE webfront_structure
        (
            accession VARCHAR(4) PRIMARY KEY NOT NULL,
            name VARCHAR(512) NOT NULL,
            short_name VARCHAR(32) DEFAULT NULL,
            other_names LONGTEXT DEFAULT NULL,
            source_database VARCHAR(10) NOT NULL,
            experiment_type VARCHAR(16) NOT NULL,
            release_date DATETIME NOT NULL,
            resolution FLOAT DEFAULT NULL,
            literature LONGTEXT NOT NULL,
            chains LONGTEXT NOT NULL,
            proteins LONGTEXT NOT NULL,
            secondary_structures LONGTEXT NOT NULL,
            counts LONGTEXT DEFAULT NULL
            # CONSTRAINT fk_structure_database
            #   FOREIGN KEY (source_database)
            #   REFERENCES webfront_database (name)
        ) CHARSET=utf8 DEFAULT COLLATE=utf8_unicode_ci
        """
    )

    cur.execute(
        """
        CREATE TABLE webfront_set
        (
            accession VARCHAR(20) PRIMARY KEY NOT NULL,
            name VARCHAR(100),
            description TEXT,
            source_database VARCHAR(10) NOT NULL,
            is_set TINYINT NOT NULL,
            relationships LONGTEXT NOT NULL,
            integrated LONGTEXT DEFAULT NULL,
            counts LONGTEXT DEFAULT NULL
            # CONSTRAINT fk_set_database
            #   FOREIGN KEY (source_database)
            #   REFERENCES webfront_database (name)
        ) CHARSET=utf8 DEFAULT COLLATE=utf8_unicode_ci
        """
    )

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
            # CONSTRAINT fk_alignment_set
            #   FOREIGN KEY (set_acc)
            #   REFERENCES webfront_set (accession)
        ) CHARSET=utf8 DEFAULT COLLATE=utf8_unicode_ci
        """
    )

    cur.execute(
        """
        CREATE TABLE webfront_taxonomyperentry
        (
          id INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
          tax_id VARCHAR(20) NOT NULL,
          entry_acc VARCHAR(25) NOT NULL,
          counts LONGTEXT NULL NULL
          # CONSTRAINT fk_taxonomy_entry_tax
          #   FOREIGN KEY (tax_id)
          #   REFERENCES webfront_taxonomy (accession),
          # CONSTRAINT fk_taxonomy_entry_entry
          #   FOREIGN KEY (entry_acc)
          #   REFERENCES webfront_entry (accession)
        ) CHARSET=utf8 DEFAULT COLLATE=utf8_unicode_ci
        """
    )

    cur.execute(
        """
        CREATE TABLE webfront_taxonomyperentrydb
        (
          id INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
          tax_id VARCHAR(20) NOT NULL,
          source_database VARCHAR(10) NOT NULL,
          counts LONGTEXT NOT NULL
          # CONSTRAINT fk_taxonomy_database_tax
          #   FOREIGN KEY (tax_id)
          #   REFERENCES webfront_taxonomy (accession),
          # CONSTRAINT fk_taxonomy_database_db
          #   FOREIGN KEY (source_database)
          #   REFERENCES webfront_database (name)
        ) CHARSET=utf8 DEFAULT COLLATE=utf8_unicode_ci
        """
    )

    cur.execute(
        """
        CREATE TABLE IF NOT EXISTS webfront_release_note
        (
            version VARCHAR(20) PRIMARY KEY NOT NULL,
            release_date DATETIME NOT NULL,
            content LONGTEXT NOT NULL
        ) CHARSET=utf8 DEFAULT COLLATE=utf8_unicode_ci
        """
    )

    cur.close()
    con.close()


def add_foreign_keys(url: str):
    con = MySQLdb.connect(**parse_url(url), charset="utf8")
    cur = con.cursor()
    cur.execute(
        """
        ALTER TABLE webfront_entry 
        ADD CONSTRAINT fk_entry_entry_integrated
        FOREIGN KEY (integrated_id) REFERENCES webfront_entry (accession)
        """
    )
    cur.execute(
        """
        ALTER TABLE webfront_entry 
        ADD CONSTRAINT fk_entry_database_source
        FOREIGN KEY (source_database) REFERENCES webfront_database (name)
        """
    )
    cur.execute(
        """
        ALTER TABLE webfront_entryannotation 
        ADD CONSTRAINT fk_entryannotation_entry_accession
        FOREIGN KEY (accession_id) REFERENCES webfront_entry (accession)
        """
    )
    cur.execute(
        """
        ALTER TABLE webfront_taxonomy 
        ADD CONSTRAINT fk_taxonomy_taxonomy_parent
        FOREIGN KEY (parent_id) REFERENCES webfront_taxonomy (accession)
        """
    )
    cur.execute(
        """
        ALTER TABLE webfront_proteome 
        ADD CONSTRAINT fk_proteome_taxonomy_taxonomy
        FOREIGN KEY (taxonomy_id) REFERENCES webfront_taxonomy (accession)
        """
    )
    cur.execute(
        """
        ALTER TABLE webfront_protein 
        ADD CONSTRAINT fk_protein_taxonomy_tax
        FOREIGN KEY (tax_id) REFERENCES webfront_taxonomy (accession)
        """
    )
    cur.execute(
        """
        ALTER TABLE webfront_protein 
        ADD CONSTRAINT fk_protein_database_source
        FOREIGN KEY (source_database) REFERENCES webfront_database (name)
        """
    )
    cur.execute(
        """
        ALTER TABLE webfront_protein 
        ADD CONSTRAINT fk_protein_proteome_proteome
        FOREIGN KEY (proteome) REFERENCES webfront_proteome (accession)
        """
    )
    cur.execute(
        """
        ALTER TABLE webfront_structure 
        ADD CONSTRAINT fk_structure_database_source
        FOREIGN KEY (source_database) REFERENCES webfront_database (name)
        """
    )
    cur.execute(
        """
        ALTER TABLE webfront_set 
        ADD CONSTRAINT fk_set_database_source
        FOREIGN KEY (source_database) REFERENCES webfront_database (name)
        """
    )
    cur.execute(
        """
        ALTER TABLE webfront_alignment 
        ADD CONSTRAINT fk_alignment_set_set
        FOREIGN KEY (set_acc) REFERENCES webfront_set (accession)
        """
    )
    cur.execute(
        """
        ALTER TABLE webfront_taxonomyperentry 
        ADD CONSTRAINT fk_taxonomyentry_taxonomy_tax
        FOREIGN KEY (tax_id) REFERENCES webfront_taxonomy (accession)
        """
    )
    cur.execute(
        """
        ALTER TABLE webfront_taxonomyperentry 
        ADD CONSTRAINT fk_taxonomyentry_entry_entry
        FOREIGN KEY (entry_acc) REFERENCES webfront_entry (accession)
        """
    )
    cur.execute(
        """
        ALTER TABLE webfront_taxonomyperentrydb 
        ADD CONSTRAINT fk_taxonomydatabase_taxonomy_tax
        FOREIGN KEY (tax_id) REFERENCES webfront_taxonomy (accession)
        """
    )
    cur.execute(
        """
        ALTER TABLE webfront_taxonomyperentrydb 
        ADD CONSTRAINT fk_taxonomydatabase_database_source
        FOREIGN KEY (source_database) REFERENCES webfront_database (name)
        """
    )


def parse_url(url: str) -> dict:
    m = re.match(r'([^/]+)/([^@]+)@([^:]+):(\d+)/(\w+)', url)

    if m is None:
        raise RuntimeError(f"invalid connection string: {url}")

    return dict(
        user=m.group(1),
        passwd=m.group(2),
        host=m.group(3),
        port=int(m.group(4)),
        db=m.group(5)
    )


def reduce(src: dict) -> dict:
    dst = {}
    for k, v in src.items():
        if isinstance(v, dict):
            dst[k] = reduce(v)
        elif isinstance(v, (list, set, tuple)):
            dst[k] = len(v)
        else:
            dst[k] = v

    return dst
