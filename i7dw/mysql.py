#!/usr/bin/env python
# -*- coding: utf-8 -*-

import json
import logging
import time

from . import dbms, disk
from .ebi import interpro, pdbe, uniprot


logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s: %(levelname)s: %(message)s',
    datefmt='%y-%m-%d %H:%M:%S'
)


def init(uri):
    con, cur = dbms.connect(uri)

    cur.execute('DROP TABLE IF EXISTS webfront_structure')
    cur.execute('DROP TABLE IF EXISTS webfront_set')
    cur.execute('DROP TABLE IF EXISTS webfront_entryannotation')
    cur.execute('DROP TABLE IF EXISTS webfront_entry')
    cur.execute('DROP TABLE IF EXISTS webfront_protein')
    cur.execute('DROP TABLE IF EXISTS webfront_database')
    cur.execute('DROP TABLE IF EXISTS webfront_proteome')
    cur.execute('DROP TABLE IF EXISTS webfront_taxonomy')

    cur.execute(
        """
        CREATE TABLE webfront_database
        (
            name VARCHAR(10) NOT NULL PRIMARY KEY,
            name_long VARCHAR(20) NOT NULL,
            description LONGTEXT,
            version VARCHAR(20),
            release_date DATETIME,
            type ENUM('protein', 'entry', 'other') NOT NULL,
            prev_version VARCHAR(20),
            prev_release_date DATETIME
        ) CHARSET=utf8 DEFAULT COLLATE utf8_unicode_ci
        """
    )

    cur.execute(
        """
        CREATE TABLE webfront_entry
        (
            entry_id VARCHAR(10),
            accession VARCHAR(25) PRIMARY KEY NOT NULL,
            type VARCHAR(50) NOT NULL,
            name LONGTEXT,
            short_name VARCHAR(100) NOT NULL,
            source_database VARCHAR(10) NOT NULL,
            member_databases LONGTEXT NOT NULL,
            integrated_id VARCHAR(25),
            go_terms LONGTEXT NOT NULL,
            description LONGTEXT NOT NULL,
            wikipedia LONGTEXT NOT NULL,
            literature LONGTEXT NOT NULL,
            hierarchy LONGTEXT NOT NULL,
            cross_references LONGTEXT NOT NULL,
            entry_date DATETIME NOT NULL,
            overlaps_with LONGTEXT NOT NULL,
            is_featured TINYINT NOT NULL DEFAULT 0,
            CONSTRAINT fk_webfront_entry_webfront_entry_integrated_id FOREIGN KEY (integrated_id) REFERENCES webfront_entry (accession),
            CONSTRAINT fk_webfront_entry_webfront_database_source_database FOREIGN KEY (source_database) REFERENCES webfront_database (name)
        ) CHARSET=utf8 DEFAULT COLLATE utf8_unicode_ci
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
            mime_type VARCHAR(32) NOT NULL,
            CONSTRAINT fk_webfront_entryannotation_webfront_entry_accession_id FOREIGN KEY (accession_id) REFERENCES webfront_entry (accession)
        ) CHARSET=utf8 DEFAULT COLLATE utf8_unicode_ci
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
            left_number INT(11) NOT NULL,
            right_number INT(11) NOT NULL,
            CONSTRAINT fk_webfront_taxonomy_webfront_taxonomy_parent_id FOREIGN KEY (parent_id) REFERENCES webfront_taxonomy (accession)
        ) CHARSET=utf8 DEFAULT COLLATE utf8_unicode_ci
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
            proteomes LONGTEXT NOT NULL,
            gene VARCHAR(70),
            go_terms LONGTEXT NOT NULL,
            evidence_code INT(11) NOT NULL,
            source_database VARCHAR(10) NOT NULL,
            residues LONGTEXT NOT NULL,
            is_fragment TINYINT NOT NULL,
            structure LONGTEXT NOT NULL,
            tax_id VARCHAR(20) NOT NULL,
            extra_features LONGTEXT NOT NULL,
            CONSTRAINT fk_webfront_protein_webfront_taxonomy_tax_id FOREIGN KEY (tax_id) REFERENCES webfront_taxonomy (accession),
            CONSTRAINT fk_webfront_protein_webfront_database_source_database FOREIGN KEY (source_database) REFERENCES webfront_database (name)
        ) CHARSET=utf8 DEFAULT COLLATE utf8_unicode_ci
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
            CONSTRAINT fk_webfront_proteome_webfront_taxonomy_taxonomy_id FOREIGN KEY (taxonomy_id) REFERENCES webfront_taxonomy (accession)
        ) CHARSET=utf8 DEFAULT COLLATE utf8_unicode_ci
        """
    )

    cur.execute(
        """
        CREATE TABLE webfront_structure
        (
            accession VARCHAR(4) PRIMARY KEY NOT NULL,
            name VARCHAR(512) NOT NULL,
            short_name VARCHAR(32),
            other_names LONGTEXT NOT NULL,
            source_database VARCHAR(10) NOT NULL,
            experiment_type VARCHAR(16) NOT NULL,
            release_date DATETIME NOT NULL,
            resolution FLOAT,
            chains LONGTEXT NOT NULL,
            literature LONGTEXT NOT NULL,
            CONSTRAINT fk_webfront_structure_webfront_database_source_database FOREIGN KEY (source_database) REFERENCES webfront_database (name)
        ) CHARSET=utf8 DEFAULT COLLATE utf8_unicode_ci
        """
    )

    cur.execute(
        """
        CREATE TABLE webfront_set
        (
            accession VARCHAR(20) PRIMARY KEY NOT NULL,
            name VARCHAR(512) NOT NULL,
            description LONGTEXT NOT NULL,
            integrated LONGTEXT NOT NULL,
            relationships LONGTEXT NOT NULL,
            source_database VARCHAR(10) NOT NULL,
            is_set TINYINT NOT NULL,
            CONSTRAINT fk_webfront_set_webfront_database_source_database FOREIGN KEY (source_database) REFERENCES webfront_database (name)
        ) CHARSET=utf8 DEFAULT COLLATE utf8_unicode_ci
        """
    )

    cur.close()
    con.close()


def insert_taxa(ora_uri, my_uri, chunk_size=100000):
    taxa = interpro.get_taxa(ora_uri)

    con, cur = dbms.connect(my_uri)

    data = [(
        taxon['id'],
        taxon['sci_name'],
        taxon['full_name'],
        ' {} '.format(' '.join(taxon['lineage'])),  # leading/trailing whitespaces are important from API queries
        taxon['parent_id'],
        taxon['rank'],
        json.dumps(taxon['children']),
        taxon['left_number'],
        taxon['right_number']
    ) for taxon in taxa]

    for i in range(0, len(data), chunk_size):
        cur.executemany(
            """
            INSERT INTO webfront_taxonomy (
                accession,
                scientific_name,
                full_name,
                lineage,
                parent_id,
                rank,
                children,
                left_number,
                right_number
            ) VALUES (
              %s, %s, %s, %s, %s, %s, %s, %s, %s
            )
            """,
            data[i:i+chunk_size]
        )

    cur.close()
    con.commit()
    con.close()


def insert_proteomes(ora_uri, my_uri, chunk_size=100000):
    proteomes = uniprot.get_proteomes(ora_uri)

    con, cur = dbms.connect(my_uri)

    cur.execute('SELECT accession FROM webfront_taxonomy')
    taxa = set([row[0] for row in cur])

    data = []
    for p in proteomes:
        if p['tax_id'] not in taxa:
            # if tax_id not in taxa, it's very likely that INTERPRO.ETAXI (source for taxonomy table) is out-of-date
            logging.warning('missing taxon (ID: {})'.format(p['tax_id']))
            continue

        data.append((
            p['accession'],
            p['name'],
            1 if p['is_reference'] else 0,
            p['strain'],
            p['assembly'],
            p['tax_id']
        ))

    for i in range(0, len(data), chunk_size):
        cur.executemany(
            """
            INSERT INTO webfront_proteome (
              accession,
              name,
              is_reference,
              strain,
              assembly,
              taxonomy_id
            ) VALUES (
              %s, %s, %s, %s, %s, %s
            )
            """,
            data[i:i+chunk_size]
        )

    cur.close()
    con.commit()
    con.close()


def insert_databases(ora_uri, my_uri):
    databases = interpro.get_databases(ora_uri)

    con, cur = dbms.connect(my_uri)
    cur.executemany(
        """
        INSERT INTO webfront_database (
          name, name_long, description, version, release_date, type, prev_version, prev_release_date
        )
        VALUES (%s, %s, %s, %s, %s, %s, %s, %s)
        """,
        databases
    )
    cur.close()
    con.commit()
    con.close()


def insert_entries(ora_uri, pfam_uri, my_uri, chunk_size=100000):
    entries = interpro.get_entries(ora_uri)
    wiki = interpro.get_pfam_wiki(pfam_uri)

    data = [(
        None,  # entry_id
        e['accession'],
        e['type'],
        e['name'],
        e['short_name'],
        e['database'],
        json.dumps(e['member_databases']),
        e['integrated'],
        json.dumps(e['go_terms']),
        json.dumps(e['descriptions']),
        json.dumps(wiki.get(e['accession'], {})),
        json.dumps(e['citations']),
        json.dumps(e['hierarchy']),
        json.dumps(e['cross_references']),
        e['date'],
        json.dumps([]),  # overlapping entries
                         # requires supermatches so will be updated later (while populating Elastic)
        0,  # is_featured
    ) for e in entries]

    con, cur = dbms.connect(my_uri)

    for i in range(0, len(data), chunk_size):
        cur.executemany(
            """
              INSERT INTO webfront_entry (
                entry_id,
                accession,
                type,
                name,
                short_name,
                source_database,
                member_databases,
                integrated_id,
                go_terms,
                description,
                wikipedia,
                literature,
                hierarchy,
                cross_references,
                entry_date,
                overlaps_with,
                is_featured
              )
              VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
            """,
            data[i:i+chunk_size]
        )

    cur.close()
    con.commit()
    con.close()


def insert_annotations(pfam_uri, uri, chunk_size=10000):
    con, cur = dbms.connect(pfam_uri, sscursor=True)
    data = interpro.get_pfam_annotations(cur)
    cur.close()
    con.close()

    con, cur = dbms.connect(uri)
    for i in range(0, len(data), chunk_size):
        cur.executemany(
            """
            INSERT INTO webfront_entryannotation (
                annotation_id,
                accession_id,
                type,
                value,
                mime_type
            ) VALUES (%s, %s, %s, %s, %s)
            """,
            data[i:i+chunk_size]
        )

    cur.close()
    con.commit()
    con.close()


def insert_structures(ora_uri, uri, chunk_size=100000):
    structures = pdbe.get_structures(ora_uri, citations=True, fragments=False, by_protein=False)

    data = [(
        s['accession'],
        'pdb',
        s['date'],
        s['evidence'],
        s['name'],
        None,           # short_name
        s['resolution'],
        json.dumps(sorted([chain for chains in s['proteins'].values() for chain in chains])),
        json.dumps(s['citations']),
        json.dumps([])  # other_names
    ) for s in structures]

    con, cur = dbms.connect(uri)

    for i in range(0, len(data), chunk_size):
        cur.executemany(
            """
            INSERT INTO webfront_structure (
              accession,
              source_database,
              release_date,
              experiment_type,
              name,
              short_name,
              resolution,
              chains,
              literature,
              other_names
            ) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
            """,
            data[i:i + chunk_size]
        )

    cur.close()
    con.commit()
    con.close()


def insert_sets(pfam_uri, uri, chunk_size=100000):
    con, cur = dbms.connect(pfam_uri, sscursor=True)

    data = [(
        s['accession'].lower(),
        s['name'],
        s['description'],
        'pfam',
        json.dumps([]),
        json.dumps(s['relationships']),
        0
    ) for s in interpro.get_pfam_clans(cur)]

    cur.close()
    con.close()

    for s in interpro.get_cdd_superfamilies():
        data.append((
            s['accession'].lower(),
            s['name'],
            s['description'],
            'cdd',
            json.dumps([]),
            json.dumps(s['relationships']),
            0
        ))

    con, cur = dbms.connect(uri)
    for i in range(0, len(data), chunk_size):
        cur.executemany(
            """
            INSERT INTO webfront_set (
              accession,
              name,
              description,
              source_database,
              integrated,
              relationships,
              is_set
            ) VALUES (
              %s, %s, %s, %s, %s, %s, %s
            )
            """,
            data[i:i + chunk_size]
        )

    cur.close()
    con.commit()
    con.close()


def get_taxa(uri, slim=False):
    con, cur = dbms.connect(uri)
    taxa = {}

    if slim:
        sql = 'SELECT accession, scientific_name, full_name FROM webfront_taxonomy'
        cols = ('taxId', 'scientificName', 'fullName')
    else:
        sql = 'SELECT accession, scientific_name, full_name, lineage, rank FROM webfront_taxonomy'
        cols = ('taxId', 'scientificName', 'fullName', 'lineage', 'rank')

    cur.execute(sql)
    for row in cur:
        o = dict(zip(cols, row))
        taxa[o['taxId']] = o

    cur.close()
    con.close()

    return taxa


def get_proteomes(uri):
    con, cur = dbms.connect(uri)

    cur.execute(
        """
        SELECT accession, name, is_reference, strain, assembly
        FROM webfront_proteome
        """
    )

    proteomes = {}
    for row in cur:
        proteomes[row[0]] = {
            'name': row[1],
            'is_reference': bool(row[2]),
            'strain': row[3],
            'assembly': row[4]
        }

    cur.close()
    con.close()

    return proteomes


def insert_proteins(uri, proteins_f, evidences_f, descriptions_f, comments_f, proteomes_f, genes_f, annotations_f,
                    residues_f, struct_matches_f, prot_matches_extra_f, chunk_size=100000, limit=0):
    # MySQL data
    logging.info('loading taxa from MySQL')
    taxa = get_taxa(uri, slim=True)

    con, cur = dbms.connect(uri)
    logging.info('truncating table')
    cur.execute('TRUNCATE TABLE webfront_protein')
    logging.info('dropping indexes')
    for index in ('ui_webfront_protein_identifier', 'i_webfront_protein_length'):
        try:
            cur.execute('DROP INDEX {} ON webfront_protein'.format(index))
        except Exception:
            pass

    proteins = disk.Store(proteins_f)
    evidences = disk.Store(evidences_f)
    descriptions = disk.Store(descriptions_f)
    comments = disk.Store(comments_f)
    proteomes = disk.Store(proteomes_f)
    genes = disk.Store(genes_f)
    annotations = disk.Store(annotations_f)
    residues = disk.Store(residues_f)
    struct_matches = disk.Store(struct_matches_f)
    prot_matches_extra = disk.Store(prot_matches_extra_f)

    logging.info('inserting proteins')
    data = []
    cnt = 0
    ts = time.time()

    for acc, protein in proteins.iter():
        taxon_id = protein['taxon']
        taxon = taxa.get(taxon_id)

        if not taxon:
            logging.warning('invalid taxon ({}) for protein {}'.format(protein['taxon'], acc))
            continue

        evidence = evidences.get(acc)
        if not evidence:
            logging.warning('missing evidence for protein {}'.format(acc))
            continue

        if protein['length'] <= 100:
            size = 'small'
        elif protein['length'] <= 1000:
            size = 'medium'
        else:
            size = 'large'

        name, other_names = descriptions.get(acc, (None, None))

        # Enqueue record for protein table
        data.append((
            acc.lower(),
            protein['identifier'],
            json.dumps(taxon),
            name,
            json.dumps(other_names),
            json.dumps(comments.get(acc, [])),
            protein['sequence'],
            protein['length'],
            size,
            json.dumps(proteomes.get(acc, [])),
            genes.get(acc),
            json.dumps(annotations.get(acc, [])),
            evidence,
            'reviewed' if protein['isReviewed'] else 'unreviewed',
            json.dumps(residues.get(acc, {})),
            1 if protein['isFrag'] else 0,
            json.dumps(struct_matches.get(acc, {})),
            taxon_id,
            json.dumps(prot_matches_extra.get(acc, {}))
        ))

        if len(data) == chunk_size:
            cur.executemany(
                """
                INSERT INTO webfront_protein (
                  accession, identifier, organism, name, other_names, description, sequence, length, size,
                  proteomes, gene, go_terms, evidence_code, source_database, residues, is_fragment, structure, tax_id, 
                  extra_features
                )
                VALUES (
                  %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s
                )
                """,
                data
            )
            con.commit()
            data = []

        cnt += 1

        if not cnt % 1000000:
            logging.info('{:>12} ({:.0f} proteins/sec)'.format(cnt, cnt // (time.time() - ts)))

        if cnt == limit:
            break

    if data:
        cur.executemany(
            """
            INSERT INTO webfront_protein (
              accession, identifier, organism, name, other_names, description, sequence, length, size,
              proteomes, gene, go_terms, evidence_code, source_database, residues, is_fragment, structure, tax_id, 
              extra_features
            )
            VALUES (
              %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s
            )
            """,
            data
        )
        con.commit()
        data = []

    logging.info('{:>12} ({:.0f} proteins/sec)'.format(cnt, cnt // (time.time() - ts)))

    logging.info('indexing/analyzing table')
    cur = con.cursor()
    cur.execute('CREATE UNIQUE INDEX ui_webfront_protein_identifier ON webfront_protein (identifier)')
    cur.execute('CREATE INDEX i_webfront_protein_length ON webfront_protein (length)')
    cur.execute('ANALYZE TABLE webfront_protein')
    cur.close()
    con.close()

    logging.info('complete')


def _find_node(node, accession, relations=[]):
    """
    Find a entry (node) in its hierarchy tree and store its relations (ancestors + direct children)
    """
    if node['accession'] == accession:
        relations += [child['accession'] for child in node['children']]
        return node

    for child in node['children']:
        child = _find_node(child, accession, relations)

        if child:
            relations.append(node['accession'])
            return child

    return None


def get_entries(uri):
    con, cur = dbms.connect(uri)
    cur.execute(
        """
        SELECT
            accession, source_database, entry_date, description, integrated_id, name, type,
            short_name, member_databases, go_terms, literature, cross_references, hierarchy
        FROM webfront_entry
        """
    )

    entries = {}
    for row in cur:
        accession = row[0]
        relations = []
        hierarchy = json.loads(row[12])
        if hierarchy:
            _find_node(hierarchy, accession, relations)

        entries[accession] = {
            'accession': accession,
            'database': row[1],
            'date': row[2],
            'descriptions': json.loads(row[3]),
            'integrated': row[4],
            'name': row[5],
            'type': row[6],
            'short_name': row[7],
            'member_databases': json.loads(row[8]),
            'go_terms': json.loads(row[9]),
            'citations': json.loads(row[10]),
            'cross_references': json.loads(row[11]),
            'root': hierarchy.get('accession'),
            'relations': relations
        }

    cur.close()
    con.close()

    return entries


def get_sets(uri):
    con, cur = dbms.connect(uri)
    cur.execute(
        """
        SELECT accession, relationships, source_database
        FROM webfront_set
        """
    )

    sets = {}
    for row in cur:
        set_ac = row[0]
        rel = json.loads(row[1])
        database = row[2]

        for l in rel['links']:
            for k in ('source', 'target'):
                method_ac = l[k]

                if method_ac not in sets:
                    sets[method_ac] = {}

                sets[method_ac][set_ac] = database

    cur.close()
    con.close()

    return sets


def get_entry_databases(uri):
    con, cur = dbms.connect(uri)

    cur.execute("SELECT name, name_long FROM webfront_database WHERE type = 'entry'")
    databases = dict(cur.fetchall())

    cur.close()
    con.close()

    return databases
