#!/usr/bin/env python
# -*- coding: utf-8 -*-

import json
import logging
import time
from datetime import datetime

from . import dbms, disk
from .ebi import interpro, pdbe, uniprot


logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s: %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)


def init_tables(uri):
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
            name_long VARCHAR(25) NOT NULL,
            description LONGTEXT,
            version VARCHAR(20),
            release_date DATETIME,
            type ENUM('protein', 'entry', 'other') NOT NULL,
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
            member_databases LONGTEXT DEFAULT NULL,
            integrated_id VARCHAR(25) DEFAULT NULL,
            go_terms LONGTEXT DEFAULT NULL,
            description LONGTEXT DEFAULT NULL,
            wikipedia LONGTEXT DEFAULT NULL,
            literature LONGTEXT DEFAULT NULL,
            hierarchy LONGTEXT DEFAULT NULL,
            cross_references LONGTEXT DEFAULT NULL,
            entry_date DATETIME NOT NULL,
            overlaps_with LONGTEXT DEFAULT NULL,
            is_featured TINYINT NOT NULL DEFAULT 0,
            is_alive TINYINT NOT NULL DEFAULT 1,
            deletion_date DATETIME DEFAULT NULL,
            counts LONGTEXT DEFAULT NULL,
            CONSTRAINT fk_entry_entry
              FOREIGN KEY (integrated_id)
              REFERENCES webfront_entry (accession),
            CONSTRAINT fk_entry_database
              FOREIGN KEY (source_database)
              REFERENCES webfront_database (name)
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
            mime_type VARCHAR(32) NOT NULL,
            CONSTRAINT fk_entryannotation_entry
              FOREIGN KEY (accession_id)
              REFERENCES webfront_entry (accession)
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
            left_number INT(11) NOT NULL,
            right_number INT(11) NOT NULL,
            counts LONGTEXT DEFAULT NULL,
            CONSTRAINT fk_taxonomy_taxonomy
              FOREIGN KEY (parent_id)
              REFERENCES webfront_taxonomy (accession)
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
            counts LONGTEXT DEFAULT NULL,
            CONSTRAINT fk_protein_taxonomy
              FOREIGN KEY (tax_id)
              REFERENCES webfront_taxonomy (accession),
            CONSTRAINT fk_protein_database
              FOREIGN KEY (source_database)
              REFERENCES webfront_database (name)
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
            counts LONGTEXT DEFAULT NULL,
            CONSTRAINT fk_proteome_taxonomy
              FOREIGN KEY (taxonomy_id)
              REFERENCES webfront_taxonomy (accession)
        ) CHARSET=utf8 DEFAULT COLLATE=utf8_unicode_ci
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
            literature LONGTEXT NOT NULL,
            chains LONGTEXT NOT NULL,
            proteins LONGTEXT NOT NULL,
            counts LONGTEXT DEFAULT NULL,
            CONSTRAINT fk_structure_database
              FOREIGN KEY (source_database)
              REFERENCES webfront_database (name)
        ) CHARSET=utf8 DEFAULT COLLATE=utf8_unicode_ci
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
            counts LONGTEXT DEFAULT NULL,
            CONSTRAINT fk_set_database
              FOREIGN KEY (source_database)
              REFERENCES webfront_database (name)
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


def insert_taxa(ora_uri, my_uri, chunk_size=100000):
    taxa = interpro.get_taxa(ora_uri)

    con, cur = dbms.connect(my_uri)

    data = [(
        taxon['id'],
        taxon['sci_name'],
        taxon['full_name'],
        # leading/trailing whitespaces are important from API queries
        ' {} '.format(' '.join(taxon['lineage'])),
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
        # overlapping entries (requires supermatches: updated later)
        json.dumps([]),

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
                overlaps_with
              )
              VALUES (
                %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s
              )
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
    structures = pdbe.get_structures(ora_uri)

    data = []
    for s in structures.values():
        proteins = {}
        chains = set()

        for acc, p in s["proteins"].items():
            proteins[acc] = []
            for chain in p:
                chains.add(chain)
                proteins[acc].append(chain)

        data.append((
            s["id"],
            "pdb",
            s["data"],
            s["evidence"],
            s["name"],
            None,  # short name
            s["resolution"],
            json.dumps(sorted(chains)),
            json.dumps(s["citation"]),
            json.dumps([]),  # other names
            json.dumps(proteins)
        ))

    con, cur = dbms.connect(uri)

    for i in range(0, len(data), chunk_size):
        cur.executemany(
            """
            INSERT INTO webfront_structure (
              accession,
              name,
              short_name,
              other_names,
              source_database,
              experiment_type,
              release_date,
              resolution,
              literature,
              chains,
              proteins
            ) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
            """,
            data[i:i + chunk_size]
        )

    cur.close()
    con.commit()
    con.close()


def get_structures(uri: str) -> dict:
    structures = {}
    con, cur = dbms.connect(uri)
    cur.execute(
        """
        SELECT accession, name, experiment_type, resolution, proteins
        FROM webfront_structure
        """
    )

    for acc, name, _type, resolution, proteins in cur:
        structures[acc] = {
            "accession": acc,
            "name": name,
            "type": _type,
            "resolution": resolution,
            "proteins": json.loads(proteins)
        }

    cur.close()
    con.close()

    return structures


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


def insert_proteins(uri, proteins_f, sequences_f, evidences_f,
                    descriptions_f, comments_f, proteomes_f,
                    genes_f, annotations_f, residues_f,
                    struct_matches_f, prot_matches_extra_f,
                    chunk_size=100000, limit=0):
    logging.info('starting')

    # MySQL data
    taxa = get_taxa(uri, slim=True)

    proteins = disk.Store(proteins_f)
    sequences = disk.Store(sequences_f)
    evidences = disk.Store(evidences_f)
    descriptions = disk.Store(descriptions_f)
    comments = disk.Store(comments_f)
    proteomes = disk.Store(proteomes_f)
    genes = disk.Store(genes_f)
    annotations = disk.Store(annotations_f)
    residues = disk.Store(residues_f)
    struct_matches = disk.Store(struct_matches_f)
    prot_matches_extra = disk.Store(prot_matches_extra_f)

    con, cur = dbms.connect(uri)
    cur.execute('TRUNCATE TABLE webfront_protein')
    for index in ("ui_webfront_protein_identifier",
                  "i_webfront_protein_length"):
        try:
            cur.execute('DROP INDEX {} ON webfront_protein'.format(index))
        except Exception:
            pass

    logging.info('inserting proteins')
    data = []
    n_proteins = 0
    unknown_taxa = {}
    ts = time.time()
    for acc, protein in proteins.iter():
        tax_id = protein['taxon']

        try:
            taxon = taxa[tax_id]
        except KeyError:
            if tax_id in unknown_taxa:
                unknown_taxa[tax_id] += 1
            else:
                unknown_taxa[tax_id] = 1
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
            sequences.get(acc),
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
            tax_id,
            json.dumps(prot_matches_extra.get(acc, {}))
        ))

        if len(data) == chunk_size:
            cur.executemany(
                """
                INSERT INTO webfront_protein (
                  accession, identifier, organism, name, other_names,
                  description, sequence, length, size,
                  proteomes, gene, go_terms, evidence_code, source_database,
                  residues, is_fragment, structure, tax_id,
                  extra_features
                )
                VALUES (
                  %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s,
                  %s, %s, %s, %s, %s, %s, %s, %s
                )
                """,
                data
            )
            con.commit()
            data = []

        n_proteins += 1
        if n_proteins == limit:
            break
        elif not n_proteins % 1000000:
            logging.info('{:>12} ({:.0f} proteins/sec)'.format(
                n_proteins, n_proteins // (time.time() - ts)
            ))

    if data:
        cur.executemany(
            """
            INSERT INTO webfront_protein (
              accession, identifier, organism, name, other_names,
              description, sequence, length, size,
              proteomes, gene, go_terms, evidence_code, source_database,
              residues, is_fragment, structure, tax_id,
              extra_features
            )
            VALUES (
              %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s,
              %s, %s, %s, %s
            )
            """,
            data
        )
        con.commit()
        data = []

    logging.info('{:>12} ({:.0f} proteins/sec)'.format(
        n_proteins, n_proteins // (time.time() - ts)
    ))

    if unknown_taxa:
        logging.warning("{} unknown taxa:".format(len(unknown_taxa)))
        for tax_id in sorted(unknown_taxa):
            logging.warning("\t{:>8}\t{:>12} skipped proteins".format(
                tax_id, unknown_taxa[tax_id]
            ))

    logging.info('indexing/analyzing table')
    cur = con.cursor()
    cur.execute(
        """
        CREATE UNIQUE INDEX ui_webfront_protein_identifier
        ON webfront_protein (identifier)
        """
    )
    cur.execute(
        """
        CREATE INDEX i_webfront_protein_length
        ON webfront_protein (length)
        """
    )
    cur.execute("ANALYZE TABLE webfront_protein")
    cur.close()
    con.close()

    logging.info("complete")


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


def get_sets(uri, by_members=True):
    con, cur = dbms.connect(uri)
    cur.execute(
        """
        SELECT accession, source_database, relationships
        FROM webfront_set
        """
    )

    sets = {}
    for row in cur:
        set_ac = row[0]
        database = row[1]

        if by_members:
            rel = json.loads(row[2])

            for l in rel['links']:
                for k in ('source', 'target'):
                    method_ac = l[k]

                    # todo: can a method belong to more than one set?
                    if method_ac in sets:
                        sets[method_ac][set_ac] = database
                    else:
                        set_ac[method_ac] = {set_ac: database}
        else:
            sets[set_ac] = database

    cur.close()
    con.close()

    return sets


def get_entry_databases(uri):
    con, cur = dbms.connect(uri)

    cur.execute("SELECT name, name_long, version FROM webfront_database WHERE type = 'entry'")
    databases = {}
    for name, name_long, version in cur:
        databases[name] = {
            "name_long": name_long,
            "version": version
        }

    cur.close()
    con.close()

    return databases


def make_release_notes(stg_uri, rel_uri, proteins_f, prot_matches_f,
                       struct_matches_f, proteomes_f, version, release_date):
    con, cur = dbms.connect(stg_uri)

    # Get PDB structures
    cur.execute(
        """
        SELECT accession, release_date
        FROM webfront_structure
        ORDER BY release_date
        """
    )
    structures = set()
    pdbe_release_date = None
    for row in cur:
        structures.add(row[0])
        pdbe_release_date = row[1]

    # Get proteomes
    proteomes = set(get_proteomes(stg_uri))

    # Get taxa
    taxa = set(get_taxa(stg_uri, slim=True))

    # Get UniProtKB version
    cur.execute(
        """
        SELECT name_long, version
        FROM webfront_database
        WHERE type='protein'
        """
    )

    proteins = {
        name: {
            "version": version,
            "count": 0,
            "signatures": 0,
            "integrated_signatures": 0
        }
        for name, version in cur
    }

    # Get sets
    db2set = {}
    for set_ac, db in get_sets(stg_uri, by_members=False).items():
        if db in db2set:
            db2set[db] += 1
        else:
            db2set[db] = 1

    cur.close()
    con.close()

    proteins_s = disk.Store(proteins_f)
    prot_matches_s = disk.Store(prot_matches_f)
    struct_matches_s = disk.Store(struct_matches_f)
    proteomes_s = disk.Store(proteomes_f)

    interpro_structures = set()
    interpro_proteomes = set()
    interpro_taxa = set()

    n_proteins = 0
    ts = time.time()
    for acc, protein in proteins_s.iter():
        if protein["isReviewed"]:
            k = "UniProtKB/Swiss-Prot"
        else:
            k = "UniProtKB/TrEMBL"

        proteins[k]["count"] += 1
        matches = prot_matches_s.get(acc)
        if matches:
            for m in matches:
                if m["entry_ac"]:
                    proteins[k]["integrated_signatures"] += 1
                    interpro_taxa.add(protein["taxon"])
                    interpro_proteomes |= set(proteomes_s.get(acc, []))
                    _structures = struct_matches_s.get(acc)
                    if _structures:
                        interpro_structures |= {
                            v["domain_id"]
                            for v in
                            _structures["feature"].get("pdb", {}).values()
                        }
                    break

            proteins[k]["signatures"] += 1

        n_proteins += 1
        if not n_proteins % 1000000:
            logging.info("{:>12} ({:.0f} proteins/sec)".format(
                n_proteins, n_proteins // (time.time() - ts)
            ))

    proteins_s.close()
    prot_matches_s.close()
    struct_matches_s.close()
    proteomes_s.close()

    for k in ("count", "signatures", "integrated_signatures"):
        proteins["UniProtKB"][k] = (proteins["UniProtKB/Swiss-Prot"][k]
                                    + proteins["UniProtKB/TrEMBL"][k])

    logging.info("{:>12} ({:.0f} proteins/sec)".format(
        n_proteins, n_proteins // (time.time() - ts)
    ))

    bad = interpro_structures - structures
    if bad:
        logging.warning("structures issues: {}".format(bad))

    bad = interpro_proteomes - proteomes
    if bad:
        logging.warning("proteomes issues: {}".format(bad))

    bad = interpro_taxa - taxa
    if bad:
        logging.warning("taxonomy issues: {}".format(bad))

    notes = {
        "interpro": {},
        "member_databases": [],
        "proteins": proteins,
        "structures": {
            "total": len(structures),
            "integrated": len(interpro_structures & structures),
            "version": pdbe_release_date.strftime("%Y-%m-%d")
        },
        "proteomes": {
            "total": len(proteomes),
            "integrated": len(interpro_proteomes & proteomes),
            "version": proteins["UniProtKB"]["version"]
        },
        "taxonomy": {
            "total": len(taxa),
            "integrated": len(interpro_taxa & taxa),
            "version": proteins["UniProtKB"]["version"]
        }
    }

    rel_entries = []
    rel_interpro_entries = set()
    already_integrated = set()
    for e in get_entries(rel_uri).values():
        rel_entries.append(e)
        if e["database"] == "interpro":
            rel_interpro_entries.add(e["accession"])
        elif e["integrated"]:
            # Signature already integrated during the last release
            already_integrated.add(e["accession"])

    stg_entries = []
    new_entries = []
    for e in get_entries(stg_uri).values():
        stg_entries.append(e)
        acc = e["accession"]
        if e["database"] == "interpro" and acc not in rel_interpro_entries:
            new_entries.append(acc)

    # Member database changes
    stg_dbs = get_entry_databases(stg_uri)
    rel_dbs = get_entry_databases(rel_uri)
    updated_databases = set()
    new_databases = set()
    for database in stg_dbs:
        if database not in rel_dbs:
            new_databases.add(database)
        elif stg_dbs[database]["version"] != rel_dbs[database]["version"]:
            updated_databases.add(database)

    member_databases = {}
    interpro_types = {}
    citations = set()
    n_interpro2go = 0
    latest_entry = None
    for entry in sorted(stg_entries, key=lambda x: x["accession"]):
        acc = entry["accession"]
        database = entry["database"]
        _type = entry["type"]

        citations |= {
            item["PMID"]
            for item in entry["citations"].values()
            if item["PMID"] is not None
        }

        if database == "interpro":
            if _type in interpro_types:
                interpro_types[_type] += 1
            else:
                interpro_types[_type] = 1

            n_interpro2go += len(entry["go_terms"])

            latest_entry = acc
        else:
            if database in member_databases:
                db = member_databases[database]
            else:
                db = member_databases[database] = {
                    "name": stg_dbs[database]["name_long"],
                    "version": stg_dbs[database]["version"],
                    "signatures": 0,
                    "integrated_signatures": 0,
                    "recently_integrated": [],
                    "is_new": database in new_databases,
                    "is_updated": database in updated_databases,
                    "sets": db2set.get(database, 0)
                }

            db["signatures"] += 1
            if entry["integrated"]:
                db["integrated_signatures"] += 1

                if acc not in already_integrated:
                    # Recent integration
                    db["recently_integrated"].append(acc)

    notes.update({
        "interpro": {
            "entries": sum(interpro_types.values()),
            "new_entries": new_entries,
            "latest_entry": latest_entry,
            "types": interpro_types,
            "go_terms": n_interpro2go
        },
        "member_databases": member_databases,
        "citations": len(citations)
    })

    con, cur = dbms.connect(stg_uri)
    cur.execute(
        """
        SELECT COUNT(*)
        FROM webfront_release_note
        WHERE version = %s
        """,
        (version,)
    )
    n = cur.fetchone()[0]

    if n:
        cur.execute(
            """
            UPDATE webfront_release_note
            SET content = %s
            WHERE version = %s
            """,
            (json.dumps(notes), version)
        )
    else:
        cur.execute(
            """
            INSERT INTO webfront_release_note (
              version, release_date, content
            )
            VALUES (%s, %s, %s)
            """,
            (
                version,
                datetime.strptime(release_date, '%Y-%m-%d'),
                json.dumps(notes)
            )
        )

    cur.close()
    con.commit()
    con.close()
    logging.info("complete")
