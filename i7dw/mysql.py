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
            member_databases LONGTEXT,
            integrated_id VARCHAR(25),
            go_terms LONGTEXT,
            description LONGTEXT,
            wikipedia LONGTEXT,
            literature LONGTEXT,
            hierarchy LONGTEXT,
            cross_references LONGTEXT,
            entry_date DATETIME NOT NULL,
            overlaps_with LONGTEXT DEFAULT NULL,
            is_featured TINYINT NOT NULL DEFAULT 0,
            is_alive TINYINT NOT NULL DEFAULT 1,
            deletion_date DATETIME,
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
            counts LONGTEXT DEFAULT NULL,
            CONSTRAINT fk_protein_taxonomy
              FOREIGN KEY (tax_id)
              REFERENCES webfront_taxonomy (accession),
            CONSTRAINT fk_protein_database
              FOREIGN KEY (source_database)
              REFERENCES webfront_database (name),
            CONSTRAINT fk_protein_proteome
              FOREIGN KEY (proteome)
              REFERENCES webfront_proteome (accession)
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
            name VARCHAR(100),
            description TEXT,
            source_database VARCHAR(10) NOT NULL,
            is_set TINYINT NOT NULL,
            relationships LONGTEXT NOT NULL,
            integrated LONGTEXT DEFAULT NULL,
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
    taxa = get_taxa(my_uri, "basic")

    data = []
    con, cur = dbms.connect(my_uri)
    for p in proteomes.values():
        if p["tax_id"] not in taxa:
            """
            If tax_id not in taxa, it's very likely that INTERPRO.ETAXI
            (source for taxonomy table) is out-of-date
            """
            logging.warning("missing taxon (ID: {})".format(p["tax_id"]))
            continue

        data.append((
            p["accession"],
            p["name"],
            1 if p["is_reference"] else 0,
            p["strain"],
            p["assembly"],
            p["tax_id"]
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
          name, name_long, description, version, release_date, type,
          prev_version, prev_release_date
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
        e["accession"],
        e["type"],
        e["name"],
        e["short_name"],
        e["database"],
        json.dumps(e["member_databases"]),
        e["integrated"],
        json.dumps(e["go_terms"]),
        json.dumps(e["descriptions"]),
        json.dumps(wiki.get(e["accession"])),
        json.dumps(e["citations"]),
        json.dumps(e["hierarchy"]),
        json.dumps(e["cross_references"]),
        e["date"],
        1,      # is alive
        None    # deletion date
    ) for e in entries]

    for e in interpro.get_deleted_entries(ora_uri):
        if e["deletion_date"] is None:
            e["deletion_date"] = e["creation_date"]

        data.append((
            e["accession"],
            e["type"],
            e["name"],
            e["short_name"],
            e["database"],
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            e["creation_date"],
            0,
            e["deletion_date"]
        ))

    con, cur = dbms.connect(my_uri)

    for i in range(0, len(data), chunk_size):
        cur.executemany(
            """
              INSERT INTO webfront_entry (
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
                is_alive,
                deletion_date
              )
              VALUES (
                %s, %s, %s, %s, %s, %s, %s, %s, %s, %s,
                %s, %s, %s, %s, %s, %s
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
            s["name"],
            "pdb",
            s["evidence"],
            s["date"],
            s["resolution"],
            json.dumps(s["citations"]),
            json.dumps(sorted(chains)),
            json.dumps(proteins)
        ))

    con, cur = dbms.connect(uri)

    for i in range(0, len(data), chunk_size):
        cur.executemany(
            """
            INSERT INTO webfront_structure (
              accession,
              name,
              source_database,
              experiment_type,
              release_date,
              resolution,
              literature,
              chains,
              proteins
            ) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s)
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


def insert_sets(ora_uri, pfam_uri, my_uri, chunk_size=100000):
    data = []

    con, cur = dbms.connect(pfam_uri, sscursor=True)
    sets = interpro.get_profile_alignments(ora_uri, "pfam")
    for clan in interpro.get_pfam_clans(cur):
        set_ac = clan["accession"].lower()

        try:
            rel = sets[set_ac]["relationships"]
        except KeyError:
            # TODO: warning/error
            continue
        else:
            # Use nodes from Pfam database (for the score...)
            rel["nodes"] = clan["relationships"]["nodes"]

            data.append((
                set_ac,
                clan["name"],
                clan["description"],
                "pfam",
                1,
                json.dumps(rel)
            ))

    cur.close()
    con.close()

    sets = interpro.get_profile_alignments(ora_uri, "cdd")
    for supfam in interpro.get_cdd_superfamilies():
        set_ac = supfam["accession"].lower()

        try:
            s = sets[set_ac]
        except KeyError:
            # TODO: warning/error
            continue
        else:
            data.append((
                set_ac,
                supfam["name"],
                supfam["description"],
                "cdd",
                1,
                json.dumps(s["relationships"])
            ))

    for s in interpro.get_profile_alignments(ora_uri, "panther").values():
        data.append((
            s["accession"],
            s["name"],          # None
            s["description"],   # None
            "panther",
            1,
            json.dumps(s["relationships"])
        ))

    for s in interpro.get_profile_alignments(ora_uri, "pirsf").values():
        data.append((
            s["accession"],
            s["name"],          # None
            s["description"],   # None
            "pirsf",
            1,
            json.dumps(s["relationships"])
        ))

    con, cur = dbms.connect(my_uri)
    for i in range(0, len(data), chunk_size):
        cur.executemany(
            """
            INSERT INTO webfront_set (
              accession, name, description, source_database, is_set,
              relationships
            ) VALUES (%s, %s, %s, %s, %s, %s)
            """,
            data[i:i + chunk_size]
        )

    cur.close()
    con.commit()
    con.close()


def get_taxa(uri: str, method: str=None):
    if method is None:
        method = "default"
    elif method not in ("basic", "default", "complete"):
        raise ValueError("cannot find context for {}".format(method))

    con, cur = dbms.connect(uri)
    if method == "basic":
        cur.execute("SELECT accession FROM webfront_taxonomy")
        taxa = {row[0] for row in cur}
    elif method == "default":
        cur.execute(
            """
            SELECT accession, scientific_name, full_name
            FROM webfront_taxonomy
            """
        )
        cols = ("taxId", "scientificName", "fullName")
        taxa = {row[0]: dict(zip(cols, row)) for row in cur}
    else:
        cur.execute(
            """
            SELECT accession, scientific_name, full_name, lineage, rank
            FROM webfront_taxonomy
            """
        )
        cols = ("taxId", "scientificName", "fullName", "lineage", "rank")
        taxa = {row[0]: dict(zip(cols, row)) for row in cur}

    cur.close()
    con.close()

    return taxa


def get_proteomes(uri: str) -> dict:
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


def insert_proteins(uri, src_proteins, src_sequences, src_misc,
                    src_names, src_comments, src_proteomes,
                    src_residues, src_structures, src_features,
                    src_matches, chunk_size=100000, limit=0):
    logging.info("starting")

    # MySQL data
    taxa = get_taxa(uri, method="default")
    entries = get_entries(uri)
    integrated = {
        acc: e["integrated"]
        for acc, e in entries.items()
        if e["integrated"]
    }

    protein2pdb = {}
    for pdb_id, s in get_structures(uri).items():
        for acc in s["proteins"]:
            if acc in protein2pdb:
                protein2pdb[acc] += 1
            else:
                protein2pdb[acc] = 1

    entry2set = {}
    for set_ac, s in get_sets(uri).items():
        for acc in s["members"]:
            entry2set[acc] = set_ac

    proteins = disk.Store(src_proteins)
    protein2sequence = disk.Store(src_sequences)
    protein2misc = disk.Store(src_misc)
    protein2names = disk.Store(src_names)
    protein2comments = disk.Store(src_comments)
    protein2proteome = disk.Store(src_proteomes)
    protein2residues = disk.Store(src_residues)
    protein2structures = disk.Store(src_structures)
    protein2features = disk.Store(src_features)
    protein2matches = disk.Store(src_matches)

    con, cur = dbms.connect(uri)
    cur.execute("TRUNCATE TABLE webfront_protein")
    for index in ("ui_webfront_protein_identifier",
                  "i_webfront_protein_length"):
        try:
            cur.execute("DROP INDEX {} ON webfront_protein".format(index))
        except Exception:
            pass

    logging.info("inserting proteins")
    data = []
    n_proteins = 0
    ts = time.time()
    for acc, protein in proteins:
        if protein["length"] <= 100:
            size = "small"
        elif protein["length"] <= 1000:
            size = "medium"
        else:
            size = "large"

        evidence, gene = protein2misc.get(acc, (None, None))
        name, other_names = protein2names.get(acc, (None, None))

        # InterPro2GO + InterPro matches -> UniProt-GO
        go_terms = {}
        _entries = set()
        for m in protein2matches.get(acc, []):
            method_ac = m["method_ac"]
            _entries.add(method_ac)

            if method_ac in integrated:
                entry_ac = integrated[method_ac]

                _entries.add(entry_ac)
                for term in entries[entry_ac]["go_terms"]:
                    go_terms[term["identifier"]] = term

        protein2entries = {"total": len(_entries)}
        protein2sets = set()
        for entry_ac in _entries:
            db_name = entries[entry_ac]["database"]
            if db_name in protein2entries:
                protein2entries[db_name] += 1
            else:
                protein2entries[db_name] = 1

            if entry_ac in entry2set:
                protein2sets.add(entry2set[entry_ac])

        # Enqueue record for protein table
        data.append((
            acc.lower(),
            protein["identifier"],
            json.dumps(taxa[protein["taxon"]]),
            name,
            json.dumps(other_names),
            json.dumps(protein2comments.get(acc, [])),
            protein2sequence[acc],
            protein["length"],
            size,
            protein2proteome.get(acc),
            gene,
            json.dumps(list(go_terms.values())),
            evidence,
            "reviewed" if protein["isReviewed"] else "unreviewed",
            json.dumps(protein2residues.get(acc, {})),
            1 if protein["isFrag"] else 0,
            json.dumps(protein2structures.get(acc, {})),
            protein["taxon"],
            json.dumps(protein2features.get(acc, {})),
            json.dumps({
                "entries": protein2entries,
                "structures": protein2pdb.get(acc, 0),
                "sets": len(protein2sets),
                "proteomes": 1 if upid else 0,
                "taxa": 1
            })
        ))

        if len(data) == chunk_size:
            cur.executemany(
                """
                INSERT INTO webfront_protein (
                  accession, identifier, organism, name, other_names,
                  description, sequence, length, size,
                  proteome, gene, go_terms, evidence_code, source_database,
                  residues, is_fragment, structure, tax_id,
                  extra_features, counts
                )
                VALUES (
                  %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s,
                  %s, %s, %s, %s, %s, %s, %s, %s, %s
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
            logging.info('{:>12,} ({:.0f} proteins/sec)'.format(
                n_proteins, n_proteins / (time.time() - ts)
            ))

    if data:
        cur.executemany(
            """
            INSERT INTO webfront_protein (
              accession, identifier, organism, name, other_names,
              description, sequence, length, size,
              proteome, gene, go_terms, evidence_code, source_database,
              residues, is_fragment, structure, tax_id,
              extra_features, counts
            )
            VALUES (
              %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s,
              %s, %s, %s, %s,  %s
            )
            """,
            data
        )
        con.commit()

    logging.info('{:>12,} ({:.0f} proteins/sec)'.format(
        n_proteins, n_proteins / (time.time() - ts)
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
    Find a entry (node) in its hierarchy tree
    and store its relations (ancestors + direct children)
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


def get_entries(uri: str) -> dict:
    con, cur = dbms.connect(uri)
    cur.execute(
        """
        SELECT
            accession, source_database, entry_date, description,
            integrated_id, name, type, short_name, member_databases,
            go_terms, literature, cross_references, hierarchy
        FROM webfront_entry
        WHERE is_alive = 1
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
            "accession": accession,
            "database": row[1],
            "date": row[2],
            "descriptions": json.loads(row[3]),
            "integrated": row[4],
            "name": row[5],
            "type": row[6],
            "short_name": row[7],
            "member_databases": json.loads(row[8]),
            "go_terms": json.loads(row[9]),
            "citations": json.loads(row[10]),
            "cross_references": json.loads(row[11]),
            "root": hierarchy.get("accession"),
            "relations": relations
        }

    cur.close()
    con.close()

    return entries


def get_sets(uri: str) -> dict:
    con, cur = dbms.connect(uri)
    cur.execute(
        """
        SELECT accession, source_database, relationships
        FROM webfront_set
        """
    )

    sets = {}
    for acc, database, relationships in cur:
        sets[acc] = {
            "database": database,
            "members": [
                n["accession"]
                for n in json.loads(relationships)["nodes"]
            ]
        }

    cur.close()
    con.close()

    return sets


def get_entry_databases(uri: str) -> dict:
    con, cur = dbms.connect(uri)

    cur.execute(
        """
        SELECT name, name_long, version
        FROM webfront_database WHERE type = 'entry'
        """
    )
    databases = {}
    for name, name_long, version in cur:
        databases[name] = {
            "name_long": name_long,
            "version": version
        }

    cur.close()
    con.close()

    return databases


def make_release_notes(stg_uri, rel_uri, src_proteins, src_matches,
                       src_structures, src_proteomes, version, release_date):
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
    taxa = get_taxa(stg_uri, method="basic")

    # Integrated signatures
    integrated = {
        acc
        for acc, e in get_entries(stg_uri).items()
        if e["integrated"]
    }

    # Get UniProtKB version
    cur.execute(
        """
        SELECT name_long, version
        FROM webfront_database
        WHERE type='protein'
        """
    )

    uniprot = {
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
    for set_ac, s in get_sets(stg_uri).items():
        db = s["database"]
        if db in db2set:
            db2set[db] += 1
        else:
            db2set[db] = 1

    cur.close()
    con.close()

    proteins = disk.Store(src_proteins)
    protein2matches = disk.Store(src_matches)
    protein2structures = disk.Store(src_structures)
    protein2proteome = disk.Store(src_proteomes)

    interpro_structures = set()
    interpro_proteomes = set()
    interpro_taxa = set()

    n_proteins = 0
    ts = time.time()
    for acc, protein in proteins:
        if protein["isReviewed"]:
            k = "UniProtKB/Swiss-Prot"
        else:
            k = "UniProtKB/TrEMBL"

        uniprot[k]["count"] += 1
        upid = protein2proteome.get(acc)
        matches = protein2matches.get(acc)
        if matches:
            # Protein has a list one signature
            uniprot[k]["signatures"] += 1

            # Search if the protein has at least one integrated signature
            for m in matches:
                if m["method_ac"] in integrated:
                    # It has!
                    uniprot[k]["integrated_signatures"] += 1

                    # Add taxon, proteome, and structures
                    interpro_taxa.add(protein["taxon"])

                    if upid:
                        interpro_proteomes.add(upid)

                    _structures = protein2structures.get(acc)
                    if _structures:
                        interpro_structures |= {
                            v["domain_id"]
                            for v in
                            _structures["feature"].get("pdb", {}).values()
                        }

                    break

        n_proteins += 1
        if not n_proteins % 1000000:
            logging.info("{:>12,} ({:.0f} proteins/sec)".format(
                n_proteins, n_proteins / (time.time() - ts)
            ))

    proteins.close()
    protein2matches.close()
    protein2structures.close()
    protein2proteome.close()

    for k in ("count", "signatures", "integrated_signatures"):
        uniprot["UniProtKB"][k] = (uniprot["UniProtKB/Swiss-Prot"][k]
                                    + uniprot["UniProtKB/TrEMBL"][k])

    logging.info("{:>12,} ({:.0f} proteins/sec)".format(
        n_proteins, n_proteins / (time.time() - ts)
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
        "proteins": uniprot,
        "structures": {
            "total": len(structures),
            "integrated": len(interpro_structures & structures),
            "version": pdbe_release_date.strftime("%Y-%m-%d")
        },
        "proteomes": {
            "total": len(proteomes),
            "integrated": len(interpro_proteomes & proteomes),
            "version": uniprot["UniProtKB"]["version"]
        },
        "taxonomy": {
            "total": len(taxa),
            "integrated": len(interpro_taxa & taxa),
            "version": uniprot["UniProtKB"]["version"]
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
            # Signature already integrated during the previous release
            already_integrated.add(e["accession"])

    stg_entries = []
    new_entries = []
    for acc, e in get_entries(stg_uri).items():
        stg_entries.append(e)
        if e["database"] == "interpro" and acc not in rel_interpro_entries:
            new_entries.append(acc)

    # Member database changes
    stg_dbs = get_entry_databases(stg_uri)
    rel_dbs = get_entry_databases(rel_uri)
    updated_databases = set()
    new_databases = set()
    for name, info in stg_dbs.items():
        if name not in rel_dbs:
            new_databases.add(name)
        elif info["version"] != rel_dbs[name]["version"]:
            updated_databases.add(name)

    member_databases = {}
    interpro_types = {}
    citations = set()
    n_interpro2go = 0
    latest_entry = None
    for entry in sorted(stg_entries, key=lambda x: x["accession"]):
        acc = entry["accession"]
        db_name = entry["database"]
        _type = entry["type"]

        citations |= {
            item["PMID"]
            for item in entry["citations"].values()
            if item["PMID"] is not None
        }

        if db_name == "interpro":
            if _type in interpro_types:
                interpro_types[_type] += 1
            else:
                interpro_types[_type] = 1

            n_interpro2go += len(entry["go_terms"])

            latest_entry = acc
        else:
            if db_name in member_databases:
                db = member_databases[db_name]
            else:
                db = member_databases[db_name] = {
                    "name": stg_dbs[db_name]["name_long"],
                    "version": stg_dbs[db_name]["version"],
                    "signatures": 0,
                    "integrated_signatures": 0,
                    "recently_integrated": [],
                    "is_new": db_name in new_databases,
                    "is_updated": db_name in updated_databases,
                    "sets": db2set.get(db_name, 0)
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
