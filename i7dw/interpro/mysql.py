#!/usr/bin/env python
# -*- coding: utf-8 -*-

import json
import logging
import time
from datetime import datetime

from . import oracle
from .. import cdd, dbms, io, pdbe, pfam, uniprot


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
    taxa = oracle.get_taxa(ora_uri)

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
    taxa = set(get_taxa(my_uri, lineage=False))

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


def insert_databases(ora_uri: str, my_uri: str, version: str,
                     release_date: str):
    data = []
    for db in oracle.get_databases(ora_uri):
        if db["name"] == "interpro" and db["version"]["code"] != version:
            """
            Happens when Oracle hasn't been updated yet
            (DB_VERSION still on the previous release)

            --> move the `version` to `previous_version`
            """
            db["previous_version"] = db["version"]
            db["version"] = {
                "code": version,
                "date": datetime.strptime(release_date, "%Y-%m-%d"),
            }

        data.append((
            db["name"],
            db["name_long"],
            db["description"],
            db["type"],
            db["version"]["code"],
            db["version"]["date"],
            db["previous_version"]["code"],
            db["previous_version"]["date"]
        ))

    con, cur = dbms.connect(my_uri)
    cur.executemany(
        """
        INSERT INTO webfront_database (
          name, name_long, description, type,
          version, release_date, prev_version, prev_release_date
        )
        VALUES (%s, %s, %s, %s, %s, %s, %s, %s)
        """,
        data
    )
    con.commit()
    cur.close()
    con.close()


def insert_entries(ora_uri, pfam_uri, my_uri, chunk_size=100000):
    entries = oracle.get_entries(ora_uri)
    wiki = pfam.get_wiki(pfam_uri)

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

    for e in oracle.get_deleted_entries(ora_uri):
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
    data = pfam.get_annotations(pfam_uri)

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


def insert_sets(ora_uri, pfam_uri, my_uri, chunk_size=100000):
    data = []

    logging.info("loading Pfam clans")
    clans = pfam.get_clans(pfam_uri)
    for s in oracle.get_profile_alignments(ora_uri, "pfam"):
        set_ac = s["accession"]

        try:
            clan = clans[set_ac]
        except KeyError:
            continue  # TODO: warning

        # Use nodes from Pfam DB for the score
        s["relationships"]["nodes"] = clan["relationships"]["nodes"]

        data.append((
            set_ac,
            clan["name"],
            clan["description"],
            "pfam",
            1,
            json.dumps(s["relationships"])
        ))

    logging.info("loading CDD superfamilies")
    supfams = cdd.get_superfamilies()
    for s in oracle.get_profile_alignments(ora_uri, "cdd"):
        set_ac = s["accession"]

        try:
            supfam = supfams[set_ac]
        except KeyError:
            continue  # TODO: warning

        data.append((
            set_ac,
            supfam["name"],
            supfam["description"],
            "cdd",
            1,
            json.dumps(s["relationships"])
        ))

    # logging.info("loading PANTHER superfamilies")
    # for s in oracle.get_profile_alignments(ora_uri, "panther"):
    #     data.append((
    #         s["accession"],
    #         s["name"],          # None
    #         s["description"],   # None
    #         "panther",
    #         1,
    #         json.dumps(s["relationships"])
    #     ))

    logging.info("loading PIRSF superfamilies")
    for s in oracle.get_profile_alignments(ora_uri, "pirsf"):
        data.append((
            s["accession"],
            s["name"],          # None
            s["description"],   # None
            "pirsf",
            1,
            json.dumps(s["relationships"])
        ))

    logging.info("inserting sets")
    con, cur = dbms.connect(my_uri)
    for i in range(0, len(data), chunk_size):
        cur.executemany(
            """
            INSERT INTO webfront_set (
              accession, name, description, source_database, is_set,
              relationships
            ) VALUES (%s, %s, %s, %s, %s, %s)
            """,
            data[i:i+chunk_size]
        )

    cur.close()
    con.commit()
    con.close()


def insert_proteins(uri, src_proteins, src_sequences, src_misc,
                    src_names, src_comments, src_proteomes,
                    src_residues, src_structures, src_features,
                    src_matches, chunk_size=100000, limit=0):
    logging.info("starting")

    # MySQL data
    taxa = get_taxa(uri, lineage=False)
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

    proteins = io.Store(src_proteins)
    protein2sequence = io.Store(src_sequences)
    protein2misc = io.Store(src_misc)
    protein2names = io.Store(src_names)
    protein2comments = io.Store(src_comments)
    protein2proteome = io.Store(src_proteomes)
    protein2residues = io.Store(src_residues)
    protein2structures = io.Store(src_structures)
    protein2features = io.Store(src_features)
    protein2matches = io.Store(src_matches)

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
        tax_id = protein["taxon"]
        try:
            taxon = taxa[tax_id]
        except KeyError:
            continue

        if protein["length"] <= 100:
            size = "small"
        elif protein["length"] <= 1000:
            size = "medium"
        else:
            size = "large"

        evidence, gene = protein2misc.get(acc, (None, None))
        if not evidence:
            continue
        name, other_names = protein2names.get(acc, (None, None))
        upid = protein2proteome.get(acc)

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
            json.dumps(taxon),
            name,
            json.dumps(other_names),
            json.dumps(protein2comments.get(acc, [])),
            protein2sequence[acc],
            protein["length"],
            size,
            upid,
            gene,
            json.dumps(list(go_terms.values())),
            evidence,
            "reviewed" if protein["isReviewed"] else "unreviewed",
            json.dumps(protein2residues.get(acc, {})),
            1 if protein["isFrag"] else 0,
            json.dumps(protein2structures.get(acc, {})),
            tax_id,
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


def get_taxa(uri: str, lineage: bool=False) -> dict:
    con, cur = dbms.connect(uri)
    if lineage:
        cur.execute(
            """
            SELECT accession, scientific_name, full_name, lineage, rank
            FROM webfront_taxonomy
            """
        )
        cols = ("taxId", "scientificName", "fullName", "lineage", "rank")
        taxa = {row[0]: dict(zip(cols, row)) for row in cur}
    else:
        cur.execute(
            """
            SELECT accession, scientific_name, full_name
            FROM webfront_taxonomy
            """
        )
        cols = ("taxId", "scientificName", "fullName")
        taxa = {row[0]: dict(zip(cols, row)) for row in cur}

    cur.close()
    con.close()

    return taxa


def get_proteomes(uri: str) -> dict:
    con, cur = dbms.connect(uri)

    cur.execute(
        """
        SELECT accession, name, is_reference, strain, assembly, taxonomy_id
        FROM webfront_proteome
        """
    )

    proteomes = {}
    for row in cur:
        proteomes[row[0]] = {
            'name': row[1],
            'is_reference': bool(row[2]),
            'strain': row[3],
            'assembly': row[4],
            'taxon': row[5]
        }

    cur.close()
    con.close()

    return proteomes


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


def find_node(node, accession, relations=[]):
    """
    Find a entry (node) in its hierarchy tree
    and store its relations (ancestors + direct children)
    """
    if node['accession'] == accession:
        relations += [child['accession'] for child in node['children']]
        return node

    for child in node['children']:
        child = find_node(child, accession, relations)

        if child:
            relations.append(node['accession'])
            return child

    return None


def parse_json(value):
    if value is None:
        return value
    else:
        json.loads(value)


def get_entries(uri: str, has_is_alive: bool=True) -> dict:
    query = """
        SELECT
            accession, source_database, entry_date, description,
            integrated_id, name, type, short_name, member_databases,
            go_terms, literature, cross_references, hierarchy
        FROM webfront_entry
    """

    if has_is_alive:
        query += "WHERE is_alive = 1"

    con, cur = dbms.connect(uri)
    cur.execute(query)

    entries = {}
    for row in cur:
        accession = row[0]
        relations = []
        hierarchy = parse_json(row[12])
        if hierarchy:
            find_node(hierarchy, accession, relations)
            _hierarchy = hierarchy.get("accession")
        else:
            _hierarchy = None

        entries[accession] = {
            "accession": accession,
            "database": row[1],
            "date": row[2],
            "descriptions": parse_json(row[3]),
            "integrated": row[4],
            "name": row[5],
            "type": row[6],
            "short_name": row[7],
            "member_databases": parse_json(row[8]),
            "go_terms": parse_json(row[9]),
            "citations": parse_json(row[10]),
            "cross_references": parse_json(row[11]),
            "root": _hierarchy,
            "relations": relations
        }

    cur.close()
    con.close()

    return entries


def get_sets(uri: str) -> dict:
    con, cur = dbms.connect(uri, sscursor=True)
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
                       src_proteomes, version, release_date):
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
    taxa = set(get_taxa(stg_uri, lineage=False))

    # Integrated signatures
    integrated = {
        acc
        for acc, e in get_entries(stg_uri).items()
        if e["integrated"]
    }

    # Protein to PDBe structures
    protein2pdb = {}
    for pdb_id, s in get_structures(stg_uri).items():
        for acc in s["proteins"]:
            if acc in protein2pdb:
                protein2pdb[acc].add(pdb_id)
            else:
                protein2pdb[acc] = {pdb_id}

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

    proteins = io.Store(src_proteins)
    protein2matches = io.Store(src_matches)
    protein2proteome = io.Store(src_proteomes)

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

                    upid = protein2proteome.get(acc)
                    if upid:
                        interpro_proteomes.add(upid)

                    interpro_structures |= protein2pdb.get(acc, set())
                    break

        n_proteins += 1
        if not n_proteins % 1000000:
            logging.info("{:>12,} ({:.0f} proteins/sec)".format(
                n_proteins, n_proteins / (time.time() - ts)
            ))

    proteins.close()
    protein2matches.close()
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
    for e in get_entries(rel_uri, False).values():
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
                datetime.strptime(release_date, "%Y-%m-%d"),
                json.dumps(notes)
            )
        )

    cur.close()
    con.commit()
    con.close()
    logging.info("complete")


def reduce(src: dict):
    dst = {}
    for k, v in src.items():
        if isinstance(v, dict):
            dst[k] = reduce(v)
        elif isinstance(v, set):
            dst[k] = len(v)
        else:
            dst[k] = v

    return dst


def update_taxa_counts(uri: str, src_taxa: str, processes: int=1):
    logging.info("loading taxa")
    taxa = {}
    with io.Store(src_taxa) as store:
        for tax_id, xrefs in store.iter(processes):
            for e in xrefs["proteins"]:
                # xrefs["proteins"] is a set of one item
                xrefs["proteins_total"] = e

            taxa[tax_id] = xrefs

    logging.info("propagating cross-references to taxa lineage")
    cnt = 0
    all_taxa = set()
    for tax_id, t in get_taxa(uri, lineage=True).items():
        all_taxa.add(tax_id)

        try:
            taxon = taxa[tax_id]
        except KeyError:
            continue

        n_proteins = taxon["proteins"].pop()

        # lineage stored as a string in MySQL (string include the taxon)
        # -2: first item to include (second to last; last is current taxon)
        # -1: negative step (reverse list)
        lineage = t["lineage"].strip().split()[-2::-1]
        for parent_id in lineage:
            try:
                parent = taxa[parent_id]
            except KeyError:
                parent = {
                    "domains": set(),
                    "entries": {},
                    "proteomes": set(),
                    "proteins": {0},
                    "proteins_total": 0,
                    "sets": set(),
                    "structures": set()
                }

            parent["proteins_total"] += n_proteins

            for _type in ("domains", "proteomes", "sets", "structures"):
                try:
                    accessions = taxon[_type]
                except KeyError:
                    # Type absent in taxon (e.g. no cross-refs)
                    accessions = taxon[_type] = set()
                finally:
                    if _type in parent:
                        parent[_type] |= set(accessions)
                    else:
                        parent[_type] = set(accessions)

            try:
                entries = taxon["entries"]
            except KeyError:
                entries = taxon["entries"] = {}
            finally:
                if "entries" not in parent:
                    parent["entries"] = {}

                for entry_db, db_entries in entries.items():
                    if entry_db in parent["entries"]:
                        parent["entries"][entry_db] |= set(db_entries)
                    else:
                        parent["entries"][entry_db] = set(db_entries)

            taxa[parent_id] = parent

        taxa[tax_id] = taxon

        cnt += 1
        if not cnt % 100000:
            logging.info(cnt)

    logging.info("updating webfront_taxonomy")
    con, cur = dbms.connect(uri)
    for tax_id, taxon in taxa.items():
        try:
            all_taxa.remove(tax_id)  # todo: check if needed
        except KeyError:
            continue

        counts = reduce(taxon)
        counts["proteins"] = counts.pop("proteins_total")

        try:
            counts["entries"]["total"] = sum(counts["entries"].values())
        except KeyError:
            counts["entries"] = {"total": 0}

        cur.execute(
            """
            UPDATE webfront_taxonomy
            SET counts = %s
            WHERE accession = %s
            """,
            (json.dumps(counts), tax_id)
        )

    for tax_id in all_taxa:
        cur.execute(
            """
            UPDATE webfront_taxonomy
            SET counts = %s
            WHERE accession = %s
            """,
            (json.dumps({
                "domains": 0,
                "entries": {"total": 0},
                "proteomes": 0,
                "proteins": 0,
                "sets": 0,
                "structures": 0
            }), tax_id)
        )

    con.commit()
    cur.close()
    con.close()


def _update_taxa_counts(uri: str, src_taxa: str, processes: int=1):
    with io.KVdb(cache_size=10000) as taxa:
        logging.info("loading taxa")
        with io.Store(src_taxa) as store:
            for tax_id, xrefs in store.iter(processes):
                for e in xrefs["proteins"]:
                    # xrefs["proteins"] is a set of one item
                    xrefs["proteins_total"] = e
                taxa[tax_id] = xrefs

        logging.info("propagating cross-references to taxa lineage")
        all_taxa = set()
        cnt = 0
        for tax_id, t in get_taxa(uri, lineage=True).items():
            all_taxa.add(tax_id)

            try:
                taxon = taxa[tax_id]
            except KeyError:
                continue

            n_proteins = taxon["proteins"].pop()

            # lineage stored as a string in MySQL (string include the taxon)
            # -2: first item to include (second to last; last is current taxon)
            # -1: negative step (reverse list)
            lineage = t["lineage"].strip().split()[-2::-1]
            for parent_id in lineage:
                try:
                    parent = taxa[parent_id]
                except KeyError:
                    parent = {
                        "domains": set(),
                        "entries": {},
                        "proteomes": set(),
                        "proteins": {0},
                        "proteins_total": 0,
                        "sets": set(),
                        "structures": set()
                    }

                parent["proteins_total"] += n_proteins

                for _type in ("domains", "proteomes", "sets", "structures"):
                    try:
                        accessions = taxon[_type]
                    except KeyError:
                        # Type absent in taxon (e.g. no cross-refs)
                        accessions = taxon[_type] = set()
                    finally:
                        if _type in parent:
                            parent[_type] |= set(accessions)
                        else:
                            parent[_type] = set(accessions)

                try:
                    entries = taxon["entries"]
                except KeyError:
                    entries = taxon["entries"] = {}
                finally:
                    if "entries" not in parent:
                        parent["entries"] = {}

                    for entry_db, db_entries in entries.items():
                        if entry_db in parent["entries"]:
                            parent["entries"][entry_db] |= set(db_entries)
                        else:
                            parent["entries"][entry_db] = set(db_entries)

                # Write back parent to DB
                taxa[parent_id] = parent

            # Write back taxon to DB
            taxa[tax_id] = taxon

            cnt += 1
            if not cnt % 100000:
                logging.info(cnt)

        logging.info("updating webfront_taxonomy")
        con, cur = dbms.connect(uri)
        for tax_id, taxon in taxa:
            try:
                all_taxa.remove(tax_id)  # todo: check if needed
            except KeyError:
                continue

            counts = reduce(taxon)
            counts["proteins"] = counts.pop("proteins_total")

            try:
                counts["entries"]["total"] = sum(counts["entries"].values())
            except KeyError:
                counts["entries"] = {"total": 0}

            cur.execute(
                """
                UPDATE webfront_taxonomy
                SET counts = %s
                WHERE accession = %s
                """,
                (json.dumps(counts), tax_id)
            )

        logging.info("database size: {:,}".format(taxa.getsize()))

    for tax_id in all_taxa:
        cur.execute(
            """
            UPDATE webfront_taxonomy
            SET counts = %s
            WHERE accession = %s
            """,
            (json.dumps({
                "domains": 0,
                "entries": {"total": 0},
                "proteomes": 0,
                "proteins": 0,
                "sets": 0,
                "structures": 0
            }), tax_id)
        )

    con.commit()
    cur.close()
    con.close()


def update_proteomes_counts(uri: str, src_proteomes: str, processes: int=1):
    con, cur = dbms.connect(uri)

    logging.info("updating webfront_proteome")
    proteomes = set(get_proteomes(uri))
    with io.Store(src_proteomes) as store:
        for upid, xrefs in store.iter(processes):
            try:
                proteomes.remove(upid)
            except KeyError:
                continue

            counts = reduce(xrefs)
            try:
                counts["entries"]["total"] = sum(counts["entries"].values())
            except KeyError:
                counts["entries"] = {"total": 0}

            cur.execute(
                """
                UPDATE webfront_proteome
                SET counts = %s
                WHERE accession = %s
                """,
                (json.dumps(counts), upid)
            )

    for upid in proteomes:
        cur.execute(
            """
            UPDATE webfront_proteome
            SET counts = %s
            WHERE accession = %s
            """,
            (json.dumps({
                "domains": 0,
                "entries": {"total": 0},
                "proteins": 0,
                "sets": 0,
                "structures": 0,
                "taxa": 1,
            }), upid)
        )

    con.commit()
    cur.close()
    con.close()


def update_structures_counts(uri: str, src_structures: str, processes: int=1):
    con, cur = dbms.connect(uri)
    logging.info("updating webfront_structure")
    structures = set(get_structures(uri))
    with io.Store(src_structures) as store:
        for pdb_id, xrefs in store.iter(processes):
            try:
                structures.remove(pdb_id)
            except KeyError:
                continue

            counts = reduce(xrefs)
            try:
                counts["entries"]["total"] = sum(counts["entries"].values())
            except KeyError:
                counts["entries"] = {"total": 0}

            cur.execute(
                """
                UPDATE webfront_structure
                SET counts = %s
                WHERE accession = %s
                """,
                (json.dumps(counts), pdb_id)
            )

    for pdb_id in structures:
        cur.execute(
            """
            UPDATE webfront_structure
            SET counts = %s
            WHERE accession = %s
            """,
            (json.dumps({
                "domains": 0,
                "entries": {"total": 0},
                "proteins": 0,
                "proteomes": 0,
                "sets": 0,
                "taxa": 0,
            }), pdb_id)
        )

    con.commit()
    cur.close()
    con.close()


def update_entries_sets_counts(uri: str, src_entries: str, processes: int=1):
    logging.info("updating webfront_entry")
    sets = get_sets(uri)
    entry2set = {}
    for set_ac, s in sets.items():
        for entry_ac in s["members"]:
            entry2set[entry_ac] = set_ac

    entries = {}
    all_entries = set(get_entries(uri, has_is_alive=False))

    con, cur = dbms.connect(uri)
    with io.Store(src_entries) as store:
        cnt = 0
        for entry_ac, xrefs in store.iter(processes):
            all_entries.remove(entry_ac)

            counts = reduce(xrefs)
            counts["matches"] = xrefs["matches"].pop()  # set of one item
            if entry2set.get(entry_ac):
                entries[entry_ac] = xrefs
                counts["sets"] = 1
            else:
                counts["sets"] = 0

            cur.execute(
                """
                UPDATE webfront_entry
                SET counts = %s
                WHERE accession = %s
                """,
                (json.dumps(counts), entry_ac)
            )

            cnt += 1
            if not cnt % 10000:
                logging.info(cnt)

    for entry_ac in all_entries:
        cur.execute(
            """
            UPDATE webfront_entry
            SET counts = %s
            WHERE accession = %s
            """,
            (json.dumps({
                "matches": 0,
                "proteins": 0,
                "proteomes": 0,
                "sets": 1 if entry_ac in entry2set else 0,
                "structures": 0,
                "taxa": 0
            }), entry_ac)
        )

        cnt += 1
        if not cnt % 10000:
            logging.info(cnt)

    logging.info("updating webfront_set")
    for set_ac, s in sets.items():
        xrefs = {
            "domains": set(),
            "entries": {
                s["database"]: len(s["members"]),
                "total": len(s["members"])
            },
            "proteins": set(),
            "proteomes": set(),
            "structures": set(),
            "taxa": set()
        }

        for entry_ac in s["members"]:
            try:
                entry = entries.pop(entry_ac)
            except KeyError:
                continue

            for _type in ("domains", "proteins", "proteomes", "structures", "taxa"):
                try:
                    accessions = entry[_type]
                except KeyError:
                    # Type absent
                    continue
                else:
                    xrefs[_type] |= accessions

        cur.execute(
            """
            UPDATE webfront_set
            SET counts = %s
            WHERE accession = %s
            """,
            (json.dumps(reduce(xrefs)), set_ac)
        )

    con.commit()
    cur.close()
    con.close()


def _update_entries_sets_counts(uri: str, src_entries: str, processes: int=1):
    logging.info("updating webfront_entry")
    sets = get_sets(uri)
    entry2set = {}
    for set_ac, s in sets.items():
        for entry_ac in s["members"]:
            entry2set[entry_ac] = set_ac

    all_entries = set(get_entries(uri, has_is_alive=False))
    with io.KVdb(cache_size=10) as entries:
        con, cur = dbms.connect(uri)

        with io.Store(src_entries) as store:
            cnt = 0
            for entry_ac, xrefs in store.iter(processes):
                all_entries.remove(entry_ac)

                counts = reduce(xrefs)
                counts["matches"] = xrefs["matches"].pop()  # set of one item
                if entry2set.get(entry_ac):
                    entries[entry_ac] = xrefs
                    counts["sets"] = 1
                else:
                    counts["sets"] = 0

                cur.execute(
                    """
                    UPDATE webfront_entry
                    SET counts = %s
                    WHERE accession = %s
                    """,
                    (json.dumps(counts), entry_ac)
                )

                cnt += 1
                if not cnt % 10000:
                    logging.info(cnt)

        for entry_ac in all_entries:
            cur.execute(
                """
                UPDATE webfront_entry
                SET counts = %s
                WHERE accession = %s
                """,
                (json.dumps({
                    "matches": 0,
                    "proteins": 0,
                    "proteomes": 0,
                    "sets": 1 if entry_ac in entry2set else 0,
                    "structures": 0,
                    "taxa": 0
                }), entry_ac)
            )

            cnt += 1
            if not cnt % 10000:
                logging.info(cnt)

        logging.info("updating webfront_set")
        for set_ac, s in sets.items():
            xrefs = {
                "domains": set(),
                "entries": {
                    s["database"]: len(s["members"]),
                    "total": len(s["members"])
                },
                "proteins": set(),
                "proteomes": set(),
                "structures": set(),
                "taxa": set()
            }

            for entry_ac in s["members"]:
                try:
                    entry = entries[entry_ac]
                except KeyError:
                    continue

                for _type in ("domains", "proteins", "proteomes", "structures", "taxa"):
                    try:
                        accessions = entry[_type]
                    except KeyError:
                        # Type absent
                        continue
                    else:
                        xrefs[_type] |= accessions

            cur.execute(
                """
                UPDATE webfront_set
                SET counts = %s
                WHERE accession = %s
                """,
                (json.dumps(reduce(xrefs)), set_ac)
            )

        con.commit()
        cur.close()
        con.close()
        logging.info("database size: {:,}".format(entries.getsize()))


def update_counts(uri: str, src_entries: str, src_proteomes: str,
                  src_structures: str, src_taxa: str, processes: int=1):
    update_taxa_counts(uri, src_taxa, processes)
    update_proteomes_counts(uri, src_proteomes, processes)
    update_structures_counts(uri, src_structures, processes)
    update_entries_sets_counts(uri, src_entries, processes)
    logging.info("complete")
