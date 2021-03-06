# -*- coding: utf-8 -*-

import os
from typing import Optional

import MySQLdb

from interpro7dw import logger
from interpro7dw.ebi.interpro.utils import Table, overlaps_pdb_chain
from interpro7dw.utils import DumpFile, DirectoryTree, Store, merge_dumps
from interpro7dw.utils import loadobj, deepupdate, url2dict
from .utils import jsonify, reduce


def init_xrefs() -> dict:
    return {
        "entries": {},
        "proteins": {
            "all": 0,
            "databases": {},  # grouped by entry database
            "entries": {}     # grouped by entry
        },
        "proteomes": set(),
        "structures": {
            "all": set(),
            "entries": {}     # overlapping with these entries
        }
    }


def dump_xrefs(xrefs: dict, taxonomy: dict, output: str):
    # Init all taxa
    final_xrefs = {}
    for taxon_id in taxonomy:
        final_xrefs[taxon_id] = init_xrefs()

    while xrefs:
        taxon_id, taxon_xrefs = xrefs.popitem()

        for node_id in taxonomy[taxon_id]["lineage"]:
            deepupdate(taxon_xrefs, final_xrefs[node_id], replace=False)

    with DumpFile(output, compress=True) as f:
        for taxon_id in sorted(final_xrefs):
            f.dump((taxon_id, final_xrefs[taxon_id]))


def insert_taxonomy(p_entries: str, p_proteins: str, p_structures: str,
                    p_taxonomy: str, p_uniprot2matches: str,
                    p_uniprot2proteome: str, stg_url: str,
                    p_interpro2taxonomy: str, tmpdir: Optional[str] = None):
    logger.info("preparing data")
    dt = DirectoryTree(tmpdir)
    entries = loadobj(p_entries)
    taxonomy = loadobj(p_taxonomy)
    uniprot2pdbe = {}
    for pdb_id, entry in loadobj(p_structures).items():
        for uniprot_acc, chains in entry["proteins"].items():
            try:
                uniprot2pdbe[uniprot_acc][pdb_id] = chains
            except KeyError:
                uniprot2pdbe[uniprot_acc] = {pdb_id: chains}

    proteins = Store(p_proteins)
    u2matches = Store(p_uniprot2matches)
    u2proteome = Store(p_uniprot2proteome)

    logger.info("starting")
    i = 0
    xrefs = {}
    files = []
    for uniprot_acc, info in proteins.items():
        taxon_id = info["taxid"]

        try:
            taxon = xrefs[taxon_id]
        except KeyError:
            taxon = xrefs[taxon_id] = init_xrefs()

        try:
            proteome_id = u2proteome[uniprot_acc]
        except KeyError:
            pass
        else:
            taxon["proteomes"].add(proteome_id)

        taxon["proteins"]["all"] += 1

        protein_structures = uniprot2pdbe.get(uniprot_acc, {})

        # Add structures to taxon, regardless of entry matches
        taxon["structures"]["all"] |= set(protein_structures.keys())

        databases = set()
        for entry_acc, locations in u2matches.get(uniprot_acc, {}).items():
            entry = entries[entry_acc]
            database = entry.database

            try:
                taxon["entries"][database].add(entry_acc)
            except KeyError:
                taxon["entries"][database] = {entry_acc}

            if database not in databases:
                # Counting the protein *once* per database
                databases.add(database)
                try:
                    taxon["proteins"]["databases"][database] += 1
                except KeyError:
                    taxon["proteins"]["databases"][database] = 1

            try:
                taxon["proteins"]["entries"][entry_acc] += 1
            except KeyError:
                taxon["proteins"]["entries"][entry_acc] = 1

            for pdb_id, chains in protein_structures.items():
                for chain_id, segments in chains.items():
                    if overlaps_pdb_chain(locations, segments):
                        try:
                            taxon["structures"]["entries"][entry_acc].add(pdb_id)
                        except KeyError:
                            taxon["structures"]["entries"][entry_acc] = {pdb_id}

                        break  # Skip other chains

        i += 1
        if not i % 1000000:
            output = dt.mktemp()
            dump_xrefs(xrefs, taxonomy, output)
            files.append(output)
            xrefs = {}

            if not i % 10000000:
                logger.info(f"{i:>12,}")

    if xrefs:
        output = dt.mktemp()
        dump_xrefs(xrefs, taxonomy, output)
        files.append(output)
        xrefs = {}

    logger.info(f"{i:>12,}")
    logger.info(f"temporary files: "
                f"{sum(map(os.path.getsize, files))/1024/1024:.0f} MB")

    proteins.close()
    u2matches.close()
    u2proteome.close()

    logger.info("populating taxonomy tables")
    con = MySQLdb.connect(**url2dict(stg_url), charset="utf8mb4")
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
            children LONGTEXT,
            counts LONGTEXT NOT NULL
        ) CHARSET=utf8mb4 DEFAULT COLLATE=utf8mb4_unicode_ci
        """
    )
    cur.execute("DROP TABLE IF EXISTS webfront_taxonomyperentry")
    cur.execute(
        """
        CREATE TABLE webfront_taxonomyperentry
        (
          id INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
          tax_id VARCHAR(20) NOT NULL,
          entry_acc VARCHAR(25) NOT NULL,
          counts LONGTEXT NULL NULL
        ) CHARSET=utf8mb4 DEFAULT COLLATE=utf8mb4_unicode_ci
        """
    )
    cur.execute("DROP TABLE IF EXISTS webfront_taxonomyperentrydb")
    cur.execute(
        """
        CREATE TABLE webfront_taxonomyperentrydb
        (
          id INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
          tax_id VARCHAR(20) NOT NULL,
          source_database VARCHAR(10) NOT NULL,
          counts LONGTEXT NOT NULL
        ) CHARSET=utf8mb4 DEFAULT COLLATE=utf8mb4_unicode_ci
        """
    )
    cur.close()

    table = Table(con, query="""
        INSERT INTO webfront_taxonomy VALUES (%s, %s, %s, %s, %s, %s, %s, %s) 
    """)
    per_entry = Table(con, query="""
        INSERT INTO webfront_taxonomyperentry (tax_id,entry_acc,counts)
        VALUES (%s, %s, %s) 
    """)
    per_database = Table(con, query="""
        INSERT INTO webfront_taxonomyperentrydb (tax_id,source_database,counts)
        VALUES (%s, %s, %s) 
    """)

    with DumpFile(p_interpro2taxonomy, compress=True) as interpro2taxonomy:
        interpro_entries = {
            entry.accession
            for entry in entries.values()
            if entry.database == "interpro" and not entry.is_deleted
        }

        i = 0
        for taxon_id, taxon_xrefs in merge_dumps(files):
            taxon = taxonomy[taxon_id]

            protein_counts = taxon_xrefs.pop("proteins")
            structure_counts = taxon_xrefs.pop("structures")
            counts = reduce(taxon_xrefs)

            # Add total protein count (not grouped by database/entry)
            counts["proteins"] = protein_counts["all"]

            # Add total structure count
            counts["structures"] = len(structure_counts["all"])

            # Add total entry count (not grouped by database)
            counts["entries"]["total"] = sum(counts["entries"].values())

            table.insert((
                taxon_id,
                taxon["sci_name"],
                taxon["full_name"],
                f" {' '.join(taxon['lineage'])} ",
                taxon["parent"],
                taxon["rank"],
                jsonify(taxon["children"]),
                jsonify(counts)
            ))

            # Remove the 'entry' property
            # (no needed for webfront_taxonomyperentry)
            entry_counts = counts.pop("entries")

            database_structures = {}
            for entry_acc, count in protein_counts["entries"].items():
                if entry_acc in interpro_entries:
                    interpro2taxonomy.dump((entry_acc, taxon_id, count))

                counts["proteins"] = count

                try:
                    entry_structures = structure_counts["entries"][entry_acc]
                except KeyError:
                    counts["structures"] = 0
                else:
                    counts["structures"] = len(entry_structures)

                    database = entries[entry_acc].database
                    try:
                        database_structures[database] |= entry_structures
                    except KeyError:
                        database_structures[database] = entry_structures.copy()
                finally:
                    per_entry.insert((taxon_id, entry_acc, jsonify(counts)))

            for database, count in protein_counts["databases"].items():
                counts.update({
                    "entries": entry_counts[database],
                    "proteins": count,
                    "structures": len(database_structures.get(database, []))
                })
                per_database.insert((taxon_id, database, jsonify(counts)))

            i += 1
            if not i % 100000:
                logger.info(f"{i:>12,}")

        logger.info(f"{i:>12,}")

    table.close()
    per_entry.close()
    per_database.close()
    con.commit()

    dt.remove()

    logger.info("indexing")
    cur = con.cursor()
    cur.execute(
        """
        CREATE INDEX i_webfront_taxonomyperentry_tax
        ON webfront_taxonomyperentry (tax_id)
        """
    )
    cur.execute(
        """
        CREATE INDEX i_webfront_taxonomyperentry_entry
        ON webfront_taxonomyperentry (entry_acc)
        """
    )
    cur.execute(
        """
        CREATE INDEX i_webfront_taxonomyperentrydb_tax
        ON webfront_taxonomyperentrydb (tax_id)
        """
    )
    cur.execute(
        """
        CREATE INDEX i_webfront_taxonomyperentrydb_database
        ON webfront_taxonomyperentrydb (source_database)
        """
    )
    cur.close()
    con.close()
    logger.info("complete")
