# -*- coding: utf-8 -*-

import json
import os
from typing import Generator, Optional

import MySQLdb
import MySQLdb.cursors

from i7dw import io, logger
from i7dw.interpro import Table, mysql, oracle
from .utils import parse_url, reduce


def insert_taxa(my_url: str, ora_url: str):
    query = """
        INSERT INTO webfront_taxonomy (accession, scientific_name, full_name,
                                       lineage, parent_id, rank, children)
        VALUES (%s, %s, %s, %s, %s, %s, %s)
    """

    con = MySQLdb.connect(**parse_url(my_url), charset="utf8")
    with Table(con, query) as table:
        for taxon in oracle.get_taxa(ora_url):
            table.insert((
                taxon["id"],
                taxon["scientific_name"],
                taxon["full_name"],
                # leading/trailing whitespaces are important from API queries
                ' ' + ' '.join(taxon["lineage"]) + ' ',
                taxon["parent_id"],
                taxon["rank"],
                json.dumps(taxon["children"])
            ))

    con.commit()
    con.close()


def iter_taxa(url: str, lineage: bool=False) -> Generator[dict, None, None]:
    con = MySQLdb.connect(**parse_url(url), charset="utf8")
    cur = MySQLdb.cursors.SSCursor(con)
    if lineage:
        cur.execute(
            """
            SELECT accession, scientific_name, full_name, lineage, rank
            FROM webfront_taxonomy
            """
        )
        for row in cur:
            yield {
                "id": row[0],
                "scientific_name": row[1],
                "full_name": row[2],
                "lineage": row[3].strip().split(),
                "rank": row[4]
            }
    else:
        cur.execute(
            """
            SELECT accession, scientific_name, full_name
            FROM webfront_taxonomy
            """
        )
        for row in cur:
            yield {
                "id": row[0],
                "scientific_name": row[1],
                "full_name": row[2]
            }

    cur.close()
    con.close()


def update_counts(url: str, src_proteins: str, src_proteomes:str,
                  src_matches: str, buffer_size: int=100000,
                  tmpdir: Optional[str]=None):
    # Get required MySQL data
    logger.info("loading data")
    entries = mysql.entries.get_entries(url)
    entries = {k: v["database"] for k, v in entries.items()}

    structures = {}
    for s in mysql.structures.iter_structures(url):
        pdbe_id = s["accession"]
        for protein_acc, chains in s["proteins"].items():
            try:
                structures[protein_acc].add(pdbe_id)
            except KeyError:
                structures[protein_acc] = {pdbe_id}

    taxa_db = io.mktemp(suffix=".db", dir=tmpdir)
    with io.KVdb(taxa_db, writeback=True) as kvdb:
        lineages = {}
        for taxon in iter_taxa(url, lineage=True):
            # `lineage` contains the taxon ID: do not consider it
            lineages[taxon["id"]] = taxon["lineage"][:-1]

            # Initiate Kvdb so all taxa are there
            kvdb[taxon["id"]] = {
                "entries": {},
                "proteins": {"all": 0, "databases": {}, "entries": {}},
                "proteomes": set(),
                "structures": set()
            }

        kvdb.sync()

        logger.info("starting")
        proteins = io.Store(src_proteins)
        proteomes = io.Store(src_proteomes)
        matches = io.Store(src_matches)

        i_progress = 0
        for protein_acc, protein_info in proteins:
            taxon_id = protein_info["taxon"]
            protein_counts = {"all": 1, "databases": {}, "entries": {}}
            protein_entries = {}

            for entry_acc in matches.get(protein_acc, {}):
                database = entries[entry_acc]

                try:
                    protein_entries[database].add(entry_acc)
                except KeyError:
                    protein_entries[database] = {entry_acc}

                protein_counts["databases"][database] = 1
                protein_counts["entries"][entry_acc] = 1

            upid = proteomes.get(protein_acc)

            # Update cross-reference for protein's taxon
            node = kvdb[taxon_id]
            io.traverse({
                "entries": protein_entries,
                "proteins": protein_counts,
                "proteomes": {upid} if upid else set(),
                "structures": structures.get(protein_acc, set())
            }, node, replace=False)
            kvdb[taxon_id] = node

            # Update protein counts for taxon's ancestors
            for tax_id in lineages[taxon_id]:
                node = kvdb[tax_id]
                io.traverse({"proteins": protein_counts}, node, replace=False)
                kvdb[tax_id] = node

            i_progress += 1
            if not i_progress % 10000000:
                logger.info(f"{i_progress:>12,}")

            if len(kvdb.cache) >= buffer_size:
                kvdb.sync()

        proteins.close()
        proteomes.close()
        matches.close()
        logger.info(f"{i_progress:>12,}")

        logger.info("propagating to lineage")
        i_progress = 0
        for taxon_id, ancestors in lineages.items():
            taxon = kvdb[taxon_id]
            del taxon["proteins"]  # proteins were already propagated

            for tax_id in ancestors:
                node = kvdb[tax_id]
                io.traverse(taxon, node)
                kvdb[tax_id] = node

            if len(kvdb.cache) >= buffer_size:
                kvdb.sync()

            i_progress += 1
            if not i_progress % 100000:
                logger.info(f"{i_progress:>10,}")

        kvdb.sync()
        lineages.clear()
        logger.info(f"{i_progress:>10,}")

        logger.info("updating MySQL tables")
        con = MySQLdb.connect(**parse_url(url), charset="utf8")
        table1 = Table(con, query="UPDATE webfront_taxonomy SET counts = %s "
                                  "WHERE accession = %s")
        table2 = Table(con, query="INSERT INTO webfront_taxonomyperentry  "
                                  "(tax_id, entry_acc, counts)  "
                                  "VALUES (%s, %s, %s)")
        table3 = Table(con, query="INSERT INTO webfront_taxonomyperentrydb "
                                  "(tax_id, source_database, counts) "
                                  "VALUES (%s, %s, %s)")

        i_progress = 0
        for tax_id, xrefs in kvdb:
            protein_counts = xrefs.pop("proteins")

            # Counts for `webfront_taxonomy`
            counts = reduce(xrefs)
            # All proteins (not filtered by database/entry)
            counts["proteins"] = protein_counts["all"]
            # Total number of entries (all databases)
            counts["entries"]["total"] = sum(counts["entries"].values())
            table1.update((json.dumps(counts), tax_id))

            # We do not need `entries` for other tables
            del counts["entries"]

            # Counts for `webfront_taxonomyperentry`
            for entry_acc, cnt in protein_counts["entries"].items():
                counts["proteins"] = cnt
                table2.insert((tax_id, entry_acc, json.dumps(counts)))

            # Counts for `webfront_taxonomyperentrydb`
            for database, cnt in protein_counts["databases"].items():
                counts["proteins"] = cnt
                table3.insert((tax_id, database, json.dumps(counts)))

            i_progress += 1
            if not i_progress % 100000:
                logger.info(f"{i_progress:>10,}")

        for t in (table1, table2, table3):
            t.close()

        logger.info(f"{i_progress:>10,}")
        con.commit()

        logger.info("creating indexes")
        cur = con.cursor()
        cur.execute(
            """
            CREATE UNIQUE INDEX i_webfront_taxonomyperentry_tax
            ON webfront_taxonomyperentry (tax_id)
            """
        )
        cur.execute(
            """
            CREATE UNIQUE INDEX i_webfront_taxonomyperentry_entry
            ON webfront_taxonomyperentry (entry_acc)
            """
        )
        cur.execute(
            """
            CREATE UNIQUE INDEX i_webfront_taxonomyperentrydb_tax
            ON webfront_taxonomyperentrydb (tax_id)
            """
        )
        cur.execute(
            """
            CREATE UNIQUE INDEX i_webfront_taxonomyperentrydb_database
            ON webfront_taxonomyperentrydb (source_database)
            """
        )
        cur.close()
        con.close()

    logger.info(f"disk usage: {os.path.getsize(taxa_db)/1024/1024:.0f} MB")
    os.remove(taxa_db)
