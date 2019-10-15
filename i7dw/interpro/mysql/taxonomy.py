# -*- coding: utf-8 -*-

import json
from typing import Generator, Optional

import MySQLdb
import MySQLdb.cursors

from i7dw import io, logger
from i7dw.interpro import DomainArchitecture, Table, mysql, oracle
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


def _export(url: str, src_proteins: str, src_proteomes: str, src_matches: str,
            sync_frequency: int, tmpdir: Optional[str]) -> io.Store:
    # Get required MySQL data
    logger.info("loading data")
    entries = mysql.entries.get_entries(url)
    dom_arch = DomainArchitecture(entries)

    entry_set = {}
    for s in mysql.entries.iter_sets(url):
        set_acc = s["accession"]
        for entry_acc in s["members"]:
            entry_set[entry_acc] = set_acc

    lineages = {}
    for taxon in iter_taxa(url, lineage=True):
        lineages[taxon["id"]] = taxon["lineage"]

    structures = {}
    for s in mysql.structures.iter_structures(url):
        pdbe_id = s["accession"]
        for protein_acc, chains in s["proteins"].items():
            try:
                structures[protein_acc].add(pdbe_id)
            except KeyError:
                structures[protein_acc] = {pdbe_id}

    logger.info("starting")
    proteins = io.Store(src_proteins)
    proteomes = io.Store(src_proteomes)
    matches = io.Store(src_matches)
    taxa = io.Store(keys=io.Store.chunk_keys(lineages, 100), tmpdir=tmpdir)

    cnt_proteins = 0
    cnt_updates = 0
    for protein_acc, protein_info in proteins:
        cnt_proteins += 1
        if not cnt_proteins % 1000000:
            logger.info(f"{cnt_proteins:>12,}")

        taxon_id = protein_info["taxon"]
        protein_counts = {"all": 1, "entries": {}, "databases": {}}
        protein_entries = {}
        protein_matches = matches.get(protein_acc, {})
        protein_sets = set()
        for entry_acc in protein_matches:
            database = entries[entry_acc]["database"]

            try:
                protein_entries[database].add(entry_acc)
            except KeyError:
                protein_entries[database] = {entry_acc}

            try:
                set_acc = entry_set[entry_acc]
            except KeyError:
                pass
            else:
                protein_sets.add(set_acc)

            protein_counts["databases"][database] = 1
            protein_counts["entries"][entry_acc] = 1

        dom_arch.update(protein_matches)
        upid = proteomes.get(protein_acc)

        xrefs = {
            "domain_architectures": {dom_arch.identifier},
            "entries": protein_entries,
            "proteins": protein_counts,
            "proteomes": {upid} if upid else set(),
            "sets": protein_sets,
            "structures": structures.get(protein_acc, set())
        }

        for tax_id in lineages[taxon_id]:
            taxa.update(tax_id, xrefs, replace=False)
        # taxa.update(taxon_id, xrefs, replace=False)

        cnt_updates += 1
        if not cnt_updates % sync_frequency:
            taxa.sync()

    proteins.close()
    proteomes.close()
    matches.close()
    logger.info(f"{cnt_proteins:>12}")

    for tax_id in lineages:
        xrefs = {
            "domain_architectures": set(),
            "entries": {},
            "proteins": {"all": 0, "databases": {}, "entries": {}},
            "proteomes": set(),
            "sets": set(),
            "structures": set()
        }
        taxa.update(tax_id, xrefs, replace=False)

    return taxa


def update_counts(url: str, src_proteins: str, src_proteomes:str,
                  src_matches: str, processes: int=1,
                  sync_frequency: int=1000000, tmpdir: Optional[str]=None):
    taxa = _export(url, src_proteins, src_proteomes, src_matches,
                   sync_frequency, tmpdir)
    _ = input("press any key ")
    size = taxa.merge(processes=processes)
    logger.info(f"disk usage: {size/1024/1024:.0f} MB")
    print(taxa.filepath)
    _ = input("press any key ")

    logger.info("updating MySQL tables")
    con = MySQLdb.connect(**parse_url(url), charset="utf8")
    table1 = Table(con, query="UPDATE webfront_taxonomy SET counts = %s "
                              "WHERE accession = %s")
    table2 = Table(con, query="INSERT INTO webfront_taxonomy_database "
                              "(tax_id, source_database, counts) "
                              "VALUES (%s, %s, %s)")
    table3 = Table(con, query="INSERT INTO webfront_taxonomy_entry  "
                              "(tax_id, entry_acc, counts)  "
                              "VALUES (%s, %s, %s)")

    for tax_id, xrefs in taxa:
        protein_counts = xrefs.pop("proteins")

        # Counts for `webfront_taxonomy`
        counts = reduce(xrefs)
        # All proteins (not filtered by database/entry)
        counts["proteins"] = protein_counts["all"]
        # Total number of entries (all databases)
        counts["entries"]["total"] = sum(counts["entries"].values())
        table1.update((json.dumps(counts), tax_id))

        # Remove elements we do not need for other tables
        del counts["entries"]
        del counts["sets"]

        # Counts for `webfront_taxonomy_database`
        for database, cnt in protein_counts["databases"].items():
            counts["proteins"] = cnt
            table2.insert((tax_id, database, json.dumps(counts)))

        # Counts for `webfront_taxonomy_database`
        for entry_acc, cnt in protein_counts["entries"].items():
            counts["proteins"] = cnt
            table3.insert((tax_id, entry_acc, json.dumps(counts)))

    for t in (table1, table2, table3):
        t.close()

    con.commit()
    con.close()
    taxa.close()
    logger.info("complete")


    # logger.info("creating taxonomy database")
    # fd, database = mkstemp(dir=tmpdir, suffix=".db")
    # os.close(fd)
    # os.remove(database)
    #
    # with io.KVdb(database, insertonly=True) as kvdb:
    #     for tax_id, xrefs in taxa:
    #         kvdb[tax_id] = xrefs
    #
    # taxa.close()  # delete temporary Store
    #
    # with io.KVdb(database, writeback=True) as kvdb:
    #     logger.info("propagating to lineage")
    #     cnt_taxa = 0
    #     for taxon in iter_taxa(url, lineage=True):
    #         node = kvdb[taxon["id"]]
    #
    #         for tax_id in taxon["lineage"]:
    #             if tax_id == taxon["id"]:
    #                 continue
    #
    #             _node = kvdb[tax_id]
    #             _node["domain_architectures"] |= node["domain_architectures"]
    #             _node["proteomes"] |= node["proteomes"]
    #             _node["structures"] |= node["structures"]
    #             _node["sets"] |= node["sets"]
    #
    #             for database, accessions in node["entries"].items():
    #                 try:
    #                     _node["entries"][database] |= accessions
    #                 except KeyError:
    #                     _node["entries"][database] = accessions.copy()
    #
    #             kvdb[tax_id] = _node
    #
    #         cnt_taxa += 1
    #         if not cnt_taxa % sync_frequency:
    #             kvdb.sync()
    #             logger.info(f"{cnt_taxa:>8,}")
    #
    #     logger.info(f"{cnt_taxa:>8,}")
    #     logger.info("updating MySQL tables")
    #     con = MySQLdb.connect(**parse_url(url), charset="utf8")
    #     table1 = Table(con,
    #                    query="UPDATE webfront_taxonomy "
    #                          "SET counts = %s "
    #                          "WHERE accession = %s")
    #     table2 = Table(con,
    #                    query="INSERT INTO webfront_taxonomy_database "
    #                          "(tax_id, source_database, counts) "
    #                          "VALUES (%s, %s, %s)")
    #     table3 = Table(con,
    #                    query="INSERT INTO webfront_taxonomy_entry "
    #                          "(tax_id, entry_acc, counts) "
    #                          "VALUES (%s, %s, %s)")
    #
    #     for tax_id, xrefs in kvdb:
    #         protein_counts = xrefs.pop("proteins")
    #
    #         # Counts for `webfront_taxonomy`
    #         counts = reduce(xrefs)
    #         # All proteins (not filtered by database/entry)
    #         counts["proteins"] = protein_counts["all"]
    #         # Total number of entries (all databases)
    #         counts["entries"]["total"] = sum(counts["entries"].values())
    #         table1.update((json.dumps(counts), tax_id))
    #
    #         # Remove elements we do not need anymore
    #         del counts["entries"]
    #         del counts["sets"]
    #
    #         # Counts for `webfront_taxonomy_database`
    #         for database, cnt in protein_counts["databases"].items():
    #             counts["proteins"] = cnt
    #             table2.insert((tax_id, database, json.dumps(counts)))
    #
    #         # Counts for `webfront_taxonomy_database`
    #         for entry_acc, cnt in protein_counts["entries"].items():
    #             counts["proteins"] = cnt
    #             table3.insert((tax_id, entry_acc, json.dumps(counts)))
    #
    #     for t in (table1, table2, table3):
    #         t.close()
    #
    #     con.commit()
    #     con.close()
    #     size += kvdb.size
    #
    # os.remove(database)
    # logger.info(f"disk usage: {size/1024/1024:.0f} MB")
