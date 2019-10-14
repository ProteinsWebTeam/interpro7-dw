# -*- coding: utf-8 -*-

import json
import os
from typing import Dict, Generator, Optional

import MySQLdb
import MySQLdb.cursors

from i7dw import io, logger, pdbe
from i7dw.interpro import DomainArchitecture, Table, mysql
from .utils import parse_url, reduce


def insert_structures(my_url: str, ora_url: str):
    query = """
        INSERT INTO webfront_structure (accession, name, source_database, 
                                        experiment_type, release_date, 
                                        resolution, literature, chains, 
                                        proteins, secondary_structures) 
        VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
    """

    sec_structures = pdbe.get_secondary_structures(ora_url)

    con = MySQLdb.connect(**parse_url(my_url), charset="utf8")
    with Table(con, query) as table:
        for s in pdbe.get_structures(ora_url):
            pdbe_id = s["id"]

            chains = set()
            for protein_chains in s["proteins"].values():
                for chain_id in protein_chains:
                    chains.add(chain_id)

            table.insert((
                pdbe_id, s["name"], "pdb", s["evidence"], s["date"],
                s["resolution"], json.dumps(s["citations"]),
                json.dumps(sorted(chains)), json.dumps(s["proteins"]),
                json.dumps(sec_structures.get(pdbe_id, []))
            ))

    con.commit()
    con.close()


def iter_structures(url: str) -> Generator[Dict, None, None]:
    con = MySQLdb.connect(**parse_url(url), charset="utf8")
    cur = MySQLdb.cursors.SSCursor(con)
    cur.execute(
        """
        SELECT accession, name, experiment_type, release_date, 
               resolution, literature, proteins
        FROM webfront_structure
        """
    )

    for row in cur:
        yield {
            "accession": row[0],
            "name": row[1],
            "evidence": row[2],
            "date": row[3],
            "resolution": row[4],
            "citations": json.loads(row[5]),
            "proteins": json.loads(row[6])
        }

    cur.close()
    con.close()


def update_counts(my_url: str, src_proteins: str, src_proteomes: str,
                  src_matches: str, processes: int=1,
                  sync_frequency: int=100000, tmpdir: Optional[str]=None):
    if tmpdir:
        os.makedirs(tmpdir, exist_ok=True)

    # Get required MySQL data
    logger.info("loading data")
    entries = mysql.entries.get_entries(my_url)
    entry_set = {}
    for s in mysql.entries.iter_sets(my_url):
        set_acc = s["accession"]
        for entry_acc in s["members"]:
            entry_set[entry_acc] = set_acc

    pdbe_ids = set()
    structures = {}
    for s in iter_structures(my_url):
        pdbe_id = s["accession"]
        pdbe_ids.add(pdbe_id)
        for protein_acc, chains in s["proteins"].items():
            try:
                structures[protein_acc][pdbe_id] = chains
            except KeyError:
                structures[protein_acc] = {pdbe_id: chains}

    dom_arch = DomainArchitecture(entries)

    # Open existing stores containing protein-related info
    proteins = io.Store(src_proteins)
    proteomes = io.Store(src_proteomes)
    matches = io.Store(src_matches)

    protein_counts = {}
    cnt_proteins = 0
    cnt_updates = 0
    with io.Store(keys=io.Store.chunk_keys(pdbe_ids, 100), tmpdir=tmpdir) as store:
        for protein_acc, protein_info in proteins:
            cnt_proteins += 1
            if not cnt_proteins % 10000000:
                logger.info(f"{cnt_proteins:>12}")

            try:
                protein_structures = structures[protein_acc]
            except KeyError:
                continue

            protein_matches = matches.get(protein_acc, {})
            upid = proteomes.get(protein_acc)

            dom_arch.update(protein_matches)

            base_xrefs = {
                "domain_architectures": {dom_arch.identifier},
                "proteomes": {upid} if upid else set(),
                "taxa": {protein_info["taxon"]}
            }

            for pdbe_id, chains in protein_structures.items():
                protein_entries = {}
                for entry_acc, locations in protein_matches.items():
                    database = entries[entry_acc]["database"]

                    if mysql.entries.overlaps_with_structure(locations, chains):
                        try:
                            protein_entries[database].add(entry_acc)
                        except KeyError:
                            protein_entries[database] = {entry_acc}

                # Adding sets for overlapping entries
                protein_sets = set()
                for accessions in protein_entries.values():
                    for entry_acc in accessions:
                        try:
                            set_acc = entry_set[entry_acc]
                        except KeyError:
                            pass
                        else:
                            protein_sets.add(set_acc)

                xrefs = base_xrefs.copy()
                xrefs["entries"] = protein_entries
                xrefs["sets"] = protein_sets
                store.update(pdbe_id, xrefs)

                try:
                    protein_counts[pdbe_id] += 1
                except KeyError:
                    protein_counts[pdbe_id] = 1

            cnt_updates += 1
            if not cnt_updates % sync_frequency:
                store.sync()

        proteins.close()
        proteomes.close()
        matches.close()
        logger.info(f"{cnt_proteins:>12}")

        for pdbe_id in pdbe_ids:
            xrefs = {
                "domain_architectures": set(),
                "entries": {},
                "proteomes": set(),
                "proteins": 0,
                "sets": set(),
                "taxa": set()
            }

            try:
                xrefs["proteins"] = protein_counts[pdbe_id]
            except KeyError:
                pass
            finally:
                store.update(pdbe_id, xrefs)

        size = store.merge(processes=processes)
        logger.info(f"disk usage: {size/1024/1024:.0f}MB")

        con = MySQLdb.connect(**parse_url(my_url), charset="utf8")
        query = "UPDATE webfront_structure SET counts = %s WHERE accession = %s"
        with Table(con, query) as table:
            for pdbe_id, xrefs in store:
                counts = reduce(xrefs)
                counts["entries"]["total"] = sum(counts["entries"].values())
                table.update((json.dumps(counts), pdbe_id))

        con.commit()
        con.close()

    logger.info("complete")
