# -*- coding: utf-8 -*-

import os
import json
from typing import Generator, Optional

import MySQLdb
import MySQLdb.cursors

from i7dw import io, logger, uniprot
from i7dw.interpro import Table, mysql
from .utils import parse_url


def insert_proteomes(my_url: str, ora_url: str):
    query = """
        INSERT INTO webfront_proteome (accession, name, is_reference, strain,
                                       assembly, taxonomy_id)
        VALUES (%s, %s, %s, %s, %s, %s)
    """

    taxa = {tax["id"] for tax in mysql.taxonomy.iter_taxa(my_url)}

    con = MySQLdb.connect(**parse_url(my_url), charset="utf8")
    with Table(con, query) as table:
        for p in uniprot.get_proteomes(ora_url):
            if p["tax_id"] in taxa:
                table.insert((
                    p["id"],
                    p["name"],
                    1 if p["is_reference"] else 0,
                    p["strain"],
                    p["assembly"],
                    p["tax_id"]
                ))
            else:
                # INTERPRO.ETAXI might be out-of-date
                logger.warning(f"missing taxon (ID: {p['tax_id']})")

    con.commit()
    con.close()


def iter_proteomes(url: str) -> Generator[dict, None, None]:
    con = MySQLdb.connect(**parse_url(url), charset="utf8")
    cur = MySQLdb.cursors.SSCursor(con)
    cur.execute(
        """
        SELECT accession, name, is_reference, strain, assembly, taxonomy_id
        FROM webfront_proteome
        """
    )

    for row in cur:
        yield {
            "accession": row[0],
            "name": row[1],
            "is_reference": bool(row[2]),
            "strain": row[3],
            "assembly": row[4],
            "taxon": row[5]
        }

    cur.close()
    con.close()


def update_counts(url: str, src_proteins: str, src_proteomes:str,
                  src_matches: str, processes: int=1,
                  tmpdir: Optional[str]=None, sync_frequency: int=100000):
    # Get required MySQL data
    logger.info("loading data")
    entries = mysql.entries.get_entries(url)
    entry_set = {}
    for s in mysql.entries.iter_sets(url):
        for entry_acc in s["members"]:
            entry_set[entry_acc] = s["accession"]

    proteomes = {p["accession"]: p for p in iter_proteomes(url)}
    structures = {}
    for s in mysql.structures.iter_structures(url):
        pdbe_id = s["accession"]
        for protein_acc, chains in s["proteins"].items():
            try:
                structures[protein_acc][pdbe_id] = chains
            except KeyError:
                structures[protein_acc] = {pdbe_id: chains}

    # Open existing stores containing protein-related info
    proteins = io.Store(src_proteins)
    protein_proteome = io.Store(src_proteomes)
    protein_matches = io.Store(src_matches)

    protein_counts = {}
    cnt_proteins = 0
    cnt_updates = 0
    with io.Store(keys=io.Store.chunk_keys(proteomes, 100), tmpdir=tmpdir) as xrefs:
        for protein_acc, protein_info in proteins:
            cnt_proteins += 1
            if not cnt_proteins % 10000000:
                logger.info(f"{cnt_proteins:>12}")

            tax_id = protein_info["taxon"]
            try:
                upid = protein_proteome[protein_acc]
            except KeyError:
                continue

            # TODO: resume refactoring there
            protein_entries = set(protein_matches.get(protein_acc, []))

            entry_databases = {}
            entry_sets = set()
            for entry_acc in prot_entries:
                database = entries[entry_acc]["database"]
                if database in entry_databases:
                    entry_databases[database].add(entry_acc)
                else:
                    entry_databases[database] = {entry_acc}

                try:
                    set_acc = entry2set[entry_acc]
                except KeyError:
                    pass
                else:
                    entry_sets.add(set_acc)

            _xrefs = {
                "domain_architectures": set(),
                "entries": entry_databases,
                "sets": entry_sets,
                "structures": set(),
                "taxa": {tax_id}
            }
            try:
                ida, ida_id = protein2ida[protein_acc]
            except KeyError:
                pass
            else:
                _xrefs["domain_architectures"] = {ida}

            try:
                pdb_ids = protein2structures[protein_acc]
            except KeyError:
                pass
            else:
                _xrefs["structures"] = pdb_ids

            xrefs.update(upid, _xrefs)
            cnt_updates += 1
            if not cnt_updates % sync_frequency:
                xrefs.sync()

            if upid in protein_counts:
                protein_counts[upid] += 1
            else:
                protein_counts[upid] = 1

        proteins.close()
        protein2proteome.close()
        protein2matches.close()
        protein2ida.close()
        logger.info(f"{cnt_proteins:>12}")

        for upid, cnt in protein_counts.items():
            proteomes.pop(upid)
            xrefs.update(upid, {"proteins": cnt})

        # Remaining proteomes
        for upid in proteomes:
            xrefs.update(upid, {
                "domain_architectures": set(),
                "entries": {},
                "proteins": set(),
                "sets": set(),
                "structures": set(),
                "taxa": {proteomes[upid]["taxon"]},
            })

        entries = None
        protein2structures = None
        entry2set = None
        proteomes = None
        protein_counts = None

        size = xrefs.merge(processes=processes)
        logger.info("Disk usage: {:.0f}MB".format(size/1024**2))

        con = MySQLdb.connect(**parse_url(my_url), charset="utf8")
        query = "UPDATE webfront_proteome SET counts = %s WHERE accession = %s"
        with Table(con, query) as table:
            for upid, _xrefs in xrefs:
                counts = reduce(_xrefs)
                counts["entries"]["total"] = sum(counts["entries"].values())
                table.update((json.dumps(counts), upid))

        con.commit()
        con.close()

    logger.info("complete")
