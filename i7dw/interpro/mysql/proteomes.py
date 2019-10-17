# -*- coding: utf-8 -*-

import os
import json
from typing import Generator, Optional

import MySQLdb
import MySQLdb.cursors

from i7dw import io, logger, uniprot
from i7dw.interpro import DomainArchitecture, Table, mysql
from .utils import parse_url, reduce


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
    dom_arch = DomainArchitecture(entries)

    entry_set = {}
    for s in mysql.entries.iter_sets(url):
        for entry_acc in s["members"]:
            entry_set[entry_acc] = s["accession"]

    proteomes = [p["accession"] for p in iter_proteomes(url)]
    structures = {}
    for s in mysql.structures.iter_structures(url):
        pdbe_id = s["accession"]
        for protein_acc, chains in s["proteins"].items():
            try:
                structures[protein_acc].add(pdbe_id)
            except KeyError:
                structures[protein_acc] = {pdbe_id}

    logger.info("starting")
    tmp_proteomes = io.mktemp(dir=tmpdir)
    with io.Store(tmp_proteomes, keys=io.Store.chunk_keys(proteomes, 100),
                  tmpdir=tmpdir) as store:
        # Init all proteomes
        for upid in proteomes:
            store.update(upid, {
                "domain_architectures": set(),
                "entries": {},
                "proteins": 0,
                "sets": set(),
                "structures": set(),
                "taxa": set()
            }, replace=False)

        proteins = io.Store(src_proteins)
        proteomes = io.Store(src_proteomes)
        matches = io.Store(src_matches)

        cnt_proteins = 0
        for protein_acc, protein_info in proteins:
            try:
                upid = proteomes[protein_acc]
            except KeyError:
                continue

            protein_matches = matches.get(protein_acc, {})
            protein_entries = {}
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

            dom_arch.update(protein_matches)
            store.update(upid, {
                "domain_architectures": {dom_arch.identifier},
                "entries": protein_entries,
                "proteins": 1,
                "sets": protein_sets,
                "structures": structures.get(protein_acc, set()),
                "taxa": {protein_info["taxon"]}
            }, replace=False)

            cnt_proteins += 1
            if not cnt_proteins % 10000000:
                logger.info(f"{cnt_proteins:>12}")

            if not cnt_proteins % sync_frequency:
                store.sync()

        proteins.close()
        proteomes.close()
        matches.close()
        logger.info(f"{cnt_proteins:>12}")

        entries = None
        entry_set = None
        dom_arch = None
        structures = None

        size = store.merge(processes=processes)
        size += os.path.getsize(tmp_proteomes)

        con = MySQLdb.connect(**parse_url(url), charset="utf8")
        query = "UPDATE webfront_proteome SET counts = %s WHERE accession = %s"
        with Table(con, query) as table:
            for upid, xrefs in store:
                counts = reduce(xrefs)
                counts["entries"]["total"] = sum(counts["entries"].values())
                table.update((json.dumps(counts), upid))

        con.commit()
        con.close()

    os.remove(tmp_proteomes)
    logger.info(f"disk usage: {size/1024/1024:.0f} MB")
