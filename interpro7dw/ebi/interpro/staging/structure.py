# -*- coding: utf-8 -*-

from typing import Optional

import MySQLdb

from interpro7dw import logger
from interpro7dw.ebi.interpro.utils import Table, overlaps_pdb_chain
from interpro7dw.utils import Store, dataload, url2dict
from .utils import jsonify, reduce


def insert_structures(entries: str, structures: str, proteins: str,
                      uniprot2ida: str, uniprot2matches: str,
                      uniprot2proteome: str, stg_url: str,
                      dir: Optional[str]=None, processes: int=1):
    logger.info("preparing data")
    entries = dataload(entries)
    uniprot2pdbe = {}
    pdb_ids = []
    for pdb_id, entry in dataload(structures).items():
        pdb_ids.append(pdb_id)
        for uniprot_acc, chains in entry["proteins"].items():
            try:
                uniprot2pdbe[uniprot_acc][pdb_id] = chains
            except KeyError:
                uniprot2pdbe[uniprot_acc] = {pdb_id: chains}

    with Store(None, Store.chunk(pdb_ids, 10000), dir) as store:
        # Init all structures
        for pdb_id in pdb_ids:
            store.update(pdb_id, {
                "domain_architectures": set(),
                "entries": {},
                "proteomes": set(),
                "proteins": 0,
                "sets": set(),
                "structures": set()
            }, replace=False)

        proteins = Store(proteins)
        u2ida = Store(uniprot2ida)
        u2matches = Store(uniprot2matches)
        u2proteome = Store(uniprot2proteome)

        logger.info("starting")
        i =0
        for uniprot_acc, info in proteins.items():
            try:
                pdb_entries = uniprot2pdbe[uniprot_acc]
            except KeyError:
                continue  # UniProt entry has no PDBe structures

            try:
                dom_arch, dom_arch_id = u2ida[uniprot_acc]
            except KeyError:
                ida_obj = set()
            else:
                ida_obj = {dom_arch_id}

            try:
                proteome_obj = {u2proteome[uniprot_acc]}
            except KeyError:
                proteome_obj = set()

            xrefs = {
                "domain_architectures": ida_obj,
                "proteomes": proteome_obj,
                "proteins": 1,
                "taxa": {info["taxid"]}
            }

            matches = u2matches.get(uniprot_acc, {})
            for pdb_id, chains in pdb_entries.items():
                databases = {}
                clans = set()

                for entry_acc, locations in matches.items():
                    entry = entries[entry_acc]

                    for chain_id, segments in chains.items():
                        if overlaps_pdb_chain(locations, segments):
                            try:
                                databases[entry.database].add(entry_acc)
                            except KeyError:
                                databases[entry.database] = {entry_acc}

                            if entry.clan:
                                clans.add(entry.clan["accession"])

                            break  # Skip other chains: next entry

                _xrefs = xrefs.copy()
                _xrefs["entries"] = databases
                _xrefs["sets"] = clans
                store.update(pdb_id, _xrefs, replace=False)

            i += 1
            if not i % 1000000:
                store.sync()

                if not i % 10000000:
                    logger.info(f"{i:>12,}")

        logger.info(f"{i:>12,}")

        proteins.close()
        u2ida.close()
        u2matches.close()
        u2proteome.close()

        size = store.merge(processes=processes)
        logger.info(f"temporary files: {size/1024/1024:.0f} MB")

        con = MySQLdb.connect(**url2dict(stg_url))
        cur = con.cursor()
        cur.execute("DROP TABLE IF EXISTS webfront_structure")
        cur.execute(
            """
            CREATE TABLE webfront_structure
            (
                accession VARCHAR(4) PRIMARY KEY NOT NULL,
                name VARCHAR(512) NOT NULL,
                source_database VARCHAR(10) NOT NULL,
                experiment_type VARCHAR(16) NOT NULL,
                release_date DATETIME NOT NULL,
                resolution FLOAT DEFAULT NULL,
                literature LONGTEXT,
                chains LONGTEXT NOT NULL,
                proteins LONGTEXT NOT NULL,
                secondary_structures LONGTEXT NOT NULL,
                counts LONGTEXT DEFAULT NULL
            ) CHARSET=utf8 DEFAULT COLLATE=utf8_unicode_ci
            """
        )
        cur.close()

        structures = dataload(structures)

        sql = """
            INSERT INTO webfront_structure 
            VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s) 
        """
        with Table(con, sql) as table:
            for pdb_id, xrefs in store.items():
                info = structures[pdb_id]
                table.insert((
                    pdb_id,
                    info["name"],
                    "pdb",
                    info["evidence"],
                    info["date"],
                    info["resolution"],
                    jsonify(info["citations"]),
                    # Sorted list of unique chain (e.g. 'A', 'B', ...)
                    jsonify(sorted({chain_id
                                    for chains in info["proteins"].values()
                                    for chain_id in chains})),
                    jsonify(info["proteins"]),
                    jsonify(info["secondary_structures"]),
                    jsonify(reduce(xrefs))
                ))

        con.commit()
        con.close()

    logger.info("complete")
