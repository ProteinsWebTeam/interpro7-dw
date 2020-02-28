# -*- coding: utf-8 -*-

from typing import Optional

import MySQLdb

from interpro7dw import logger
from interpro7dw.ebi.interpro.utils import Table
from interpro7dw.utils import Store, dataload, url2dict
from .utils import jsonify, reduce


def insert_proteomes(entries: str, proteomes: str, structures: str,
                     proteins: str, uniprot2ida: str, uniprot2matches: str,
                     uniprot2proteome: str, stg_url: str,
                     dir: Optional[str]=None, processes: int=1):
    logger.info("preparing data")
    entries = dataload(entries)
    proteomes = dataload(proteomes)
    uniprot2pdbe = {}
    for pdb_id, entry in dataload(structures).items():
        for uniprot_acc in entry["proteins"]:
            try:
                uniprot2pdbe[uniprot_acc].add(pdb_id)
            except KeyError:
                uniprot2pdbe[uniprot_acc] = {pdb_id}

    with Store(None, Store.chunk(proteomes.keys(), 100), dir) as store:
        # Init all proteomes
        for proteome_id in proteomes:
            store.update(proteome_id, {
                "domain_architectures": set(),
                "entries": {},
                "proteins": 0,
                "sets": set(),
                "structures": set(),
                "taxa": set()
            }, replace=False)

        proteins = Store(proteins)
        u2ida = Store(uniprot2ida)
        u2matches = Store(uniprot2matches)
        u2proteome = Store(uniprot2proteome)

        logger.info("starting")
        i =0
        for uniprot_acc, info in proteins.items():
            try:
                proteome_id = u2proteome[uniprot_acc]
            except KeyError:
                continue

            try:
                dom_arch, dom_arch_id = u2ida[uniprot_acc]
            except KeyError:
                ida_obj = set()
            else:
                ida_obj = {dom_arch_id}

            databases = {}
            clans = set()
            for entry_acc in u2matches.get(uniprot_acc, []):
                entry = entries[entry_acc]

                try:
                    databases[entry.database].add(entry_acc)
                except KeyError:
                    databases[entry.database] = {entry_acc}

                if entry.clan:
                    clans.add(entry.clan["accession"])

            store.update(proteome_id, {
                "domain_architectures": ida_obj,
                "entries": databases,
                "proteins": 1,
                "sets": clans,
                "structures": uniprot2pdbe.get(uniprot_acc, set()),
                "taxa": {info["taxid"]}
            }, replace=False)

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
        cur.execute("DROP TABLE IF EXISTS webfront_proteome")
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
                counts LONGTEXT DEFAULT NULL
            ) CHARSET=utf8 DEFAULT COLLATE=utf8_unicode_ci
            """
        )
        cur.close()

        sql = """
            INSERT INTO webfront_proteome VALUES (%s, %s, %s, %s, %s, %s) 
        """
        with Table(con, sql) as table:
            for proteome_id, xrefs in store.items():
                info = proteomes[proteome_id]

                table.insert((
                    proteome_id,
                    info["name"],
                    1 if info["is_reference"] else 0,
                    info["strain"],
                    info["assembly"],
                    info["taxon_id"],
                    jsonify(reduce(xrefs))
                ))

        con.commit()
        con.close()

    logger.info("complete")
