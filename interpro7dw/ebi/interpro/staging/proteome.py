# -*- coding: utf-8 -*-

import MySQLdb

from interpro7dw import logger
from interpro7dw.ebi.interpro.utils import Table
from interpro7dw.utils import Store, loadobj, url2dict
from .utils import jsonify, reduce


def insert_proteomes(p_proteomes: str, p_structures: str, p_proteins: str,
                     p_uniprot2ida: str, p_uniprot2entries: str,
                     p_uniprot2proteome: str, stg_url: str):
    logger.info("preparing data")
    proteomes = loadobj(p_proteomes)
    uniprot2pdbe = {}
    for pdb_id, entry in loadobj(p_structures).items():
        for uniprot_acc in entry["proteins"]:
            try:
                uniprot2pdbe[uniprot_acc].add(pdb_id)
            except KeyError:
                uniprot2pdbe[uniprot_acc] = {pdb_id}

    # Init all proteomes
    xrefs = {}
    for proteome_id in proteomes:
        xrefs[proteome_id] = {
            "domain_architectures": set(),
            "entries": {},
            "proteins": 0,
            "sets": set(),
            "structures": set(),
            "taxa": set()
        }

    proteins = Store(p_proteins)
    u2entries = Store(p_uniprot2entries)
    u2ida = Store(p_uniprot2ida)
    u2proteome = Store(p_uniprot2proteome)

    logger.info("starting")
    i = 0
    for uniprot_acc, proteome_id in u2proteome.items():
        proteome = xrefs[proteome_id]
        proteome["proteins"] += 1

        info = proteins[uniprot_acc]
        proteome["taxa"].add(info["taxid"])

        try:
            dom_members, dom_arch, dom_arch_id = u2ida[uniprot_acc]
        except KeyError:
            pass
        else:
            proteome["domain_architectures"].add(dom_arch_id)

        entries = u2entries.get(uniprot_acc, [])
        for entry_acc, database, clan_acc, go_terms in entries:
            try:
                proteome["entries"][database].add(entry_acc)
            except KeyError:
                proteome["entries"][database] = {entry_acc}

            if clan_acc:
                proteome["sets"].add(clan_acc)

        try:
            pdb_ids = uniprot2pdbe[uniprot_acc]
        except KeyError:
            pass
        else:
            proteome["structures"] |= pdb_ids

        i += 1
        if not i % 10000000:
            logger.info(f"{i:>12,}")

    logger.info(f"{i:>12,}")

    proteins.close()
    u2ida.close()
    u2entries.close()
    u2proteome.close()

    con = MySQLdb.connect(**url2dict(stg_url), charset="utf8mb4")
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
            counts LONGTEXT NOT NULL
        ) CHARSET=utf8mb4 DEFAULT COLLATE=utf8mb4_unicode_ci
        """
    )
    cur.close()

    sql = """
        INSERT INTO webfront_proteome VALUES (%s, %s, %s, %s, %s, %s, %s) 
    """
    with Table(con, sql) as table:
        for proteome_id, info in proteomes.items():
            counts = reduce(xrefs[proteome_id])
            counts["entries"]["total"] = sum(counts["entries"].values())
            table.insert((
                proteome_id,
                info["name"],
                1 if info["is_reference"] else 0,
                info["strain"],
                info["assembly"],
                info["taxon_id"],
                jsonify(counts)
            ))

    con.commit()
    con.close()

    logger.info("complete")
