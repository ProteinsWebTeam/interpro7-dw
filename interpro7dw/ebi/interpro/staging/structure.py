# -*- coding: utf-8 -*-

import MySQLdb

from interpro7dw import logger
from interpro7dw.ebi.interpro.utils import Table, overlaps_pdb_chain
from interpro7dw.utils import Store, dataload, url2dict
from .utils import jsonify, reduce


def insert_structures(p_entries: str, p_proteins: str, p_structures: str,
                      p_uniprot2ida: str, p_uniprot2matches: str,
                      p_uniprot2proteome: str, stg_url: str):
    logger.info("preparing data")
    entries = {}
    for entry in dataload(p_entries).values():
        entries[entry.accession] = (entry.database, entry.clan)

    uniprot2pdbe = {}
    xrefs = {}
    for pdb_id, entry in dataload(p_structures).items():
        for uniprot_acc, chains in entry["proteins"].items():
            try:
                uniprot2pdbe[uniprot_acc][pdb_id] = chains
            except KeyError:
                uniprot2pdbe[uniprot_acc] = {pdb_id: chains}

        xrefs[pdb_id] = {
            "domain_architectures": set(),
            "entries": {},
            "proteomes": set(),
            "proteins": 0,
            "sets": set(),
            "taxa": set()
        }

    proteins = Store(p_proteins)
    u2ida = Store(p_uniprot2ida)
    u2matches = Store(p_uniprot2matches)
    u2proteome = Store(p_uniprot2proteome)

    logger.info("starting")
    i = 0
    for uniprot_acc in sorted(uniprot2pdbe):
        info = proteins[uniprot_acc]

        try:
            dom_arch, dom_arch_id = u2ida[uniprot_acc]
        except KeyError:
            dom_arch_id = None

        proteome_id = u2proteome.get(uniprot_acc)
        matches = u2matches.get(uniprot_acc, {})

        for pdb_id, chains in uniprot2pdbe[uniprot_acc].items():
            _xrefs = xrefs[pdb_id]

            if dom_arch_id:
                _xrefs["domain_architectures"].add(dom_arch_id)

            if proteome_id:
                _xrefs["proteomes"].add(proteome_id)

            _xrefs["proteins"] += 1
            _xrefs["taxa"].add(info["taxid"])

            for entry_acc, locations in matches.items():
                database, clan = entries[entry_acc]

                for chain_id, segments in chains.items():
                    if overlaps_pdb_chain(locations, segments):
                        try:
                            _xrefs["entries"][database].add(entry_acc)
                        except KeyError:
                            _xrefs["entries"][database] = {entry_acc}

                        if clan:
                            _xrefs["sets"].add(clan["accession"])

                        break  # Skip other chains

        i += 1
        if not i % 10000:
            logger.info(f"{i:>12,}")

    logger.info(f"{i:>12,}")

    proteins.close()
    u2ida.close()
    u2matches.close()
    u2proteome.close()

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
            resolution FLOAT,
            literature LONGTEXT,
            chains LONGTEXT NOT NULL,
            proteins LONGTEXT NOT NULL,
            secondary_structures LONGTEXT,
            counts LONGTEXT NOT NULL
        ) CHARSET=utf8 DEFAULT COLLATE=utf8_unicode_ci
        """
    )
    cur.close()

    sql = """
        INSERT INTO webfront_structure 
        VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s) 
    """
    with Table(con, sql) as table:
        for pdb_id, info in dataload(p_structures).items():
            counts = reduce(xrefs[pdb_id])
            counts["entries"]["total"] = sum(counts["entries"].values())
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
                jsonify(counts)
            ))

    con.commit()
    con.close()

    logger.info("complete")
