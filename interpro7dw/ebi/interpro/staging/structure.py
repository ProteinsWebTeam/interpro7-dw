from typing import Dict, List

import cx_Oracle
import MySQLdb

from interpro7dw import logger
from interpro7dw.ebi.interpro.utils import Table, blob_as_str
from interpro7dw.ebi.interpro.utils import overlaps_pdb_chain
from interpro7dw.utils import DumpFile, Store, loadobj, url2dict
from .utils import jsonify, reduce


def insert_structures(p_entries: str, p_proteins: str, p_structures: str,
                      p_uniprot2ida: str, p_uniprot2matches: str,
                      p_uniprot2proteome: str, stg_url: str):
    logger.info("preparing data")
    entries = {}
    for entry in loadobj(p_entries).values():
        entries[entry.accession] = (entry.database, entry.clan)

    uniprot2pdbe = {}
    xrefs = {}
    for pdb_id, entry in loadobj(p_structures).items():
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
            dom_members, dom_arch, dom_arch_id = u2ida[uniprot_acc]
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

    con = MySQLdb.connect(**url2dict(stg_url), charset="utf8mb4")
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
        ) CHARSET=utf8mb4 DEFAULT COLLATE=utf8mb4_unicode_ci
        """
    )
    cur.close()

    sql = """
        INSERT INTO webfront_structure 
        VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s) 
    """
    with Table(con, sql) as table:
        for pdb_id, info in loadobj(p_structures).items():
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
                                for chain_id in chains}), nullable=False),
                jsonify(info["proteins"], nullable=False),
                jsonify(info["secondary_structures"]),
                jsonify(counts)
            ))

    con.commit()
    con.close()

    logger.info("complete")


def insert_structural_models(pro_url: str, stg_url: str, p_entries: str):
    entries = loadobj(p_entries)

    my_con = MySQLdb.connect(**url2dict(stg_url))
    my_cur = my_con.cursor()
    my_cur.execute("DROP TABLE IF EXISTS webfront_structuralmodel")
    my_cur.execute(
        """
        CREATE TABLE webfront_structuralmodel
        (
            model_id INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
            accession VARCHAR(25) NOT NULL,
            algorithm VARCHAR(20) NOT NULL,
            alignment LONGBLOB NOT NULL,
            contacts LONGBLOB NOT NULL,
            plddt LONGBLOB NOT NULL,
            structure LONGBLOB NOT NULL
        ) CHARSET=utf8mb4 DEFAULT COLLATE=utf8mb4_unicode_ci
        """
    )

    logger.info("inserting structure predictions")
    ora_con = cx_Oracle.connect(pro_url)
    ora_cur = ora_con.cursor()
    ora_cur.outputtypehandler = blob_as_str
    ora_cur.execute(
        """
        SELECT METHOD_AC, ALGORITHM, ALIGNMENTS, CONTACTS, PLDDT, STRUCTURE
        FROM INTERPRO.STRUCT_MODEL
        """
    )

    req = """
        INSERT INTO webfront_structuralmodel (
          accession, algorithm, alignment, contacts, plddt, structure
        )
        VALUES (%s, %s, %s, %s, %s, %s)
    """

    for entry_acc, algorithm, msa_gz, cmap_gz, plddt_gz, pdb_gz in ora_cur:
        try:
            entry = entries[entry_acc]
        except KeyError:
            continue

        my_cur.execute(req, (entry_acc, algorithm, msa_gz, cmap_gz, plddt_gz,
                             pdb_gz))

        if entry.integrated_in:
            # Integrated signature: add prediction for InterPro entry
            my_cur.execute(req, (entry.integrated_in, algorithm, msa_gz,
                                 cmap_gz, plddt_gz, pdb_gz))

    ora_cur.close()
    ora_con.close()

    my_con.commit()

    logger.info("indexing")
    my_cur.execute(
        """
        CREATE INDEX i_structuralmodel_entry
        ON webfront_structuralmodel (accession)
        """
    )
    my_cur.close()
    my_con.close()

    logger.info("complete")
