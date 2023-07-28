import pickle
import shelve

import MySQLdb

from interpro7dw.utils import logger
from interpro7dw.utils.store import BasicStore
from interpro7dw.utils.mysql import uri2dict
from .utils import jsonify


def populate_rosettafold(uri: str, models_file: str):
    logger.info("creating webfront_structuralmodel")
    con = MySQLdb.connect(**uri2dict(uri))
    cur = con.cursor()
    cur.execute("DROP TABLE IF EXISTS webfront_structuralmodel")
    cur.execute(
        """
        CREATE TABLE webfront_structuralmodel
        (
            model_id INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
            accession VARCHAR(30) NOT NULL,
            algorithm VARCHAR(20) NOT NULL,
            contacts LONGBLOB NOT NULL,
            plddt LONGBLOB NOT NULL,
            structure LONGBLOB NOT NULL
        ) CHARSET=utf8mb4 DEFAULT COLLATE=utf8mb4_unicode_ci
        """
    )

    query = """
        INSERT INTO webfront_structuralmodel (
          accession, algorithm, contacts, plddt, structure
        ) VALUES (%s, 'RoseTTAFold', %s, %s, %s)
    """

    with BasicStore(models_file, mode="r") as models:
        for s_acc, e_acc, cmap_gz, plddt_gz, pdb_gz in models:
            cur.execute(query, (s_acc, cmap_gz, plddt_gz, pdb_gz))

            if e_acc:
                # Integrated signature: add prediction for InterPro entry
                cur.execute(query, (e_acc, cmap_gz, plddt_gz, pdb_gz))

    con.commit()

    logger.info("creating index")
    cur.execute(
        """
        CREATE INDEX i_structuralmodel_entry
        ON webfront_structuralmodel (accession)
        """
    )
    cur.close()
    con.close()

    logger.info("done")


def populate_structures(uri: str, structures_file: str,
                        uniprot2pdb_file: str, pdbmatches_file: str,
                        xrefs_file: str):
    logger.info("loading UniProt-PDB mapping")
    pdb2chains = {}
    pdb2segments = {}
    with open(uniprot2pdb_file, "rb") as fh:
        for protein_acc, structures in pickle.load(fh).items():
            for pdb_chain, segments in structures.items():
                pdb_id, chain = pdb_chain.split("_")

                try:
                    proteins = pdb2segments[pdb_id]
                except KeyError:
                    proteins = pdb2segments[pdb_id] = {}
                    pdb2chains[pdb_id] = set()

                pdb2chains[pdb_id].add(chain)

                try:
                    chains = proteins[protein_acc]
                except KeyError:
                    chains = proteins[protein_acc] = {}

                try:
                    chains[chain] += segments
                except KeyError:
                    chains[chain] = segments

    con = MySQLdb.connect(**uri2dict(uri), charset="utf8mb4")
    cur = con.cursor()

    logger.info("creating webfront_chain_sequence")
    cur.execute("DROP TABLE IF EXISTS webfront_chain_sequence")
    cur.execute(
        """
        CREATE TABLE webfront_chain_sequence
        (
            id INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
            structure_acc VARCHAR(4) NOT NULL,
            chain_acc VARCHAR(15) COLLATE utf8mb4_bin NOT NULL,
            sequence LONGBLOB NOT NULL,
            length INT(11) NOT NULL
        ) CHARSET=utf8mb4 DEFAULT COLLATE=utf8mb4_unicode_ci
        """
    )

    with shelve.open(pdbmatches_file, writeback=False) as d:
        query = """
            INSERT INTO webfront_chain_sequence (
                structure_acc, chain_acc, sequence, length
            )
            VALUES (%s, %s, %s, %s)
        """
        for pdb_chain, pdb_entry in d.items():
            pdb_id, chain = pdb_chain.split("_")

            try:
                pdb2chains[pdb_id].add(chain)
            except KeyError:
                pdb2chains[pdb_id] = {chain}

            cur.execute(query, (pdb_id, chain, pdb_entry["sequence"],
                                pdb_entry["length"]))

    con.commit()

    logger.info("indexing webfront_chain_sequence")
    cur.execute(
        """
        CREATE UNIQUE INDEX ui_chain_sequence
        ON webfront_chain_sequence (structure_acc, chain_acc)        
        """
    )

    logger.info("loading structures")
    with open(structures_file, "rb") as fh:
        structures = pickle.load(fh)

    logger.info("creating webfront_structure")
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

    query = """
        INSERT INTO webfront_structure 
        VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s) 
    """
    params = []

    with BasicStore(xrefs_file, mode="r") as store:
        for pdb_id, xrefs in store:
            # Adds total number of entries
            num_entries = {"total": 0}
            for database, entries in xrefs["entries"].items():
                num_entries[database.lower()] = len(entries)
                num_entries["total"] += len(entries)

            structure = structures[pdb_id]

            proteins = pdb2segments.get(pdb_id, {})
            for protein_acc in proteins:
                for chain, segments in proteins[protein_acc].items():
                    segments.sort(key=lambda x: (x["protein_start"],
                                                 x["protein_end"]))

            params.append((
                pdb_id,
                structure["name"],
                "pdb",
                structure["evidence"],
                structure["date"],
                structure["resolution"],
                jsonify(structure["citations"], nullable=True),
                jsonify(sorted(pdb2chains.get(pdb_id, [])), nullable=False),
                jsonify(proteins, nullable=False),
                jsonify(structure["secondary_structures"], nullable=True),
                jsonify({
                    "domain_architectures": len(xrefs["dom_orgs"]),
                    "entries": num_entries,
                    "proteomes": len(xrefs["proteomes"]),
                    "proteins": xrefs["proteins"],
                    "sets": len(xrefs["sets"]),
                    "taxa": len(xrefs["taxa"])
                })
            ))

            if len(params) == 1000:
                cur.executemany(query, params)
                params = []

    if params:
        cur.executemany(query, params)

    con.commit()
    cur.close()
    con.close()

    logger.info("done")
