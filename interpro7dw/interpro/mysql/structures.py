import pickle

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
            accession VARCHAR(25) NOT NULL,
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
                        protein2structures_file: str, xrefs_file: str):
    logger.info("loading structures")
    with open(structures_file, "rb") as fh:
        data = pickle.load(fh)

    structure2segments = {}
    with open(protein2structures_file, "rb") as fh:
        for protein_acc, protein_structures in pickle.load(fh).items():
            for pdb_id, chains in protein_structures.items():
                if pdb_id in structure2segments:
                    proteins = structure2segments[pdb_id]
                else:
                    proteins = structure2segments[pdb_id] = {}

                if protein_acc in proteins:
                    struct_chains = proteins[protein_acc]
                else:
                    struct_chains = proteins[protein_acc] = {}

                for chain_id, segments in chains.items():
                    if chain_id in struct_chains:
                        struct_chains[chain_id] += segments
                    else:
                        struct_chains[chain_id] = segments

    logger.info("creating webfront_structure")
    con = MySQLdb.connect(**uri2dict(uri), charset="utf8mb4")
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

    query = """
        INSERT INTO webfront_structure 
        VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s) 
    """
    params = []

    with BasicStore(xrefs_file, mode="r") as store:
        pdb_entries = data["entries"]

        for pdb_id, xrefs in store:
            # Adds total number of entries
            num_entries = {"total": 0}
            for database, entries in xrefs["entries"].items():
                num_entries[database.lower()] = len(entries)
                num_entries["total"] += len(entries)

            structure = pdb_entries[pdb_id]

            proteins = structure2segments.get(pdb_id, {})
            chains = set()
            for protein_acc in proteins:
                for chain_id, segments in proteins[protein_acc].items():
                    chains.add(chain_id)
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
                jsonify(sorted(chains), nullable=False),
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
