import pickle
import shelve

from interpro7dw.utils import logger
from interpro7dw.utils.store import BasicStore
from .utils import connect, jsonify


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

    con = connect(uri)
    cur = con.cursor()

    logger.info("creating webfront_chain_sequence")
    cur.execute("DROP TABLE IF EXISTS webfront_chain_sequence")
    cur.execute(
        """
        CREATE TABLE webfront_chain_sequence
        (
            id SERIAL NOT NULL PRIMARY KEY,
            structure_acc VARCHAR(4) NOT NULL,
            chain_acc VARCHAR(15) COLLATE "C" NOT NULL,
            sequence TEXT NOT NULL,
            length INTEGER NOT NULL
        )
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
            release_date TIMESTAMP NOT NULL,
            resolution FLOAT,
            literature JSONB,
            chains JSONB NOT NULL,
            proteins JSONB NOT NULL,
            secondary_structures JSONB,
            counts JSONB NOT NULL
        )
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
                    "proteins": len(xrefs["proteins"]),
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
