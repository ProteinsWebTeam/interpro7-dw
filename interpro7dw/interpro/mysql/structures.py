import MySQLdb

from interpro7dw.utils import logger
from interpro7dw.utils.store import loadobj, SimpleStore
from interpro7dw.utils.mysql import url2dict
from .utils import jsonify


def insert_structural_models(url: str, entries_file: str, models_file: str):
    logger.info("loading entries")
    entries = loadobj(entries_file)

    logger.info("creating webfront_structuralmodel")
    con = MySQLdb.connect(**url2dict(url))
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
        ) VALUES (%s, %s, %s, %s, %s)
    """

    with SimpleStore(models_file) as models:
        for entry_acc, algorithm, cmap_gz, errs_gz, plddt_gz, pdb_gz in models:
            try:
                entry = entries[entry_acc]
            except KeyError:
                continue

            cur.execute(query, (entry_acc, algorithm, cmap_gz, plddt_gz,
                                pdb_gz))

            if entry.integrated_in:
                # Integrated signature: add prediction for InterPro entry
                cur.execute(query, (entry.integrated_in, algorithm, cmap_gz,
                                    plddt_gz, pdb_gz))

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


def insert_structures(url: str, structures_file: str, xrefs_file: str):
    logger.info("creating webfront_structure")
    con = MySQLdb.connect(**url2dict(url), charset="utf8mb4")
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
    args = []

    structures = loadobj(structures_file)
    with SimpleStore(xrefs_file) as store:
        for pdbe_id, xrefs in store:
            structure = structures[pdbe_id]

            # Adds total number of entries
            entries = {}
            total = 0
            for db, accessions in xrefs["entries"].items():
                entries[db] = len(accessions)
                total += entries[db]

            xrefs["entries"]["total"] = sum(xrefs["entries"].values())

            args.append((
                pdbe_id,
                structure["name"],
                "pdb",
                structure["evidence"],
                structure["date"],
                structure["resolution"],
                jsonify(structure["citations"], nullable=True),
                # Sorted list of unique chain (e.g. 'A', 'B', ...)
                jsonify(sorted({chain_id
                                for chains in structure["proteins"].values()
                                for chain_id in chains}), nullable=False),
                jsonify(structure["proteins"], nullable=False),
                jsonify(structure["secondary_structures"]),
                jsonify({
                    "domain_architectures": len(xrefs["domain_architectures"]),
                    "entries": entries,
                    "proteomes": len(xrefs["proteomes"]),
                    "proteins": xrefs["proteins"],
                    "sets": len(xrefs["sets"]),
                    "taxa": len(xrefs["taxa"])
                })
            ))

            if len(args) == 1000:
                cur.executemany(query, args)
                args.clear()

    if args:
        cur.executemany(query, args)
        args.clear()

    con.commit()
    cur.close()
    con.close()

    logger.info("done")
