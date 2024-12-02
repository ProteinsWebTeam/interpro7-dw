import pickle

import MySQLdb

from interpro7dw.utils import logger
from interpro7dw.utils.mysql import uri2dict
from interpro7dw.utils.store import BasicStore, KVStore

from .utils import create_index, jsonify


def populate_features(uri: str, features_file: str):
    logger.info("creating webfront_proteinfeature")

    con = MySQLdb.connect(**uri2dict(uri), charset="utf8mb4")
    cur = con.cursor()
    cur.execute("DROP TABLE IF EXISTS webfront_proteinfeature")
    cur.execute(
        """
        CREATE TABLE webfront_proteinfeature
        (
            feature_id INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
            protein_acc VARCHAR(15) NOT NULL,
            entry_acc VARCHAR(30) NOT NULL,
            source_database VARCHAR(10) NOT NULL,
            location_start INT NOT NULL,
            location_end INT NOT NULL,
            sequence_feature VARCHAR(35)
        ) CHARSET=utf8mb4 DEFAULT COLLATE=utf8mb4_unicode_ci
        """
    )

    query = """
        INSERT INTO webfront_proteinfeature (
          protein_acc, entry_acc, source_database, location_start,
          location_end, sequence_feature
        )
        VALUES (%s, %s, %s, %s, %s, %s)
    """
    params = []

    i = 0
    with KVStore(features_file) as store:
        for i, (protein_acc, features) in enumerate(store.items()):
            for feature in features:
                database = feature["database"].lower()

                if database == "antifam":
                    # AntiFam matches are in Elastic (like member databases)
                    continue

                for pos_start, pos_end, seq_feature in feature["locations"]:
                    if database == "elm":
                        seq_feature = feature["name"]

                    params.append((
                        protein_acc,
                        feature["accession"],
                        database,
                        pos_start,
                        pos_end,
                        seq_feature
                    ))

                    if len(params) == 1000:
                        cur.executemany(query, params)
                        params = []

            if (i + 1) % 1e7 == 0:
                logger.info(f"{i + 1:>15,}")

    if params:
        cur.executemany(query, params)

    logger.info(f"{i + 1:>15,}")
    con.commit()
    cur.close()
    con.close()

    logger.info("done")


def index_features(uri: str):
    con = MySQLdb.connect(**uri2dict(uri), charset="utf8mb4")
    cur = con.cursor()
    create_index(
        cur,
        """
        CREATE INDEX i_proteinfeature
        ON webfront_proteinfeature (protein_acc)
        """
    )
    cur.close()
    con.close()


def populate_toad_matches(uri: str, toad_file: str):
    logger.info("creating webfront_interpro_n")

    con = MySQLdb.connect(**uri2dict(uri), charset="utf8mb4")
    cur = con.cursor()
    cur.execute("DROP TABLE IF EXISTS webfront_interpro_n")
    cur.execute(
        """
        CREATE TABLE webfront_interpro_n
        (
            match_id INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
            protein_acc VARCHAR(15) NOT NULL,
            entry_acc VARCHAR(30) NOT NULL,
            locations LONGTEXT NOT NULL
        ) CHARSET=utf8mb4 DEFAULT COLLATE=utf8mb4_unicode_ci
        """
    )

    query = """
            INSERT INTO webfront_interpro_n (protein_acc, entry_acc, locations)
            VALUES (%s, %s, %s)
        """
    params = []

    i = 0
    with KVStore(toad_file) as store:
        for i, (protein_acc, (signatures, entries)) in enumerate(store.items()):
            for signature_acc, signature in signatures.items():
                params.append((
                    protein_acc,
                    signature_acc,
                    jsonify(signature["locations"], nullable=False)
                ))

                if len(params) == 1000:
                    cur.executemany(query, params)
                    params = []

            if (i + 1) % 1e7 == 0:
                logger.info(f"{i + 1:>15,}")

    if params:
        cur.executemany(query, params)

    logger.info(f"{i + 1:>15,}")
    con.commit()
    cur.close()
    con.close()

    logger.info("done")


def index_toad_matches(uri: str):
    con = MySQLdb.connect(**uri2dict(uri), charset="utf8mb4")
    cur = con.cursor()
    create_index(
        cur,
        """
        CREATE INDEX i_interpro_n
        ON webfront_interpro_n (protein_acc)
        """
    )
    cur.close()
    con.close()


def populate_isoforms(uri: str, isoforms_file: str):
    logger.info("creating webfront_varsplic")

    con = MySQLdb.connect(**uri2dict(uri))
    cur = con.cursor()
    cur.execute("DROP TABLE IF EXISTS webfront_varsplic")
    cur.execute(
        """
        CREATE TABLE webfront_varsplic
        (
            accession VARCHAR(20) PRIMARY KEY NOT NULL,
            protein_acc VARCHAR(15) NOT NULL,
            length INT(11) NOT NULL,
            sequence LONGTEXT NOT NULL,
            features LONGTEXT
        ) CHARSET=utf8mb4 DEFAULT COLLATE=utf8mb4_unicode_ci
        """
    )

    query = "INSERT INTO webfront_varsplic VALUES (%s, %s, %s, %s, %s)"
    params = []

    with BasicStore(isoforms_file, mode="r") as store:
        for isoform in store:
            signatures, entries = isoform["matches"]

            features = {}
            for obj in [signatures, entries]:
                for entry_acc, entry in obj.items():
                    features[entry_acc] = {
                        "accession": entry_acc,
                        "integrated": entry.get("entry"),
                        "name": entry["name"],
                        "type": entry["type"].lower(),
                        "source_database": entry["database"].lower(),
                        "locations": entry["locations"]
                    }

            params.append((
                isoform["accession"],
                isoform["protein"],
                isoform["length"],
                isoform["sequence"],
                jsonify(features)
            ))

            if len(params) == 1000:
                cur.executemany(query, params)
                params = []

    if params:
        cur.executemany(query, params)

    con.commit()

    logger.info("creating index")
    cur.execute(
        """
        CREATE INDEX i_varsplic
        ON webfront_varsplic (protein_acc)
        """
    )
    cur.close()
    con.close()

    logger.info("done")


def populate_proteins(uri: str, clans_file: str, entries_file: str,
                      isoforms_file: str, cath_scop_file: str,
                      uniprot2pdb_file: str, taxa_file: str, proteins_file: str,
                      domorgs_file: str, evidences_file: str,
                      functions_file: str, matches_file: str, names_file: str,
                      proteomes_file: str, sequences_file: str,
                      alphafold_file: str):
    """Creates and populates the MySQL webfront_protein table.

    :param uri: InterPro MySQL connection string.
    :param clans_file: File of clan information.
    :param entries_file: File of entries information.
    :param isoforms_file: BasicStore file of protein isoforms.
    :param cath_scop_file: File of CATH and SCOP domains.
    :param uniprot2pdb_file: File of protein-structures mapping.
    :param taxa_file: File of taxonomic information.
    :param proteins_file: KVStore file of proteins.
    :param domorgs_file: KVStore file of domain organisations.
    :param evidences_file: KVStore file of protein evidences/genes.
    :param functions_file: KVStore file of protein functions.
    :param matches_file: KVStore file of protein matches.
    :param names_file: KVStore file of protein descriptions/names.
    :param proteomes_file: KVStore file of protein-proteome mapping.
    :param sequences_file: KVStore file of protein sequences.
    :param alphafold_file: KVStore file of AlphaFold pLDDT scores.
    """
    logger.info("loading clans and entries")
    member2clan = {}
    with open(clans_file, "rb") as fh:
        for clan_acc, clan in pickle.load(fh).items():
            for entry_acc, _, _, _, _ in clan["members"]:
                member2clan[entry_acc] = clan_acc

    entry2go = {}
    with open(entries_file, "rb") as fh:
        for entry in pickle.load(fh).values():
            if entry.go_terms:
                entry2go[entry.accession] = entry.go_terms

    logger.info("loading CATH/SCOP domains")
    with open(cath_scop_file, "rb") as fh:
        uniprot2cath, uniprot2scop = pickle.load(fh)

    logger.info("loading UniProt-PDB mapping")
    num_structures = {}
    with open(uniprot2pdb_file, "rb") as fh:
        for protein_acc, structures in pickle.load(fh).items():
            protein_structures = set()

            for pdb_chain in structures:
                pdb_id, chain = pdb_chain.split("_")
                protein_structures.add(pdb_id)

            num_structures[protein_acc] = len(protein_structures)

    logger.info("loading isoforms")
    num_isoforms = {}
    with BasicStore(isoforms_file, mode="r") as store:
        for isoform in store:
            protein_acc = isoform["protein"]
            try:
                num_isoforms[protein_acc] += 1
            except KeyError:
                num_isoforms[protein_acc] = 1

    logger.info("loading taxa")
    taxa = {}
    with open(taxa_file, "rb") as fh:
        for taxon_id, taxon in pickle.load(fh).items():
            taxa[taxon_id] = jsonify({
                "taxId": taxon_id,
                "scientificName": taxon["sci_name"],
                "fullName": taxon["full_name"]
            })

    logger.info("creating webfront_protein")
    proteins_store = KVStore(proteins_file)
    functions_store = KVStore(functions_file)
    names_store = KVStore(names_file)
    evidences_store = KVStore(evidences_file)
    domorgs_store = KVStore(domorgs_file)
    matches_store = KVStore(matches_file)
    proteomes_store = KVStore(proteomes_file)
    sequences_store = KVStore(sequences_file)
    alphafold_store = KVStore(alphafold_file)

    con = MySQLdb.connect(**uri2dict(uri), charset="utf8")
    cur = con.cursor()
    cur.execute("DROP TABLE IF EXISTS webfront_protein")
    cur.execute(
        """
        CREATE TABLE webfront_protein
        (
            accession VARCHAR(15) PRIMARY KEY NOT NULL,
            identifier VARCHAR(16) NOT NULL,
            organism LONGTEXT NOT NULL,
            name VARCHAR(255) NOT NULL,
            description LONGTEXT,
            sequence LONGBLOB NOT NULL,
            length INT(11) NOT NULL,
            proteome VARCHAR(20),
            gene VARCHAR(70),
            go_terms LONGTEXT,
            evidence_code INT(11) NOT NULL,
            source_database VARCHAR(10) NOT NULL,
            is_fragment TINYINT NOT NULL,
            structure LONGTEXT,
            tax_id VARCHAR(20) NOT NULL,
            ida_id VARCHAR(40),
            ida TEXT,
            in_alphafold TINYINT NOT NULL,
            counts LONGTEXT NOT NULL
        ) CHARSET=utf8 DEFAULT COLLATE=utf8_unicode_ci
        """
    )

    query = """
        INSERT into webfront_protein
        VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)
    """
    params = []
    i = 0
    for i, (protein_acc, protein) in enumerate(proteins_store.items()):
        taxon_id = protein["taxid"]

        try:
            taxon = taxa[taxon_id]
        except KeyError:
            cur.close()
            con.close()
            raise KeyError(f"{protein_acc}: invalid taxon {taxon_id}")

        try:
            name = names_store[protein_acc]
        except KeyError:
            cur.close()
            con.close()
            raise KeyError(f"{protein_acc}: missing name")

        try:
            evidence, gene = evidences_store[protein_acc]
        except KeyError:
            cur.close()
            con.close()
            raise RuntimeError(f"{protein_acc}: missing evidence/gene")

        try:
            gzipped_sequence = sequences_store[protein_acc]
        except KeyError:
            cur.close()
            con.close()
            raise RuntimeError(f"{protein_acc}: missing sequence")

        proteome_id = proteomes_store.get(protein_acc)
        sig_matches, entry_matches = matches_store.get(protein_acc, ({}, {}))

        clans = []
        databases = {}
        go_terms = {}
        for obj in [sig_matches, entry_matches]:
            for entry_acc, match in obj.items():
                if entry_acc in member2clan:
                    clans.append(member2clan[entry_acc])

                database = match["database"].lower()
                if database in databases:
                    databases[database] += 1
                else:
                    databases[database] = 1

                # if database == "panther":
                #     # Check for a PANTHER subfamily
                #     for loc in match["locations"]:
                #         if loc["model"]:
                #             entry_acc = loc["model"]
                #             break

                for term in entry2go.get(entry_acc, []):
                    term_id = term["identifier"]

                    if term_id not in go_terms:
                        go_terms[term_id] = term.copy()
                    #     go_terms[term_id]["sources"] = set()
                    #
                    # go_terms[term_id]["sources"].add(database)

        # # Convert sets to lists (JSON does not support sets)
        # for term in go_terms.values():
        #     term["sources"] = list(term["sources"])

        # Adds CATH/SCOP structures
        cath_scop_domains = {}
        for key, obj in [("cath", uniprot2cath), ("scop", uniprot2scop)]:
            domains = obj.get(protein_acc)
            if domains:
                cath_scop_domains[key] = {}

                for dom in domains.values():
                    dom_id = dom["id"]

                    cath_scop_domains[key][dom_id] = {
                        "domain_id": dom["superfamily"]["id"],
                        "coordinates": dom["locations"]
                    }

        try:
            domain = domorgs_store[protein_acc]
        except KeyError:
            dom_id = dom_key = None
            dom_proteins = 0
        else:
            dom_id = domain["id"]
            dom_key = domain["key"]
            dom_proteins = domain["count"]

        params.append((
            protein_acc,
            protein["identifier"],
            taxon,
            name,
            jsonify(functions_store.get(protein_acc), nullable=True),
            gzipped_sequence,
            protein["length"],
            proteome_id,
            gene,
            jsonify(list(go_terms.values()), nullable=True),
            evidence,
            "reviewed" if protein["reviewed"] else "unreviewed",
            1 if protein["fragment"] else 0,
            jsonify(cath_scop_domains, nullable=True),
            taxon_id,
            dom_id,
            dom_key,
            len(alphafold_store.get(protein_acc, [])) != 0,
            jsonify({
                "domain_architectures": dom_proteins,
                "entries": databases,
                "isoforms": num_isoforms.get(protein_acc, 0),
                "proteomes": 1 if proteome_id else 0,
                "sets": len(set(clans)),
                "structures": num_structures.get(protein_acc, 0),
                "taxa": 1
            })
        ))

        if len(params) == 1000:
            cur.executemany(query, params)
            params = []

        if (i + 1) % 1e7 == 0:
            logger.info(f"{i + 1:>15,} / {len(proteins_store):,}")

    if params:
        cur.executemany(query, params)

    logger.info(f"{i + 1:>15,} / {len(proteins_store):,}")
    con.commit()
    cur.close()
    con.close()

    proteins_store.close()
    functions_store.close()
    names_store.close()
    evidences_store.close()
    domorgs_store.close()
    matches_store.close()
    proteomes_store.close()
    sequences_store.close()

    logger.info("done")


def index_proteins(uri: str):
    con = MySQLdb.connect(**uri2dict(uri), charset="utf8mb4")
    cur = con.cursor()
    logger.info("ui_protein_identifier")
    create_index(
        cur,
        """
        CREATE UNIQUE INDEX ui_protein_identifier
        ON webfront_protein (identifier)
        """
    )
    logger.info("u_protein_name")
    create_index(
        cur,
        """
        CREATE INDEX u_protein_name
        ON webfront_protein (name)
        """
    )
    logger.info("i_protein_proteome")
    create_index(
        cur,
        """
        CREATE INDEX i_protein_proteome
        ON webfront_protein (proteome)
        """
    )
    logger.info("i_protein_database")
    create_index(
        cur,
        """
        CREATE INDEX i_protein_database
        ON webfront_protein (source_database)
        """
    )
    logger.info("i_protein_taxon")
    create_index(
        cur,
        """
        CREATE INDEX i_protein_taxon
        ON webfront_protein (tax_id)
        """
    )
    logger.info("i_protein_taxon_gene")
    create_index(
        cur,
        """
        CREATE INDEX i_protein_taxon_gene
        ON webfront_protein (tax_id, gene)
        """
    )
    logger.info("i_protein_ida")
    create_index(
        cur,
        """
        CREATE INDEX i_protein_ida
        ON webfront_protein (ida_id)
        """
    )
    logger.info("i_protein_fragment")
    create_index(
        cur,
        """
        CREATE INDEX i_protein_fragment
        ON webfront_protein (is_fragment)
        """
    )
    cur.close()
    con.close()
    logger.info("done")


def populate_residues(uri: str, residues_file: str):
    logger.info("creating webfront_proteinresidue")

    con = MySQLdb.connect(**uri2dict(uri), charset="utf8mb4")
    cur = con.cursor()
    cur.execute("DROP TABLE IF EXISTS webfront_proteinresidue")
    cur.execute(
        """
        CREATE TABLE webfront_proteinresidue
        (
            residue_id INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
            protein_acc VARCHAR(15) NOT NULL,
            entry_acc VARCHAR(30) NOT NULL,
            entry_name VARCHAR(100),
            source_database VARCHAR(10) NOT NULL,
            description VARCHAR(255),
            fragments LONGTEXT NOT NULL
        ) CHARSET=utf8mb4 DEFAULT COLLATE=utf8mb4_unicode_ci
        """
    )

    query = """
        INSERT INTO webfront_proteinresidue (
          protein_acc, entry_acc, entry_name, source_database, description,
          fragments
        )
        VALUES (%s, %s, %s, %s, %s, %s)
    """
    params = []

    i = 0
    with BasicStore(residues_file, mode="r") as store:
        for i, (protein_acc, entries) in enumerate(store):
            for entry_acc, entry in entries.items():
                for descr, locations in entry["descriptions"].items():
                    params.append((
                        protein_acc,
                        entry_acc,
                        entry["name"] or entry_acc,
                        entry["database"].lower(),
                        descr,
                        jsonify(locations, nullable=False)
                    ))

                    if len(params) == 1000:
                        cur.executemany(query, params)
                        params = []

            if (i + 1) % 1e7 == 0:
                logger.info(f"{i + 1:>15,}")

    if params:
        cur.executemany(query, params)

    logger.info(f"{i + 1:>15,}")
    con.commit()
    cur.close()
    con.close()

    logger.info("done")


def index_residues(uri: str):
    con = MySQLdb.connect(**uri2dict(uri), charset="utf8mb4")
    cur = con.cursor()
    create_index(
        cur,
        """
        CREATE INDEX i_proteinresidue
        ON webfront_proteinresidue (protein_acc)
        """
    )
    cur.close()
    con.close()
