import pickle

from interpro7dw.utils import logger
from interpro7dw.utils.store import BasicStore, KVStore

from .utils import connect, jsonify


def populate_features(uri: str, features_file: str):
    logger.info("creating webfront_proteinfeature")

    con = connect(uri)
    cur = con.cursor()
    cur.execute("DROP TABLE IF EXISTS interpro.webfront_proteinfeature")
    cur.execute(
        """
        CREATE TABLE interpro.webfront_proteinfeature
        (
            feature_id SERIAL NOT NULL PRIMARY KEY,
            protein_acc VARCHAR(15) COLLATE "case_insensitive" NOT NULL,
            entry_acc VARCHAR(30) NOT NULL,
            source_database VARCHAR(10) NOT NULL,
            location_start INTEGER NOT NULL,
            location_end INTEGER NOT NULL,
            sequence_feature VARCHAR(255)
        )
        """
    )

    query = """
        INSERT INTO interpro.webfront_proteinfeature (
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

                if database == "antifam" or database == "pfam-n":
                    # AntiFam matches are in Elastic (like member databases)
                    # Pfam-N to be removed from Oracle DB
                    continue

                for pos_start, pos_end, seq_feature in feature["locations"]:
                    if database == "elm":
                        seq_feature = feature["name"]
                    elif database == "cathfunfam":
                        database = "funfam"  # TODO: remove when website has been updated
                        seq_feature = feature["description"]

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
    con = connect(uri)
    cur = con.cursor()
    cur.execute(
        """
        CREATE INDEX IF NOT EXISTS i_proteinfeature
        ON interpro.webfront_proteinfeature (UPPER(protein_acc))
        """
    )
    con.commit()
    cur.close()
    con.close()


def populate_toad_matches(uri: str, matches_file: str, toad_file: str):
    logger.info("creating webfront_interpro_n")

    con = connect(uri)
    cur = con.cursor()
    cur.execute("DROP TABLE IF EXISTS interpro.webfront_interpro_n")
    cur.execute(
        """
        CREATE TABLE interpro.webfront_interpro_n
        (
            match_id SERIAL NOT NULL PRIMARY KEY,
            protein_acc VARCHAR(15) COLLATE "case_insensitive" NOT NULL,
            entry_acc VARCHAR(30) NOT NULL,
            locations JSONB NOT NULL,
            in_interpro BOOLEAN NOT NULL,
            is_preferred BOOLEAN NOT NULL
        )
        """
    )

    query = """
        INSERT INTO interpro.webfront_interpro_n (
            protein_acc, entry_acc, locations, in_interpro, is_preferred
        ) VALUES (%s, %s, %s, %s, %s)
    """
    params = []

    i = 0
    with KVStore(toad_file) as ts, KVStore(matches_file) as ms:
        for i, (protein_acc, toad_matches) in enumerate(ts.items()):
            # Keep signature matches from traditional InterPro
            trad_matches = {}
            for match in ms.get(protein_acc, []):
                if match["database"].lower() != "interpro":
                    trad_matches[match["accession"]] = match

            # Compare InterPro-N matches with InterPro ones
            for toad_match in toad_matches:
                if toad_match["database"].lower() != "interpro":
                    # We don't use InterPro entries matches for InterPro-N
                    match_acc = toad_match["accession"]

                    if match_acc in trad_matches:
                        trad_match = trad_matches[match_acc]
                        trad_cov = calc_coverage(trad_match["locations"])
                        toad_cov = calc_coverage(toad_match["locations"])

                        in_interpro = True
                        is_preferred = toad_cov > (trad_cov * 1.05)
                    else:
                        in_interpro = False
                        is_preferred = True

                    params.append((
                        protein_acc,
                        match_acc,
                        jsonify(toad_match["locations"], nullable=False),
                        in_interpro,
                        is_preferred
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


def calc_coverage(locations: list[dict]) -> int:
    cov = 0
    for loc in locations:
        for frag in loc["fragments"]:
            cov += frag["end"] - frag["start"] + 1

    return cov


def index_toad_matches(uri: str):
    con = connect(uri)
    cur = con.cursor()
    cur.execute(
        """
        CREATE INDEX IF NOT EXISTS i_interpro_n
        ON interpro.webfront_interpro_n (UPPER(protein_acc))
        """
    )
    con.commit()
    cur.close()
    con.close()


def populate_isoforms(uri: str, isoforms_file: str):
    logger.info("creating webfront_varsplic")

    con = connect(uri)
    cur = con.cursor()
    cur.execute("DROP TABLE IF EXISTS interpro.webfront_varsplic")
    cur.execute(
        """
        CREATE TABLE interpro.webfront_varsplic
        (
            accession VARCHAR(20)  PRIMARY KEY NOT NULL,
            protein_acc VARCHAR(15) COLLATE "case_insensitive" NOT NULL,
            length INTEGER NOT NULL,
            sequence TEXT NOT NULL,
            features JSONB
        )
        """
    )

    query = "INSERT INTO interpro.webfront_varsplic VALUES (%s, %s, %s, %s, %s)"
    params = []

    with BasicStore(isoforms_file, mode="r") as store:
        for isoform in store:
            features = {}
            for match in isoform["matches"]:
                match_acc = match["accession"]
                features[match_acc] = {
                    "accession": match_acc,
                    "integrated": match.get("entry"),
                    "name": match["name"],
                    "type": match["type"].lower(),
                    "source_database": match["database"].lower(),
                    "locations": match["locations"]
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

    logger.info("creating index")
    cur.execute(
        """
        CREATE INDEX i_varsplic
        ON interpro.webfront_varsplic (UPPER(protein_acc))
        """
    )
    con.commit()
    cur.close()
    con.close()

    logger.info("done")


def populate_proteins(uri: str, clans_file: str, entries_file: str,
                      isoforms_file: str, cath_scop_file: str,
                      uniprot2pdb_file: str, taxa_file: str, proteins_file: str,
                      domorgs_file: str, evidences_file: str,
                      functions_file: str, matches_file: str, names_file: str,
                      proteomes_file: str, sequences_file: str,
                      alphafold_file: str, bfvd_file: str):
    """Creates and populates the webfront_protein table.

    :param uri: InterPro PostgreSQL connection string.
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
    :param bfvd_file: KVStore file of BFVD predictions.
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
    bfvd_store = KVStore(bfvd_file)

    con = connect(uri)
    cur = con.cursor()
    cur.execute("DROP TABLE IF EXISTS interpro.webfront_protein")
    cur.execute(
        """
        CREATE TABLE interpro.webfront_protein
        (
            accession VARCHAR(15) COLLATE "case_insensitive" NOT NULL,
            identifier VARCHAR(16) COLLATE "case_insensitive" NOT NULL,
            organism TEXT NOT NULL,
            name VARCHAR(255) COLLATE "case_insensitive" NOT NULL,
            description JSONB,
            sequence TEXT NOT NULL,
            length INTEGER NOT NULL,
            proteome VARCHAR(20) COLLATE "case_insensitive",
            gene VARCHAR(70) COLLATE "case_insensitive",
            go_terms JSONB,
            evidence_code INTEGER NOT NULL,
            source_database VARCHAR(10) COLLATE "case_insensitive" NOT NULL,
            is_fragment BOOLEAN NOT NULL,
            structure JSONB,
            tax_id VARCHAR(20) COLLATE "case_insensitive" NOT NULL,
            ida_id VARCHAR(40) COLLATE "case_insensitive",
            ida TEXT,
            in_alphafold BOOLEAN NOT NULL,
            in_bfvd BOOLEAN NOT NULL,
            counts JSONB NOT NULL
        )
        """
    )

    query = """
        INSERT into interpro.webfront_protein
        VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)
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

        clans = []
        databases = {}
        go_terms = {}
        for match in matches_store.get(protein_acc, []):
            match_acc = match["accession"]

            if match_acc in member2clan:
                clans.append(member2clan[match_acc])

            database = match["database"].lower()
            if database in databases:
                databases[database] += 1
            else:
                databases[database] = 1

            for term in entry2go.get(match_acc, []):
                term_id = term["identifier"]

                if term_id not in go_terms:
                    go_terms[term_id] = term.copy()

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
            protein["fragment"],
            jsonify(cath_scop_domains, nullable=True),
            taxon_id,
            dom_id,
            dom_key,
            len(alphafold_store.get(protein_acc, [])) != 0,
            bfvd_store.get(protein_acc) is not None,
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
    con = connect(uri)
    cur = con.cursor()
    logger.info("ui_protein_accession")
    cur.execute(
        """
        CREATE UNIQUE INDEX ui_protein_accession
        ON interpro.webfront_protein (UPPER(accession))
        """
    )
    logger.info("ui_protein_identifier")
    cur.execute(
        """
        CREATE UNIQUE INDEX ui_protein_identifier
        ON interpro.webfront_protein (UPPER(identifier))
        """
    )
    logger.info("u_protein_name")
    cur.execute(
        """
        CREATE INDEX u_protein_name
        ON interpro.webfront_protein (UPPER(name))
        """
    )
    logger.info("i_protein_proteome")
    cur.execute(
        """
        CREATE INDEX i_protein_proteome
        ON interpro.webfront_protein (UPPER(proteome))
        """
    )
    logger.info("i_protein_database")
    cur.execute(
        """
        CREATE INDEX i_protein_database
        ON interpro.webfront_protein (UPPER(source_database))
        """
    )
    logger.info("i_protein_taxon")
    cur.execute(
        """
        CREATE INDEX i_protein_taxon
        ON interpro.webfront_protein (UPPER(tax_id))
        """
    )
    logger.info("i_protein_taxon_gene")
    cur.execute(
        """
        CREATE INDEX i_protein_taxon_gene
        ON interpro.webfront_protein (UPPER(tax_id), UPPER(gene))
        """
    )
    logger.info("i_protein_ida")
    cur.execute(
        """
        CREATE INDEX i_protein_ida
        ON interpro.webfront_protein (UPPER(ida_id))
        """
    )
    logger.info("i_protein_fragment")
    cur.execute(
        """
        CREATE INDEX i_protein_fragment
        ON interpro.webfront_protein (is_fragment)
        """
    )
    con.commit()
    cur.close()
    con.close()
    logger.info("done")


def populate_residues(uri: str, residues_file: str):
    logger.info("creating webfront_proteinresidue")

    con = connect(uri)
    cur = con.cursor()
    cur.execute("DROP TABLE IF EXISTS interpro.webfront_proteinresidue")
    cur.execute(
        """
        CREATE TABLE interpro.webfront_proteinresidue
        (
            residue_id SERIAL NOT NULL PRIMARY KEY,
            protein_acc VARCHAR(15) COLLATE "case_insensitive" NOT NULL,
            entry_acc VARCHAR(30) NOT NULL,
            entry_name VARCHAR(100),
            source_database VARCHAR(10) NOT NULL,
            description VARCHAR(255),
            fragments JSONB NOT NULL
        )
        """
    )

    query = """
        INSERT INTO interpro.webfront_proteinresidue (
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
    con = connect(uri)
    cur = con.cursor()
    cur.execute(
        """
        CREATE INDEX IF NOT EXISTS i_proteinresidue
        ON interpro.webfront_proteinresidue (UPPER(protein_acc))
        """
    )
    con.commit()
    cur.close()
    con.close()
