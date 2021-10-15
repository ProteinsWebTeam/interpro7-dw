import MySQLdb

from interpro7dw import pdbe
from interpro7dw.utils import logger
from interpro7dw.utils.mysql import url2dict
from interpro7dw.utils.store import SimpleStore, Store, loadobj

from .utils import jsonify


def insert_proteins(ipr_url: str, pdbe_url: str, entries_file: str,
                    isoforms_file: str, structures_file: str, taxa_file: str,
                    proteins_file: str, domorgs_file: str, evidences_file: str,
                    functions_file: str, matches_file: str, names_file: str,
                    proteomes_file: str, sequences_file: str):
    """Creates and populates the MySQL webfront_protein table.

    :param ipr_url: InterPro MySQL connection string.
    :param pdbe_url: PDBe Oracle connection string.
    :param entries_file: File of InterPro entries
        and member database signatures.
    :param isoforms_file: SimpleStore file of protein isoforms.
    :param structures_file: File of PDBe structures.
    :param taxa_file: File of taxonomic information.
    :param proteins_file: Store file of proteins.
    :param domorgs_file: Store file of domain organisations.
    :param evidences_file: Store file of protein evidences/genes.
    :param functions_file: Store file of protein functions.
    :param matches_file: Store file of protein matches.
    :param names_file: Store file of protein descriptions/names.
    :param proteomes_file: Store file of protein-proteome mapping.
    :param sequences_file: Store file of protein sequences.
    """
    logger.info("loading CATH/SCOP domains")
    protein2cath = pdbe.get_cath_domains(pdbe_url)
    protein2scop = pdbe.get_scop_domains(pdbe_url)

    logger.info("loading entries")
    entries = loadobj(entries_file)

    logger.info("loading isoforms")
    protein2isoforms = {}
    with SimpleStore(isoforms_file) as store:
        for isoform in store:
            protein_acc = isoform["protein"]
            try:
                protein2isoforms[protein_acc] += 1
            except KeyError:
                protein2isoforms[protein_acc] = 1

    logger.info("loading PDBe structures")
    protein2structures = {}
    for pdbe_id, entry in loadobj(structures_file).items():
        for protein_acc in entry["proteins"]:
            try:
                protein2structures[protein_acc] += 1
            except KeyError:
                protein2structures[protein_acc] = 1

    logger.info("loading taxa")
    taxa = {}
    for taxon_id, taxon in loadobj(taxa_file).items():
        taxa[taxon_id] = jsonify({
            "taxId": taxon_id,
            "scientificName": taxon["sci_name"],
            "fullName": taxon["full_name"]
        })

    logger.info("creating webfront_protein")
    proteins = Store(proteins_file)
    functions = Store(functions_file)
    names = Store(names_file)
    evidences = Store(evidences_file)
    domorgs = Store(domorgs_file)
    matches = Store(matches_file)
    proteomes = Store(proteomes_file)
    sequences = Store(sequences_file)

    con = MySQLdb.connect(**url2dict(ipr_url), charset="utf8mb4")
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
            counts LONGTEXT NOT NULL
        ) CHARSET=utf8mb4 DEFAULT COLLATE=utf8mb4_unicode_ci
        """
    )

    query = """
        INSERT into webfront_protein
        VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)
    """
    args = []
    i = 0
    for i, (protein_acc, protein) in enumerate(proteins.items()):
        taxon_id = protein["taxid"]

        try:
            taxon = taxa[taxon_id]
        except KeyError:
            cur.close()
            con.close()
            raise KeyError(f"{protein_acc}: invalid taxon {taxon_id}")

        try:
            name = names[protein_acc]
        except KeyError:
            cur.close()
            con.close()
            raise KeyError(f"{protein_acc}: missing name")

        try:
            evidence, gene = evidences[protein_acc]
        except KeyError:
            cur.close()
            con.close()
            raise RuntimeError(f"{protein_acc}: missing evidence/gene")

        try:
            gzipped_sequence = sequences[protein_acc]
        except KeyError:
            cur.close()
            con.close()
            raise RuntimeError(f"{protein_acc}: missing sequence")

        proteome_id = proteomes.get(protein_acc)

        clans = []
        databases = {}
        go_terms = {}
        for entry_acc in matches.get(protein_acc, []):
            entry = entries[entry_acc]

            try:
                databases[entry.database] += 1
            except KeyError:
                databases[entry.database] = 1

            if entry.clan:
                clans.append(entry.clan["accession"])

            for term in entry.go_terms:
                go_terms[term["identifier"]] = term

        # Adds CATH/SCOP structures
        protein_structures = {}
        for key, obj in [("cath", protein2cath), ("scop", protein2scop)]:
            domains = obj.get(protein_acc)
            if domains:
                protein_structures[key] = {}

                for dom in domains.values():
                    dom_id = dom["id"]

                    protein_structures[key][dom_id] = {
                        "domain_id": dom["superfamily"]["id"],
                        "coordinates": dom["locations"]
                    }

        domains = protein2scop.get(protein_acc)
        if domains:
            protein_structures["cath"] = {}

            for dom in domains.values():
                dom_id = dom["id"]

                protein_structures["cath"][dom_id] = {
                    "domain_id": dom["superfamily"]["id"],
                    "coordinates": dom["locations"]
                }

        try:
            dom_id, dom_str, _, _, dom_proteins = domorgs[protein_acc]
        except KeyError:
            dom_id = dom_str = None
            dom_proteins = 0

        args.append((
            protein_acc,
            protein["identifier"],
            taxon,
            name,
            jsonify(functions.get(protein_acc), nullable=True),
            gzipped_sequence,
            protein["length"],
            proteome_id,
            gene,
            jsonify(list(go_terms.values()), nullable=True),
            evidence,
            "reviewed" if protein["reviewed"] else "unreviewed",
            1 if protein["fragment"] else 0,
            jsonify(protein_structures, nullable=True),
            taxon_id,
            dom_id,
            dom_str,
            jsonify({
                "domain_architectures": dom_proteins,
                "entries": databases,
                "isoforms": protein2isoforms.get(protein_acc, 0),
                "proteomes": 1 if proteome_id else 0,
                "sets": len(set(clans)),
                "structures": protein2structures.get(protein_acc, 0),
                "taxa": 1
            })
        ))

        if (i + 1) % 1e3 == 0:
            cur.executemany(query, args)
            args.clear()

            if (i + 1) % 10e6 == 0:
                logger.info(f"{i + 1:>15,}")

    if args:
        cur.executemany(query, args)
        args.clear()

    con.commit()
    logger.info(f"{i + 1:>15,}")

    proteins.close()
    functions.close()
    names.close()
    evidences.close()
    domorgs.close()
    matches.close()
    proteomes.close()
    sequences.close()

    logger.info("creating indexes")
    cur.execute(
        """
        CREATE UNIQUE INDEX ui_protein_identifier
        ON webfront_protein (identifier)
        """
    )
    cur.execute(
        """
        CREATE INDEX i_protein_proteome
        ON webfront_protein (proteome)
        """
    )
    cur.execute(
        """
        CREATE INDEX i_protein_database
        ON webfront_protein (source_database)
        """
    )
    cur.execute(
        """
        CREATE INDEX i_protein_taxon
        ON webfront_protein (tax_id)
        """
    )
    cur.execute(
        """
        CREATE INDEX i_protein_ida
        ON webfront_protein (ida_id)
        """
    )
    cur.execute(
        """
        CREATE INDEX i_protein_fragment
        ON webfront_protein (is_fragment)
        """
    )
    cur.close()
    con.close()

    logger.info("done")
