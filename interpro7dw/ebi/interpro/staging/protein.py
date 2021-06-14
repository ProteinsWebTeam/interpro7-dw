# -*- coding: utf-8 -*-

import gzip
from typing import Optional

import MySQLdb

from interpro7dw import logger
from interpro7dw.ebi import pdbe
from interpro7dw.ebi.interpro import production as ippro
from interpro7dw.ebi.interpro.utils import Table
from interpro7dw.utils import loadobj, merge_dumps, url2dict
from interpro7dw.utils import DirectoryTree, Store
from .utils import jsonify


def insert_isoforms(src_entries: str, pro_url: str, stg_url: str):
    entries = loadobj(src_entries)

    con = MySQLdb.connect(**url2dict(stg_url))
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
    cur.close()

    sql = """
        INSERT INTO webfront_varsplic VALUES (%s, %s, %s, %s, %s)
    """
    with Table(con, sql) as table:
        for accession, variant in ippro.get_isoforms(pro_url).items():
            features = {}
            for entry_acc, locations in variant["matches"].items():
                entry = entries[entry_acc]

                features[entry_acc] = {
                    "accession": entry_acc,
                    "integrated": entry.integrated_in,
                    "name": entry.name,
                    "type": entry.type.lower(),
                    "source_database": entry.database,
                    "locations": locations
                }

            table.insert((
                accession,
                variant["protein_acc"],
                variant["length"],
                variant["sequence"],
                jsonify(features)
            ))

    con.commit()

    cur = con.cursor()
    cur.execute(
        """
        CREATE INDEX i_varsplic
        ON webfront_varsplic (protein_acc)
        """
    )
    cur.close()
    con.close()


def insert_proteins(p_entries: str, p_proteins: str, p_structures: str,
                    p_taxonomy: str, p_uniprot2comments: str,
                    p_uniprot2name: str, p_uniprot2evidences: str,
                    p_uniprot2ida: str, p_uniprot2matches: str,
                    p_uniprot2proteome: str, p_uniprot2sequence: str,
                    pro_url: str, stg_url: str):
    logger.info("loading CATH/SCOP domains")
    uniprot2cath = pdbe.get_cath_domains(pro_url)
    uniprot2scop = pdbe.get_scop_domains(pro_url)

    logger.info("preparing data")
    proteins = Store(p_proteins)
    u2comments = Store(p_uniprot2comments)
    u2descriptions = Store(p_uniprot2name)
    u2evidences = Store(p_uniprot2evidences)
    u2ida = Store(p_uniprot2ida)
    u2matches = Store(p_uniprot2matches)
    u2proteome = Store(p_uniprot2proteome)
    u2sequence = Store(p_uniprot2sequence)

    taxonomy = {}
    for taxid, info in loadobj(p_taxonomy).items():
        taxonomy[taxid] = jsonify({
            "taxId": taxid,
            "scientificName": info["sci_name"],
            "fullName": info["full_name"]
        })

    uniprot2pdbe = {}
    for pdb_id, entry in loadobj(p_structures).items():
        for uniprot_acc in entry["proteins"]:
            try:
                uniprot2pdbe[uniprot_acc].append(pdb_id)
            except KeyError:
                uniprot2pdbe[uniprot_acc] = [pdb_id]

    logger.info("counting proteins/IDA")
    ida_count = {}
    for dom_members, dom_arch, dom_arch_id in u2ida.values():
        try:
            ida_count[dom_arch_id] += 1
        except KeyError:
            ida_count[dom_arch_id] = 1

    logger.info("inserting proteins")
    entries = loadobj(p_entries)
    con = MySQLdb.connect(**url2dict(stg_url), charset="utf8mb4")
    cur = con.cursor()
    cur.execute(
        """
        SELECT protein_acc, COUNT(*)
        FROM webfront_varsplic
        GROUP BY protein_acc
        """
    )
    isoforms = dict(cur.fetchall())

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

    cur.execute(
        """
        SELECT protein_acc, predictor, COUNT(*)
        FROM webfront_structuralmodel
        WHERE protein_acc IS NOT NULL
        GROUP BY protein_acc, predictor
        """
    )
    u2structmodels = {}
    for uniprot_acc, predictor, cnt in cur:
        try:
            u2structmodels[uniprot_acc][predictor] = cnt
        except KeyError:
            u2structmodels[uniprot_acc] = {predictor: cnt}

    cur.close()

    i = 0
    sql = """
        INSERT into webfront_protein
        VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)
    """
    with Table(con, sql) as table:
        for uniprot_acc, protein_info in proteins.items():
            taxid = protein_info["taxid"]

            try:
                taxon = taxonomy[taxid]
            except KeyError:
                table.close()
                con.close()
                raise RuntimeError(f"{uniprot_acc}: invalid taxon {taxid}")

            try:
                name = u2descriptions[uniprot_acc]
            except KeyError:
                table.close()
                con.close()
                raise RuntimeError(f"{uniprot_acc}: missing name")

            try:
                evidence, gene = u2evidences[uniprot_acc]
            except KeyError:
                table.close()
                con.close()
                raise RuntimeError(f"{uniprot_acc}: missing evidence")

            try:
                sequence = u2sequence[uniprot_acc]
            except KeyError:
                table.close()
                con.close()
                raise RuntimeError(f"{uniprot_acc}: missing sequence")

            proteome_id = u2proteome.get(uniprot_acc)

            clans = []
            databases = {}
            go_terms = {}
            for entry_acc in u2matches.get(uniprot_acc, []):
                entry = entries[entry_acc]

                try:
                    databases[entry.database] += 1
                except KeyError:
                    databases[entry.database] = 1

                if entry.clan:
                    clans.append(entry.clan["accession"])

                for term in entry.go_terms:
                    go_terms[term["identifier"]] = term

            protein_structures = {}
            domains = uniprot2cath.get(uniprot_acc)
            if domains:
                protein_structures["cath"] = {}

                for dom in domains.values():
                    dom_id = dom["id"]

                    protein_structures["cath"][dom_id] = {
                        "domain_id": dom["superfamily"]["id"],
                        "coordinates": dom["locations"]
                    }

            domains = uniprot2scop.get(uniprot_acc)
            if domains:
                protein_structures["scop"] = {}

                for dom in domains.values():
                    dom_id = dom["id"]

                    protein_structures["scop"][dom_id] = {
                        "domain_id": dom["superfamily"]["id"],
                        "coordinates": dom["locations"]
                    }

            try:
                dom_members, dom_arch, dom_arch_id = u2ida[uniprot_acc]
            except KeyError:
                dom_arch = dom_arch_id = None
                dom_count = 0
            else:
                dom_count = ida_count[dom_arch_id]

            table.insert((
                uniprot_acc,
                protein_info["identifier"],
                taxon,
                name,
                jsonify(u2comments.get(uniprot_acc)),
                gzip.compress(sequence.encode("utf-8")),
                protein_info["length"],
                proteome_id,
                gene,
                jsonify(list(go_terms.values())),
                evidence,
                "reviewed" if protein_info["reviewed"] else "unreviewed",
                1 if protein_info["fragment"] else 0,
                jsonify(protein_structures),
                protein_info["taxid"],
                dom_arch_id,
                dom_arch,
                jsonify({
                    "domain_architectures": dom_count,
                    "entries": databases,
                    "isoforms": isoforms.get(uniprot_acc, 0),
                    "proteomes": 1 if proteome_id else 0,
                    "sets": len(set(clans)),
                    "structures": len(uniprot2pdbe.get(uniprot_acc, [])),
                    "structural_models": u2structmodels.get(uniprot_acc, {}),
                    "taxa": 1
                })
            ))

            i += 1
            if not i % 10000000:
                logger.info(f"{i:>12,}")

        logger.info(f"{i:>12,}")

    con.commit()

    proteins.close()
    u2comments.close()
    u2descriptions.close()
    u2evidences.close()
    u2ida.close()
    u2matches.close()
    u2proteome.close()
    u2sequence.close()

    logger.info("indexing")
    cur = con.cursor()
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

    logger.info("complete")


def insert_extra_features(stg_url: str, p_uniprot2features: str):
    logger.info("starting")

    con = MySQLdb.connect(**url2dict(stg_url), charset="utf8mb4")
    cur = con.cursor()
    cur.execute("DROP TABLE IF EXISTS webfront_proteinfeature")
    cur.execute(
        """
        CREATE TABLE webfront_proteinfeature
        (
            feature_id INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
            protein_acc VARCHAR(15) NOT NULL,
            entry_acc VARCHAR(25) NOT NULL,
            source_database VARCHAR(10) NOT NULL,
            location_start INT NOT NULL,
            location_end INT NOT NULL,
            sequence_feature VARCHAR(35)
        ) CHARSET=utf8mb4 DEFAULT COLLATE=utf8mb4_unicode_ci
        """
    )
    cur.close()

    sql = """
        INSERT INTO webfront_proteinfeature (
          protein_acc, entry_acc, source_database, location_start,
          location_end, sequence_feature
        )
        VALUES (%s, %s, %s, %s, %s, %s)
    """
    with Store(p_uniprot2features) as proteins, Table(con, sql) as table:
        i = 0
        for uniprot_acc, entries in proteins.items():
            for entry_acc, info in entries.items():
                for pos_start, pos_end, seq_feature in info["locations"]:
                    table.insert((
                        uniprot_acc,
                        entry_acc,
                        info["database"],
                        pos_start,
                        pos_end,
                        seq_feature
                    ))

            i += 1
            if not i % 10000000:
                logger.info(f"{i:>12,}")

        logger.info(f"{i:>12,}")
    con.commit()

    logger.info("indexing")
    cur = con.cursor()
    cur.execute(
        """
        CREATE INDEX i_proteinfeature
        ON webfront_proteinfeature (protein_acc)
        """
    )
    cur.close()
    con.close()
    logger.info("complete")


def insert_residues(pro_url: str, stg_url: str, tmpdir: Optional[str] = None):
    dt = DirectoryTree(root=tmpdir)

    logger.info("exporting residues")
    files = ippro.export_residues(pro_url, dt)

    logger.info("inserting residues")
    con = MySQLdb.connect(**url2dict(stg_url), charset="utf8mb4")
    cur = con.cursor()
    cur.execute("DROP TABLE IF EXISTS webfront_proteinresidue")
    cur.execute(
        """
        CREATE TABLE webfront_proteinresidue
        (
            residue_id INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
            protein_acc VARCHAR(15) NOT NULL,
            entry_acc VARCHAR(25) NOT NULL,
            entry_name VARCHAR(100),
            source_database VARCHAR(10) NOT NULL,
            description VARCHAR(255),
            fragments LONGTEXT NOT NULL
        ) CHARSET=utf8mb4 DEFAULT COLLATE=utf8mb4_unicode_ci
        """
    )
    cur.close()

    sql = """
        INSERT INTO webfront_proteinresidue (
          protein_acc, entry_acc, entry_name, source_database, description,
          fragments
        )
        VALUES (%s, %s, %s, %s, %s, %s)
    """
    with Table(con, sql) as table:
        i = 0
        for protein_acc, entries in merge_dumps(files, replace=True):
            for entry_acc, entry in entries.items():
                for descr, locations in entry["descriptions"].items():
                    locations.sort(key=lambda x: (x[1], x[2]))
                    table.insert((
                        protein_acc,
                        entry_acc,
                        entry["name"],
                        entry["database"],
                        descr,
                        jsonify(locations, nullable=False)
                    ))

            i += 1
            if not i % 10000000:
                logger.info(f"{i:>15,}")

        logger.info(f"{i:>15,}")
    con.commit()

    logger.info(f"temporary files: {dt.size / 1024 ** 2:.0f} MB")
    dt.remove()

    logger.info("indexing")
    cur = con.cursor()
    cur.execute(
        """
        CREATE INDEX i_proteinresidue
        ON webfront_proteinresidue (protein_acc)
        """
    )
    cur.close()
    con.close()
    logger.info("complete")
