# -*- coding: utf-8 -*-

import hashlib
import json
from typing import Optional

import MySQLdb

from interpro7dw import logger
from interpro7dw.ebi import pdbe
from interpro7dw.ebi.interpro import production as ippro
from interpro7dw.ebi.interpro.utils import Table, repr_fragment
from interpro7dw.utils import dataload, url2dict, Store


def export_ida(src_entries: str, src_matches: str, dst_ida: str,
               dir: Optional[str]=None, processes: int=1):
    logger.info("starting")
    pfam2interpro = {}
    for entry in dataload(src_entries).values():
        if entry.database == "pfam":
            pfam2interpro[entry.accession] = entry.integrated_in

    with Store(src_matches) as src, Store(dst_ida, src.chunks, dir) as dst:
        i = 0
        for protein_acc, entries in src.items():
            all_locations = []
            for entry_acc, locations in entries.items():
                try:
                    interpro_acc = pfam2interpro[entry_acc]
                except KeyError:
                    # Not a Pfam signature
                    continue

                for loc in locations:
                    all_locations.append({
                        "pfam": entry_acc,
                        "interpro": interpro_acc,
                        # We do not consider fragmented locations
                        "start": loc["fragments"][0]["start"],
                        "end": max(f["end"] for f in loc["fragments"])
                    })

            domains = []
            for loc in sorted(all_locations, key=repr_fragment):
                if loc["interpro"]:
                    domains.append(f"{loc['pfam']}:{loc['interpro']}")
                else:
                    domains.append(loc["pfam"])

            dom_arch = '-'.join(domains)
            dom_arch_id = hashlib.sha1(dom_arch.encode("utf-8")).hexdigest()
            dst[protein_acc] = (dom_arch, dom_arch_id)

            i += 1
            if not i % 1000000:
                dst.sync()

                if not i % 10000000:
                    logger.info(f"{i:>12,}")

        logger.info(f"{i:>12,}")
        size = dst.merge(processes=processes)
        logger.info(f"temporary files: {size/1024/1024:.0f} MB")


def insert_isoforms(src_entries: str, pro_url: str, stg_url: str):
    entries = dataload(src_entries)

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
            features LONGTEXT NOT NULL
        ) CHARSET=utf8 DEFAULT COLLATE=utf8_unicode_ci
        """
    )
    cur.close()

    sql = """
        INSERT INTO webfront_varsplic VALUES (%s, %s, %s, %s, %s)
    """
    with Table(con, sql) as table:
        for obj in ippro.get_isoforms(pro_url):
            accession = obj[0]
            protein_acc = obj[1]
            length = obj[2]
            sequence = obj[3]
            features = obj[4]

            enriched_features = {}
            for entry_acc, locations in features.items():
                entry = entries[entry_acc]

                enriched_features[entry_acc] = {
                    "accession": entry_acc,
                    "integrated": entry.integrated_in,
                    "name": entry.name,
                    "type": entry.type,
                    "source_database": entry.database,
                    "locations": locations
                }

            table.insert((accession, protein_acc, length, sequence,
                          json.dumps(enriched_features)))

    con.commit()

    cur = con.cursor()
    cur.execute("CREATE INDEX i_varsplic "
                "ON webfront_varsplic (protein_acc)")
    cur.close()
    con.close()


def jsonify(obj, allow_na=True, fallback=None):
    if obj or allow_na:
        return json.dumps(obj)
    else:
        return fallback


def insert_proteins(src_proteins: str, src_comments: str,
                    src_descriptions: str, src_entries: str,
                    src_evidences: str, src_features: str, src_ida: str,
                    src_matches: str, src_proteomes: str, src_residues: str,
                    src_sequences: str, src_structures: str, src_taxonomy: str,
                    pro_url: str, stg_url: str):
    logger.info("loading CATH/SCOP domains")
    cath_domains = pdbe.get_cath_domains(pro_url)
    scop_domains = pdbe.get_scop_domains(pro_url)

    logger.info("preparing data")
    proteins = Store(src_proteins)
    comments = Store(src_comments)
    descriptions = Store(src_descriptions)
    evidences = Store(src_evidences)
    features = Store(src_features)
    ida = Store(src_ida)
    matches = Store(src_matches)
    proteomes = Store(src_proteomes)
    residues = Store(src_residues)
    sequences = Store(src_sequences)

    entries = {}
    for entry in dataload(src_entries).values():
        if entry.database == "interpro" and entry.go_terms:
            go_terms = entry.go_terms
        else:
            go_terms = []

        entries[entry.accession] = (entry.database, entry.clan, go_terms)

    taxonomy = {}
    for taxid, info in dataload(src_taxonomy):
        taxonomy[taxid] = json.dumps({
            "taxId": taxid,
            "scientificName": info["sci_name"],
            "fullName": info["full_name"]
        })

    uniprot2pdbe = {}
    for pdb_id, entry in dataload(src_structures).items():
        for uniprot_acc in entry["proteins"]:
            try:
                uniprot2pdbe[uniprot_acc].append(pdb_id)
            except KeyError:
                uniprot2pdbe[uniprot_acc] = [pdb_id]

    logger.info("counting proteins/IDA")
    ida_count = {}
    for dom_arch, dom_arch_id in ida.values():
        try:
            ida_count[dom_arch_id] += 1
        except KeyError:
            ida_count[dom_arch_id] = 1

    logger.info("inserting proteins")
    con = MySQLdb.connect(**url2dict(stg_url))
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
            sequence LONGTEXT NOT NULL,
            length INT(11) NOT NULL,
            proteome VARCHAR(20),
            gene VARCHAR(70),
            go_terms LONGTEXT,
            evidence_code INT(11) NOT NULL,
            source_database VARCHAR(10) NOT NULL,
            residues LONGTEXT NOT NULL,
            is_fragment TINYINT NOT NULL,
            structure LONGTEXT NOT NULL,
            tax_id VARCHAR(20) NOT NULL,
            extra_features LONGTEXT NOT NULL,
            ida_id VARCHAR(40),
            ida TEXT,
            counts LONGTEXT NOT NULL
        ) CHARSET=utf8 DEFAULT COLLATE=utf8_unicode_ci
        """
    )
    cur.close()

    i = 0
    sql = """
        INSERT into webfront_protein
        VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)
    """
    with Table(con, sql) as table:
        for accession, info in proteins.items():
            taxid = info["taxid"]

            try:
                organism = taxonomy[taxid]
            except KeyError:
                table.close()
                con.close()
                raise RuntimeError(f"{accession}: invalid taxon {taxid}")

            try:
                name = descriptions[accession]
            except KeyError:
                table.close()
                con.close()
                raise RuntimeError(f"{accession}: missing name")

            try:
                evidence, gene = evidences[accession]
            except KeyError:
                table.close()
                con.close()
                raise RuntimeError(f"{accession}: missing evidence")

            try:
                sequence = sequences[accession]
            except KeyError:
                table.close()
                con.close()
                raise RuntimeError(f"{accession}: missing sequence")

            proteome_id = proteomes.get(accession)

            clans = []
            databases = {}
            go_terms = {}
            for entry_acc in matches.get(accession, []):
                database, clan, terms = entries[entry_acc]

                try:
                    databases[database] += 1
                except KeyError:
                    databases[database] = 1

                if clan:
                    clans.append(clan)

                for term in terms:
                    go_terms[term["identifier"]] = term

            structures = uniprot2pdbe.get(accession, [])
            extra_features = {}
            for pdb_id in structures:
                if pdb_id in cath_domains:
                    domain = cath_domains[pdb_id]
                    try:
                        extra_features["cath"].update(domain)
                    except KeyError:
                        extra_features["cath"] = domain.copy()

                if pdb_id in scop_domains:
                    domain = scop_domains[pdb_id]
                    try:
                        extra_features["scop"].update(domain)
                    except KeyError:
                        extra_features["scop"] = domain.copy()

            try:
                dom_arch, dom_arch_id = ida[accession]
            except KeyError:
                dom_arch = dom_arch_id = None
                dom_count = 0
            else:
                dom_count = ida_count[accession]

            table.insert((
                accession,
                info["identifier"],
                organism,
                name,
                jsonify(comments.get(accession), allow_na=False),
                sequence,
                info["length"],
                proteome_id,
                gene,
                jsonify(list(go_terms.values()), allow_na=False),
                evidence,
                "reviewed" if info["reviewed"] else "unreviewed",
                jsonify(residues.get(accession), allow_na=False),
                1 if info["fragment"] else 0,
                jsonify(extra_features, allow_na=False),
                info["taxid"],
                jsonify(features.get(accession), allow_na=False),
                dom_arch_id,
                dom_arch,
                jsonify({
                    "entries": databases,
                    "structures": len(structures),
                    "sets": len(set(clans)),
                    "proteomes": 1 if proteome_id else 0,
                    "taxa": 1,
                    "idas": dom_count,
                    "isoforms": isoforms.get(accession, 0)
                })
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
        CREATE UNIQUE INDEX ui_protein_identifier
        ON webfront_protein (identifier)
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