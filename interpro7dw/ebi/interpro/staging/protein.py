# -*- coding: utf-8 -*-

import hashlib
from typing import Optional

import MySQLdb

from interpro7dw import logger
from interpro7dw.ebi import pdbe
from interpro7dw.ebi.interpro import production as ippro
from interpro7dw.ebi.interpro.utils import Table, repr_fragment
from interpro7dw.utils import dataload, url2dict, Store
from .utils import jsonify


def export_ida(src_entries: str, src_matches: str, dst_ida: str,
               dir: Optional[str]=None, processes: int=1):
    logger.info("starting")
    pfam2interpro = {}
    for entry in dataload(src_entries).values():
        if entry.database == "pfam":
            pfam2interpro[entry.accession] = entry.integrated_in

    with Store(src_matches) as src, Store(dst_ida, src.get_keys(), dir) as dst:
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


def export_uniprot2entries(p_entries: str, p_uniprot2matches: str, output: str,
                           dir: Optional[str]=None, processes: int=1):
    logger.info("starting")
    entries = {}
    for entry in dataload(p_entries).values():
        if entry.database == "interpro" and entry.go_terms:
            go_terms = entry.go_terms
        else:
            go_terms = []

        entries[entry.accession] = (
            entry.accession,
            entry.database,
            entry.clan["accession"] if entry.clan else None,
            go_terms
        )

    i = 0
    uniprot2matches = Store(p_uniprot2matches)
    uniprot2entries = Store(output, uniprot2matches.get_keys(), dir)
    for uniprot_acc, matches in uniprot2matches.items():
        _entries = []
        for entry_acc in matches:
            _entries.append(entries[entry_acc])

        uniprot2entries[uniprot_acc] = _entries

        i += 1
        if not i % 1000000:
            uniprot2entries.sync()

            if not i % 10000000:
                logger.info(f"{i:>12,}")

    logger.info(f"{i:>12,}")
    uniprot2matches.close()

    size = uniprot2entries.merge(processes=processes)
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
            features LONGTEXT
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
                          jsonify(enriched_features)))

    con.commit()

    cur = con.cursor()
    cur.execute("CREATE INDEX i_varsplic "
                "ON webfront_varsplic (protein_acc)")
    cur.close()
    con.close()


def insert_proteins(p_proteins: str, p_structures: str, p_taxonomy: str,
                    p_uniprot2comments: str, p_uniprot2name: str,
                    p_uniprot2entries: str, p_uniprot2evidences: str,
                    p_uniprot2features: str, p_uniprot2ida: str,
                    p_uniprot2proteome: str, p_uniprot2residues: str,
                    p_uniprot2sequence: str, pro_url: str, stg_url: str):
    logger.info("loading CATH/SCOP domains")
    cath_domains = pdbe.get_cath_domains(pro_url)
    scop_domains = pdbe.get_scop_domains(pro_url)

    logger.info("preparing data")
    proteins = Store(p_proteins)
    u2comments = Store(p_uniprot2comments)
    u2descriptions = Store(p_uniprot2name)
    u2entries = Store(p_uniprot2entries)
    u2evidences = Store(p_uniprot2evidences)
    u2features = Store(p_uniprot2features)
    u2ida = Store(p_uniprot2ida)
    u2proteome = Store(p_uniprot2proteome)
    u2residues = Store(p_uniprot2residues)
    u2sequence = Store(p_uniprot2sequence)

    taxonomy = {}
    for taxid, info in dataload(p_taxonomy).items():
        taxonomy[taxid] = jsonify({
            "taxId": taxid,
            "scientificName": info["sci_name"],
            "fullName": info["full_name"]
        })

    uniprot2pdbe = {}
    for pdb_id, entry in dataload(p_structures).items():
        for uniprot_acc in entry["proteins"]:
            try:
                uniprot2pdbe[uniprot_acc].append(pdb_id)
            except KeyError:
                uniprot2pdbe[uniprot_acc] = [pdb_id]

    logger.info("counting proteins/IDA")
    ida_count = {}
    for dom_arch, dom_arch_id in u2ida.values():
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
            residues LONGTEXT,
            is_fragment TINYINT NOT NULL,
            structure LONGTEXT,
            tax_id VARCHAR(20) NOT NULL,
            extra_features LONGTEXT,
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
                taxon = taxonomy[taxid]
            except KeyError:
                table.close()
                con.close()
                raise RuntimeError(f"{accession}: invalid taxon {taxid}")

            try:
                name = u2descriptions[accession]
            except KeyError:
                table.close()
                con.close()
                raise RuntimeError(f"{accession}: missing name")

            try:
                evidence, gene = u2evidences[accession]
            except KeyError:
                table.close()
                con.close()
                raise RuntimeError(f"{accession}: missing evidence")

            try:
                sequence = u2sequence[accession]
            except KeyError:
                table.close()
                con.close()
                raise RuntimeError(f"{accession}: missing sequence")

            proteome_id = u2proteome.get(accession)

            clans = []
            databases = {}
            go_terms = {}
            entries = u2entries.get(accession, [])
            for entry_acc, database, clan, terms in entries:
                try:
                    databases[database] += 1
                except KeyError:
                    databases[database] = 1

                if clan:
                    clans.append(clan["accession"])

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
                dom_arch, dom_arch_id = u2ida[accession]
            except KeyError:
                dom_arch = dom_arch_id = None
                dom_count = 0
            else:
                dom_count = ida_count[dom_arch_id]

            table.insert((
                accession,
                info["identifier"],
                taxon,
                name,
                jsonify(u2comments.get(accession)),
                sequence,
                info["length"],
                proteome_id,
                gene,
                jsonify(list(go_terms.values())),
                evidence,
                "reviewed" if info["reviewed"] else "unreviewed",
                jsonify(u2residues.get(accession)),
                1 if info["fragment"] else 0,
                jsonify(extra_features),
                info["taxid"],
                jsonify(u2features.get(accession)),
                dom_arch_id,
                dom_arch,
                jsonify({
                    "domain_architectures": dom_count,
                    "entries": databases,
                    "isoforms": isoforms.get(accession, 0),
                    "proteomes": 1 if proteome_id else 0,
                    "sets": len(set(clans)),
                    "structures": len(structures),
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
    u2entries.close()
    u2evidences.close()
    u2features.close()
    u2ida.close()
    u2proteome.close()
    u2residues.close()
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
