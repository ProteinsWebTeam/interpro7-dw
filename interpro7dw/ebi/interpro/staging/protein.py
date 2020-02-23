# -*- coding: utf-8 -*-

import hashlib
import json
from typing import Optional

import MySQLdb

from interpro7dw import logger
from interpro7dw.ebi.interpro import production as ippro
from interpro7dw.ebi.interpro.utils import Table, repr_fragment
from interpro7dw.utils import dataload, url2dict, Store


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


def insert_proteins(src_proteins: str, src_comments: str, src_descriptions: str,
                    src_evidences: str, src_proteomes: str, src_residues: str,
                    src_sequences: str, src_taxonomy: str,
                    url: str):
    con = MySQLdb.connect(**url2dict(url))
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
            description LONGTEXT NOT NULL,
            sequence LONGTEXT NOT NULL,
            length INT(11) NOT NULL,
            proteome VARCHAR(20),
            gene VARCHAR(70),
            go_terms LONGTEXT NOT NULL,
            evidence_code INT(11) NOT NULL,
            source_database VARCHAR(10) NOT NULL,
            residues LONGTEXT NOT NULL,
            is_fragment TINYINT NOT NULL,
            structure LONGTEXT NOT NULL,
            tax_id VARCHAR(20) NOT NULL,
            extra_features LONGTEXT NOT NULL,
            ida_id VARCHAR(40) DEFAULT NULL,
            ida TEXT DEFAULT NULL,
            counts LONGTEXT DEFAULT NULL
        ) CHARSET=utf8 DEFAULT COLLATE=utf8_unicode_ci
        """
    )
    cur.close()

    proteins = Store(src_proteins)
    comments = Store(src_comments)
    descriptions = Store(src_descriptions)
    evidences = Store(src_evidences)
    proteomes = Store(src_proteomes)
    residues = Store(src_residues)
    sequences = Store(src_sequences)

    taxonomy = {}
    for taxid, info in dataload(src_taxonomy):
        taxonomy[taxid] = json.dumps({
            "taxId": taxid,
            "scientificName": info["sci_name"],
            "fullName": info["full_name"]
        })

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

            table.insert((
                accession,              # accession
                info["identifier"],     # identifier
                organism,   # organism
                name,   # name
                json.dumps(comments.get(accession, [])),   # description
                sequence,   # sequence
                info["length"],   # length
                proteomes.get(accession),   # proteome
                gene,   # gene
                json.dumps([]),   # go_terms
                evidence,   # evidence_code
                "reviewed" if info["reviewed"] else "unreviewed",  # source_database
                json.dumps(residues.get(protein_acc, {})),,   # residues
                1 if info["fragment"] else 0,   # is_fragment
                None,   # structure
                info["taxid"],   # tax_id
                None,   # extra_features
                None,   # ida_id
                None,   # ida
                None    # counts
            ))

    con.commit()
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
