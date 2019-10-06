# -*- coding: utf-8 -*-

import json
from datetime import datetime

import MySQLdb

from i7dw import io, logger
from i7dw.interpro import mysql
from .utils import parse_url


def make_release_notes(stg_url: str, rel_url: str, src_proteins: str,
                       src_matches: str, src_proteomes: str, version: str,
                       release_date: str):
    logger.info("loading data")
    # PDBe data (and UniProt-PDBe mappings)
    protein_structures = {}
    pdbe_ids = set()
    pdbe_date = None
    for s in mysql.structures.iter_structures(stg_url):
        pdbe_id = s["accession"]
        for protein_acc in s["proteins"]:
            try:
                protein_structures[protein_acc].add(pdbe_id)
            except KeyError:
                protein_structures[protein_acc] = {pdbe_id}

        pdbe_ids.add(pdbe_id)
        if pdbe_date is None or s["date"] > pdbe_date:
            pdbe_date = s["date"]

    # UniProt proteomes
    proteomes = {p["accession"] for p in mysql.proteomes.iter_proteomes(stg_url)}

    # UniProt taxonomy
    taxa = {t["id"] for t in mysql.taxonomy.iter_taxa(stg_url, lineage=False)}

    # Integrated member database signatures
    integrated = set()
    for entry in mysql.entries.get_entries(stg_url).values():
        if entry["integrated"]:
            integrated.add(entry["accession"])

    # Counting member database sets
    database_sets = {}
    for s in mysql.entries.iter_sets(stg_url):
        try:
            database_sets[s["database"]] += 1
        except KeyError:
            database_sets[s["database"]] = 1

    # Get UniProtKB version
    con = MySQLdb.connect(**parse_url(stg_url), charset="utf8")
    cur = con.cursor()
    uniprot = {}
    cur.execute(
        """
        SELECT name_long, version
        FROM webfront_database
        WHERE type='protein'
        """
    )
    for name, version in cur:
        uniprot[name] = {
            "version": version,
            "count": 0,
            "signatures": 0,
            "integrated_signatures": 0
        }

    cur.close()
    con.close()

    logger.info("gathering UniProtKB statistics")
    proteins = io.Store(src_proteins)
    protein_matches = io.Store(src_matches)
    protein_proteome = io.Store(src_proteomes)

    interpro_structures = set()
    interpro_proteomes = set()
    interpro_taxa = set()

    n_proteins = 0
    for protein_acc, protein_info in proteins:
        if protein_info["is_reviewed"]:
            k = "UniProtKB/Swiss-Prot"
        else:
            k = "UniProtKB/TrEMBL"

        uniprot[k]["count"] += 1
        matches = protein_matches.get(protein_acc)
        if matches:
            # Protein has a least one signature
            uniprot[k]["signatures"] += 1

            # Search if the protein has at least one integrated signature
            for signature_acc in matches:
                if signature_acc in integrated:
                    # It has!
                    uniprot[k]["integrated_signatures"] += 1

                    # Add taxon, proteome, and structures
                    interpro_taxa.add(protein_info["taxon"])

                    upid = protein_proteome.get(protein_acc)
                    if upid:
                        interpro_proteomes.add(upid)

                    interpro_structures |= protein_structures.get(protein_acc, set())
                    break

        n_proteins += 1
        if not n_proteins % 10000000:
            logger.info(f"{n_proteins:>12,}")

    proteins.close()
    protein_matches.close()
    protein_proteome.close()

    for k in ("count", "signatures", "integrated_signatures"):
        sp = uniprot["UniProtKB/Swiss-Prot"][k]
        tr = uniprot["UniProtKB/TrEMBL"][k]
        uniprot["UniProtKB"][k] = sp + tr

    logger.info(f"{n_proteins:>12,}")

    bad = interpro_structures - pdbe_ids
    if bad:
        raise RuntimeError(f"PDBe structure issues: {bad}")

    bad = interpro_proteomes - proteomes
    if bad:
        raise RuntimeError(f"UniProt proteome issues: {bad}")

    bad = interpro_taxa - taxa
    if bad:
        raise RuntimeError(f"UniProt taxonomy issues: {bad}")

    logger.info("tracking changes between releases")
    notes = {
        "interpro": {},
        "member_databases": [],
        "proteins": uniprot,
        "structures": {
            "total": len(pdbe_ids),
            "integrated": len(interpro_structures),
            "version": pdbe_date.strftime("%Y-%m-%d")
        },
        "proteomes": {
            "total": len(proteomes),
            "integrated": len(interpro_proteomes),
            "version": uniprot["UniProtKB"]["version"]
        },
        "taxonomy": {
            "total": len(taxa),
            "integrated": len(interpro_taxa),
            "version": uniprot["UniProtKB"]["version"]
        }
    }

    rel_entries = set()
    rel_integrated = set()
    for entry in mysql.entries.get_entries(rel_url).values():
        if entry["database"] == "interpro":
            rel_entries.add(entry["accession"])
        elif entry["integrated"]:
            # Signature already integrated in the previous release
            rel_integrated.add(entry["accession"])

    stg_entries = []
    new_entries = []
    for acc, entry in mysql.entries.get_entries(stg_url).items():
        stg_entries.append(entry)
        if entry["database"] == "interpro" and acc not in rel_entries:
            new_entries.append(acc)

    # Member database changes
    stg_dbs = mysql.databases.get_databases(stg_url)
    rel_dbs = mysql.databases.get_databases(rel_url)
    database_updates = set()
    new_databases = set()
    for name, info in stg_dbs.items():
        if name not in rel_dbs:
            new_databases.add(name)
        elif info["version"] != rel_dbs[name]["version"]:
            database_updates.add(name)

    member_databases = {}
    interpro_types = {}
    citations = set()
    n_interpro2go = 0
    latest_entry = None
    for entry in sorted(stg_entries, key=lambda x: x["accession"]):
        for pub in entry["citations"].values():
            if pub["PMID"] is not None:
                citations.add(pub["PMID"])

        acc = entry["accession"]
        database = entry["database"]
        entry_type = entry["type"]
        if database == "interpro":
            try:
                interpro_types[entry_type] += 1
            except KeyError:
                interpro_types[entry_type] = 1

            n_interpro2go += len(entry["go_terms"])
            latest_entry = acc
        else:
            try:
                obj = member_databases[database]
            except KeyError:
                obj = member_databases[database] = {
                    "name": stg_dbs[database]["name_long"],
                    "version": stg_dbs[database]["version"],
                    "signatures": 0,
                    "integrated_signatures": 0,
                    "recently_integrated": [],
                    "is_new": database in new_databases,
                    "is_updated": database in database_updates,
                    "sets": database_sets.get(database, 0)
                }

            obj["signatures"] += 1
            if entry["integrated"]:
                obj["integrated_signatures"] += 1

                if acc not in rel_integrated:
                    # Recent integration
                    obj["recently_integrated"].append(acc)

    notes.update({
        "interpro": {
            "entries": sum(interpro_types.values()),
            "new_entries": new_entries,
            "latest_entry": latest_entry,
            "types": interpro_types,
            "go_terms": n_interpro2go
        },
        "member_databases": member_databases,
        "citations": len(citations)
    })

    con = MySQLdb.connect(**parse_url(stg_url), charset="utf8")
    cur = con.cursor()
    cur.execute(
        """
        SELECT COUNT(*)
        FROM webfront_release_note
        WHERE version = %s
        """,
        (version,)
    )
    n_rows, = cur.fetchone()

    if n_rows:
        cur.execute(
            """
            UPDATE webfront_release_note
            SET content = %s
            WHERE version = %s
            """,
            (json.dumps(notes), version)
        )
    else:
        cur.execute(
            """
            INSERT INTO webfront_release_note
            VALUES (%s, %s, %s)
            """,
            (
                version,
                datetime.strptime(release_date, "%Y-%m-%d"),
                json.dumps(notes)
            )
        )

    cur.close()
    con.commit()
    con.close()
    logger.info("complete")
