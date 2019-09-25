# -*- coding: utf-8 -*-

import json
from datetime import datetime

import MySQLdb

from i7dw import io, logger
from i7dw.interpro import mysql


def make_release_notes(stg_url: str, rel_url: str, src_proteins: str,
                       src_matches: str, src_proteomes: str, version: str,
                       release_date: str):
    con = MySQLdb.connect(**mysql.parse_url(stg_url), use_unicode=True,
                          charset="utf8")
    cur = con.cursor()

    # Get PDB structures
    cur.execute(
        """
        SELECT accession, release_date
        FROM webfront_structure
        ORDER BY release_date
        """
    )
    structures = set()
    pdbe_release_date = None
    for row in cur:
        structures.add(row[0])
        pdbe_release_date = row[1]

    # Get proteomes
    proteomes = set(mysql.proteome.get_proteomes(stg_url))

    # Get taxa
    taxa = set(mysql.taxonomy.get_taxa(stg_url, lineage=False))

    # Integrated signatures
    integrated = {
        acc
        for acc, e in mysql.entry.get_entries(stg_url).items()
        if e["integrated"]
    }

    # Protein to PDBe structures
    protein2pdb = {}
    for pdb_id, s in mysql.structure.get_structures(stg_url).items():
        for acc in s["proteins"]:
            if acc in protein2pdb:
                protein2pdb[acc].add(pdb_id)
            else:
                protein2pdb[acc] = {pdb_id}

    # Get UniProtKB version
    cur.execute(
        """
        SELECT name_long, version
        FROM webfront_database
        WHERE type='protein'
        """
    )

    uniprot = {
        name: {
            "version": version,
            "count": 0,
            "signatures": 0,
            "integrated_signatures": 0
        }
        for name, version in cur
    }

    # Get sets
    db2set = {}
    for set_ac, s in mysql.entry.get_sets(stg_url).items():
        db = s["database"]
        if db in db2set:
            db2set[db] += 1
        else:
            db2set[db] = 1

    cur.close()
    con.close()

    proteins = io.Store(src_proteins)
    protein2matches = io.Store(src_matches)
    protein2proteome = io.Store(src_proteomes)

    interpro_structures = set()
    interpro_proteomes = set()
    interpro_taxa = set()

    n_proteins = 0
    for acc, protein in proteins:
        if protein["is_reviewed"]:
            k = "UniProtKB/Swiss-Prot"
        else:
            k = "UniProtKB/TrEMBL"

        uniprot[k]["count"] += 1
        matches = protein2matches.get(acc)
        if matches:
            # Protein has a least one signature
            uniprot[k]["signatures"] += 1

            # Search if the protein has at least one integrated signature
            for m in matches:
                if m["method_ac"] in integrated:
                    # It has!
                    uniprot[k]["integrated_signatures"] += 1

                    # Add taxon, proteome, and structures
                    interpro_taxa.add(protein["taxon"])

                    upid = protein2proteome.get(acc)
                    if upid:
                        interpro_proteomes.add(upid)

                    interpro_structures |= protein2pdb.get(acc, set())
                    break

        n_proteins += 1
        if not n_proteins % 10000000:
            logger.info("{:>12,}".format(n_proteins))

    proteins.close()
    protein2matches.close()
    protein2proteome.close()

    for k in ("count", "signatures", "integrated_signatures"):
        uniprot["UniProtKB"][k] = (uniprot["UniProtKB/Swiss-Prot"][k]
                                   + uniprot["UniProtKB/TrEMBL"][k])

    logger.info("{:>12,}".format(n_proteins))

    bad = interpro_structures - structures
    if bad:
        logger.warning("structures issues: {}".format(bad))

    bad = interpro_proteomes - proteomes
    if bad:
        logger.warning("proteomes issues: {}".format(bad))

    bad = interpro_taxa - taxa
    if bad:
        logger.warning("taxonomy issues: {}".format(bad))

    notes = {
        "interpro": {},
        "member_databases": [],
        "proteins": uniprot,
        "structures": {
            "total": len(structures),
            "integrated": len(interpro_structures & structures),
            "version": pdbe_release_date.strftime("%Y-%m-%d")
        },
        "proteomes": {
            "total": len(proteomes),
            "integrated": len(interpro_proteomes & proteomes),
            "version": uniprot["UniProtKB"]["version"]
        },
        "taxonomy": {
            "total": len(taxa),
            "integrated": len(interpro_taxa & taxa),
            "version": uniprot["UniProtKB"]["version"]
        }
    }

    rel_interpro_entries = set()
    already_integrated = set()
    for e in mysql.entry.get_entries(rel_url).values():
        if e["database"] == "interpro":
            rel_interpro_entries.add(e["accession"])
        elif e["integrated"]:
            # Signature already integrated durlng the previous release
            already_integrated.add(e["accession"])

    stg_entries = []
    new_entries = []
    for acc, e in mysql.entry.get_entries(stg_url).items():
        stg_entries.append(e)
        if e["database"] == "interpro" and acc not in rel_interpro_entries:
            new_entries.append(acc)

    # Member database changes
    stg_dbs = mysql.database.get_databases(stg_url)
    rel_dbs = mysql.database.get_databases(rel_url)
    updated_databases = set()
    new_databases = set()
    for name, info in stg_dbs.items():
        if name not in rel_dbs:
            new_databases.add(name)
        elif info["version"] != rel_dbs[name]["version"]:
            updated_databases.add(name)

    member_databases = {}
    interpro_types = {}
    citations = set()
    n_interpro2go = 0
    latest_entry = None
    for entry in sorted(stg_entries, key=lambda x: x["accession"]):
        acc = entry["accession"]
        db_name = entry["database"]
        _type = entry["type"]

        citations |= {
            item["PMID"]
            for item in entry["citations"].values()
            if item["PMID"] is not None
        }

        if db_name == "interpro":
            if _type in interpro_types:
                interpro_types[_type] += 1
            else:
                interpro_types[_type] = 1

            n_interpro2go += len(entry["go_terms"])

            latest_entry = acc
        else:
            if db_name in member_databases:
                db = member_databases[db_name]
            else:
                db = member_databases[db_name] = {
                    "name": stg_dbs[db_name]["name_long"],
                    "version": stg_dbs[db_name]["version"],
                    "signatures": 0,
                    "integrated_signatures": 0,
                    "recently_integrated": [],
                    "is_new": db_name in new_databases,
                    "is_updated": db_name in updated_databases,
                    "sets": db2set.get(db_name, 0)
                }

            db["signatures"] += 1
            if entry["integrated"]:
                db["integrated_signatures"] += 1

                if acc not in already_integrated:
                    # Recent integration
                    db["recently_integrated"].append(acc)

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

    con = MySQLdb.connect(**mysql.parse_url(stg_url), use_unicode=True,
                          charset="utf8")
    cur = con.cursor()
    cur.execute(
        """
        SELECT COUNT(*)
        FROM webfront_release_note
        WHERE version = %s
        """,
        (version,)
    )
    n = cur.fetchone()[0]

    if n:
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
            INSERT INTO webfront_release_note (
              version, release_date, content
            )
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
