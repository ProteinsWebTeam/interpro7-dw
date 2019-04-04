import json
import time
from datetime import datetime

from .database import get_databases
from .entry import get_entries, get_sets
from .proteome import get_proteomes
from .structure import get_structures
from .taxonomy import get_taxa
from ... import dbms, logger
from ...io import Store


def make_release_notes(stg_uri: str, rel_uri: str, src_proteins: str,
                       src_matches: str, src_proteomes: str, version: str,
                       release_date: str):
    con, cur = dbms.connect(stg_uri)

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
    proteomes = set(get_proteomes(stg_uri))

    # Get taxa
    taxa = set(get_taxa(stg_uri, lineage=False))

    # Integrated signatures
    integrated = {
        acc
        for acc, e in get_entries(stg_uri).items()
        if e["integrated"]
    }

    # Protein to PDBe structures
    protein2pdb = {}
    for pdb_id, s in get_structures(stg_uri).items():
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
    for set_ac, s in get_sets(stg_uri).items():
        db = s["database"]
        if db in db2set:
            db2set[db] += 1
        else:
            db2set[db] = 1

    cur.close()
    con.close()

    proteins = Store(src_proteins)
    protein2matches = Store(src_matches)
    protein2proteome = Store(src_proteomes)

    interpro_structures = set()
    interpro_proteomes = set()
    interpro_taxa = set()

    n_proteins = 0
    ts = time.time()
    for acc, protein in proteins:
        if protein["isReviewed"]:
            k = "UniProtKB/Swiss-Prot"
        else:
            k = "UniProtKB/TrEMBL"

        uniprot[k]["count"] += 1
        matches = protein2matches.get(acc)
        if matches:
            # Protein has a list one signature
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
            logger.info("{:>12,} ({:.0f} proteins/sec)".format(
                n_proteins, n_proteins / (time.time() - ts)
            ))

    proteins.close()
    protein2matches.close()
    protein2proteome.close()

    for k in ("count", "signatures", "integrated_signatures"):
        uniprot["UniProtKB"][k] = (uniprot["UniProtKB/Swiss-Prot"][k]
                                   + uniprot["UniProtKB/TrEMBL"][k])

    logger.info("{:>12,} ({:.0f} proteins/sec)".format(
        n_proteins, n_proteins / (time.time() - ts)
    ))

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
    for e in get_entries(rel_uri).values():
        if e["database"] == "interpro":
            rel_interpro_entries.add(e["accession"])
        elif e["integrated"]:
            # Signature already integrated during the previous release
            already_integrated.add(e["accession"])

    stg_entries = []
    new_entries = []
    for acc, e in get_entries(stg_uri).items():
        stg_entries.append(e)
        if e["database"] == "interpro" and acc not in rel_interpro_entries:
            new_entries.append(acc)

    # Member database changes
    stg_dbs = get_databases(stg_uri)
    rel_dbs = get_databases(rel_uri)
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

    con, cur = dbms.connect(stg_uri)
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
