import pickle

import MySQLdb

from interpro7dw.utils import logger
from interpro7dw.utils.mysql import uri2dict
from interpro7dw.utils.store import KVStore
from .utils import jsonify


def populate_databases(uri: str, databases_file: str):
    logger.info("creating webfront_database")
    with open(databases_file, "rb") as fh:
        databases = pickle.load(fh)

    params = []
    for key, info in databases.items():
        params.append((
            key.lower(),
            key,
            info["name"],
            info["description"],
            info["type"],
            info["entries"],
            info["release"]["version"],
            info["release"]["date"],
            info["previous_release"]["version"],
            info["release"]["date"]
        ))

    con = MySQLdb.connect(**uri2dict(uri), charset="utf8mb4")
    cur = con.cursor()
    cur.execute("DROP TABLE IF EXISTS webfront_database")
    cur.execute(
        """
        CREATE TABLE webfront_database
        (
            name VARCHAR(10) NOT NULL PRIMARY KEY,
            name_alt VARCHAR(10) NOT NULL,
            name_long VARCHAR(25) NOT NULL,
            description LONGTEXT,
            type ENUM('protein', 'entry', 'feature', 'other') NOT NULL,
            num_entries INTEGER,
            version VARCHAR(20),
            release_date DATETIME,
            prev_version VARCHAR(20),
            prev_release_date DATETIME
        ) CHARSET=utf8mb4 DEFAULT COLLATE=utf8mb4_unicode_ci
        """
    )

    cur.executemany(
        """
        INSERT INTO webfront_database 
        VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
        """,
        params
    )
    con.commit()
    cur.close()
    con.close()

    logger.info("done")


def populate_rel_notes(stg_uri: str, rel_uri: str, clans_file: str,
                       entries_file: str, proteomeinfo_file: str,
                       structures_file: str, structureinfo_file: str,
                       taxa_file: str, proteins_file: str, matches_file: str,
                       proteomes_file: str, **kwargs):
    notes = kwargs.get("notes", [])

    logger.info("loading clans")
    member2clan = {}
    with open(clans_file, "rb") as fh:
        for clan_acc, clan in pickle.load(fh).items():
            for entry_acc, _, _, _, _ in clan["members"]:
                member2clan[entry_acc] = clan_acc

    logger.info("loading entries")
    with open(entries_file, "rb") as fh:
        entries = pickle.load(fh)

    logger.info("loading PDBe structures")
    with open(structures_file, "rb") as fh:
        protein2structures = pickle.load(fh)

    with open(structureinfo_file, "rb") as fh:
        structures = list(pickle.load(fh)["entries"].values())

    logger.info("loading sequence databases")
    seq_databases = {}
    con = MySQLdb.connect(**uri2dict(stg_uri))
    cur = con.cursor()
    cur.execute(
        """
        SELECT name, name_long, version
        FROM webfront_database
        WHERE name in ('reviewed', 'unreviewed', 'uniprot')
        """
    )

    for name, name_long, version in cur:
        seq_databases[name] = {
            "name": name_long,
            "version": version,
            "count": 0,         # total number of proteins
            "hit": 0,           # number of proteins with at least one hit
            "integrated": 0     # number of integrated proteins
        }

    cur.close()
    con.close()

    # Number of proteins with GO terms from InterPro
    uniprot2go = 0

    # Entities found in InterPro
    integrated_proteomes = set()
    integrated_structures = set()
    integrated_taxa = set()

    proteins_store = KVStore(proteins_file)
    matches_store = KVStore(matches_file)
    proteomes_store = KVStore(proteomes_file)
    i = 0
    for i, (protein_acc, protein) in enumerate(proteins_store.items()):
        if (i + 1) % 1e7 == 0:
            logger.info(f"{i + 1:>15,}")

        if protein["reviewed"]:
            database = seq_databases["reviewed"]
        else:
            database = seq_databases["unreviewed"]

        database["count"] += 1

        try:
            signature_matches, entry_matches = matches_store[protein_acc]
        except KeyError:
            continue  # No matches

        # Protein matched by at least one signature
        database["hit"] += 1

        is_integrated = False
        for entry_acc in entry_matches:
            """
            Protein matched by at least one InterPro entry,
            i.e. at least one integrated signature
            """
            is_integrated = True

            if entries[entry_acc].go_terms:
                uniprot2go += 1
                break

        if is_integrated:
            database["integrated"] += 1

            proteome_id = proteomes_store.get(protein_acc)
            if proteome_id:
                integrated_proteomes.add(proteome_id)

            protein_structures = protein2structures.get(protein_acc, {})
            if protein_structures:
                integrated_structures |= set(protein_structures.keys())

            integrated_taxa.add(protein["taxid"])

    logger.info(f"{i + 1:>15,}")

    proteomes_store.close()
    matches_store.close()
    proteomes_store.close()

    # Sums Swiss-Prot and TrEMBL counts
    for key in ("count", "hit", "integrated"):
        seq_databases["uniprot"][key] = (seq_databases["reviewed"][key]
                                         + seq_databases["unreviewed"][key])

    logger.info("tracking changes since last releases")
    con = MySQLdb.connect(**uri2dict(rel_uri), charset="utf8mb4")
    cur = con.cursor()
    cur.execute(
        """
        SELECT accession, source_database, integrated_id
        FROM webfront_entry
        WHERE is_alive = 1
        """
    )
    public_entries = set()
    public_integrated = set()
    for entry_acc, database, integrated_in in cur:
        if database == "interpro":
            public_entries.add(entry_acc)
        elif integrated_in:
            # Signature already integrated in the previous release
            public_integrated.add(entry_acc)

    cur.execute(
        """
        SELECT name, version
        FROM webfront_database
        WHERE type = 'entry'
        """
    )
    public_databases = dict(cur.fetchall())

    cur.execute("SELECT * FROM webfront_release_note")
    prev_releases = cur.fetchall()
    cur.close()
    con.close()

    con = MySQLdb.connect(**uri2dict(stg_uri), charset="utf8mb4")
    cur = con.cursor()
    cur.execute(
        """
        SELECT name, name_long, version, release_date
        FROM webfront_database
        WHERE type = 'entry'
        """
    )
    staging_databases = {row[0]: (row[1], row[2], row[3]) for row in cur}

    interpro_new = []
    interpro_types = {}
    member_databases = {}
    pubmed_citations = set()
    interpro2go = 0
    latest_entry = None

    for entry in sorted(entries.values(), key=lambda e: e.creation_date):
        if entry.deletion_date is not None:
            continue

        dbkey = entry.database.lower()
        if dbkey == "interpro":
            for pub in entry.literature.values():
                if pub["PMID"] is not None:
                    pubmed_citations.add(pub["PMID"])

            try:
                interpro_types[entry.type.lower()] += 1
            except KeyError:
                interpro_types[entry.type.lower()] = 1

            if entry.accession not in public_entries:
                interpro_new.append(entry.accession)

            interpro2go += len(entry.go_terms)
            latest_entry = entry.accession
        else:
            try:
                obj = member_databases[dbkey]
            except KeyError:
                database, version, _ = staging_databases[dbkey]

                is_new = is_updated = False
                if dbkey not in public_databases:
                    is_new = True
                elif version != public_databases[dbkey]:
                    is_updated = True

                obj = member_databases[dbkey] = {
                    "name": database,
                    "version": version,
                    "signatures": 0,
                    "integrated_signatures": 0,
                    "recently_integrated": [],
                    "is_new": is_new,
                    "is_updated": is_updated,
                    "sets": set()
                }

            obj["signatures"] += 1
            if entry.integrated_in:
                obj["integrated_signatures"] += 1

                if entry.accession not in public_integrated:
                    # Recent integration
                    obj["recently_integrated"].append(entry.accession)

            if entry.accession in member2clan:
                obj["sets"].add(member2clan[entry.accession])

    # Transforms sets of clans to counts:
    for obj in member_databases.values():
        obj["sets"] = len(obj["sets"])

    # Checks that "integrated" proteomes and taxa are valid
    with open(proteomeinfo_file, "rb") as fh:
        proteomes = set(pickle.load(fh).keys())

    errors = integrated_proteomes - proteomes
    if errors:
        raise RuntimeError(f"invalid proteomes: {', '.join(errors)}")

    with open(taxa_file, "rb") as fh:
        taxa = set(pickle.load(fh).keys())

    errors = integrated_taxa - taxa
    if errors:
        raise RuntimeError(f"invalid taxa: {', '.join(errors)}")

    logger.info("creating webfront_release_note")
    cur.execute("DROP TABLE IF EXISTS webfront_release_note")
    cur.execute(
        """
        CREATE TABLE webfront_release_note
        (
            version VARCHAR(20) PRIMARY KEY NOT NULL,
            release_date DATETIME NOT NULL,
            content LONGTEXT NOT NULL
        ) CHARSET=utf8mb4 DEFAULT COLLATE=utf8mb4_unicode_ci
        """
    )

    # Adds previous release notes
    cur.executemany(
        """
        INSERT INTO webfront_release_note
        VALUES (%s, %s, %s)
        """, prev_releases
    )
    con.commit()

    uniprot_version = seq_databases["uniprot"]["version"]

    # Rename keys (used by API/website)
    for key in list(seq_databases.keys()):
        value = seq_databases.pop(key)
        value["signatures"] = value.pop("hit")
        value["integrated_signatures"] = value.pop("integrated")

        new_key = value.pop("name")
        seq_databases[new_key] = value

    content = {
        "notes": notes,
        "interpro": {
            "entries": sum(interpro_types.values()),
            "new_entries": interpro_new,
            "latest_entry": latest_entry,
            "types": interpro_types,
            "go_terms": interpro2go,
            "uniprot2go": uniprot2go
        },
        "member_databases": member_databases,
        "proteins": seq_databases,
        "structures": {
            "total": len(structures),
            "integrated": len(integrated_structures),
            "version": max(entry["date"]
                           for entry in structures).strftime("%Y-%m-%d")
        },
        "proteomes": {
            "total": len(proteomes),
            "integrated": len(integrated_proteomes),
            "version": uniprot_version
        },
        "taxonomy": {
            "total": len(taxa),
            "integrated": len(integrated_taxa),
            "version": uniprot_version
        },
        "citations": len(pubmed_citations)
    }

    _, version, date = staging_databases["interpro"]
    cur.execute(
        """
        SELECT COUNT(*)
        FROM webfront_release_note
        WHERE version = %s
        """, (version,)
    )
    n_rows, = cur.fetchone()

    if n_rows:
        # Release notes already in the table: update it
        cur.execute(
            """
            UPDATE webfront_release_note
            SET content = %s
            WHERE version = %s
            """, (jsonify(content), version)
        )
    else:
        # Adds new release notes
        cur.execute(
            """
            INSERT INTO webfront_release_note
            VALUES (%s, %s, %s)
            """, (version, date, jsonify(content))
        )

    con.commit()
    cur.close()
    con.close()

    logger.info("done")
