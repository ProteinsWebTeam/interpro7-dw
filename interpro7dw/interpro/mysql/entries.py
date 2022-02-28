import pickle

import MySQLdb

from interpro7dw import pfam
from interpro7dw.utils import logger
from interpro7dw.utils.mysql import uri2dict
from interpro7dw.interpro.oracle.entries import Entry
from interpro7dw.utils.store import BasicStore, KVStore
from .utils import jsonify


def populate_annotations(uri: str, hmms_file: str, pfam_alignments: str):
    logger.info("creating webfront_entryannotation")
    con = MySQLdb.connect(**uri2dict(uri))
    cur = con.cursor()
    cur.execute("DROP TABLE IF EXISTS webfront_entryannotation")
    cur.execute(
        """
        CREATE TABLE webfront_entryannotation
        (
            annotation_id INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
            accession VARCHAR(25) NOT NULL,
            type VARCHAR(20) NOT NULL,
            value LONGBLOB NOT NULL,
            mime_type VARCHAR(32) NOT NULL,
            num_sequences INT
        ) CHARSET=utf8mb4 DEFAULT COLLATE=utf8mb4_unicode_ci
        """
    )

    for file in [hmms_file, pfam_alignments]:
        with BasicStore(file, mode="r") as store:
            for accession, anno_type, anno_value, count in store:
                if anno_type == "logo":
                    mime_type = "application/json"
                else:
                    mime_type = "application/gzip"

                cur.execute(
                    """
                    INSERT INTO webfront_entryannotation (
                        accession, type, value, mime_type, num_sequences
                    ) VALUES (%s, %s, %s, %s, %s)
                    """, (accession, anno_type, anno_value, mime_type, count)
                )

    con.commit()

    logger.info("indexing")
    cur.execute(
        """
        CREATE INDEX i_entryannotation 
        ON webfront_entryannotation (accession)
        """
    )

    cur.close()
    con.close()

    logger.info("done")


def make_hierarchy(entries: dict[str, Entry]) -> dict:
    child2parent = {}
    parent2children = {}

    for entry in entries.values():
        if entry.parent:
            child2parent[entry.accession] = entry.parent
            if entry.parent in parent2children:
                parent2children[entry.parent].append(entry.accession)
            else:
                parent2children[entry.parent] = [entry.accession]

    hierarchy = {}
    for entry in entries.values():
        # Find root
        accession = entry.accession
        parent_acc = child2parent.get(accession)

        while parent_acc:
            accession = parent_acc
            parent_acc = child2parent.get(accession)

        hierarchy[entry.accession] = format_node(accession, entries,
                                                 parent2children)

    return hierarchy


def format_node(accession: str, entries: dict[str, Entry],
                parent2children: dict[str, list[str]]) -> dict:
    children = []
    for child_acc in sorted(parent2children.get(accession, [])):
        children.append(format_node(child_acc, entries, parent2children))

    entry = entries[accession]
    return {
        "accession": accession,
        "name": entry.name,
        "type": entry.type,
        "children": children
    }


def populate_entries(ipr_uri: str, pfam_uri: str, clans_file: str,
                     entries_file: str, overlapping_file: str,
                     xrefs_file: str):
    logger.info("fetching Wikipedia data for Pfam entries")
    wiki = pfam.get_wiki(pfam_uri)

    logger.info("loading Pfam curation/family details")
    pfam_details = pfam.get_details(pfam_uri)

    logger.info("loading clan members")
    entries_in_clan = set()
    with open(clans_file, "rb") as fh:
        for clan in pickle.load(fh).values():
            for entry_acc, _, _ in clan["members"]:
                entries_in_clan.add(entry_acc)

    logger.info("loading entries")
    with open(entries_file, "rb") as fh:
        entries: dict[str, Entry] = pickle.load(fh)

    overlaps_with = {}
    with open(overlapping_file, "rb") as fh:
        for acc_1, acc_2 in pickle.load(fh):
            for query, target in [(acc_1, acc_2), (acc_2, acc_1)]:
                try:
                    others = overlaps_with[query]
                except KeyError:
                    others = overlaps_with[query] = []

                other = entries[target]
                others.append({
                    "accession": other.accession,
                    "name": other.name,
                    "type": other.type.lower()
                })

    hierarchy = make_hierarchy(entries)
    integrates = {}
    for entry_acc, entry in entries.items():
        if entry.integrated_in:
            key = entry.integrated_in
            try:
                mem_db = integrates[key][entry.database.lower()]
            except KeyError:
                mem_db = integrates[key][entry.database.lower()] = {}

            mem_db[entry_acc] = entry.name or entry.short_name or entry_acc


    con = MySQLdb.connect(**uri2dict(ipr_uri), charset="utf8mb4")
    cur = con.cursor()
    cur.execute("DROP TABLE IF EXISTS webfront_entry")
    cur.execute(
        """
        CREATE TABLE webfront_entry
        (
            entry_id VARCHAR(10) DEFAULT NULL,
            accession VARCHAR(25) PRIMARY KEY NOT NULL,
            type VARCHAR(50) NOT NULL,
            name LONGTEXT,
            short_name VARCHAR(100),
            source_database VARCHAR(10) NOT NULL,
            member_databases LONGTEXT,
            integrated_id VARCHAR(25),
            go_terms LONGTEXT,
            description LONGTEXT,
            wikipedia LONGTEXT,
            details LONGTEXT,
            literature LONGTEXT,
            hierarchy LONGTEXT,
            cross_references LONGTEXT,
            interactions LONGTEXT,
            pathways LONGTEXT,
            overlaps_with LONGTEXT,
            taxa LONGTEXT,
            is_featured TINYINT NOT NULL,
            is_alive TINYINT NOT NULL,
            history LONGTEXT,
            entry_date DATETIME NOT NULL,
            deletion_date DATETIME,
            counts LONGTEXT NOT NULL
        ) CHARSET=utf8mb4 DEFAULT COLLATE=utf8mb4_unicode_ci
        """
    )

    query = """
        INSERT INTO webfront_entry
        VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s,
          %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
    """

    with BasicStore(xrefs_file, mode="r") as store:
        for entry_acc, xrefs in store:
            entry = entries.pop(entry_acc)

            if xrefs["enzymes"]:
                entry.cross_references["ec"] = sorted(xrefs["enzymes"])

            pathways = {}
            for key in ["metacyc", "reactome"]:
                if xrefs[key]:
                    pathways[key] = [dict(zip(("id", "name"), item))
                                     for item in xrefs[key]]

            if entry.old_names or entry.old_integrations:
                history = {
                    "names": entry.old_names,
                    "signatures": entry.old_integrations
                }
            else:
                history = {}

            record = (
                entry_acc,
                entry.type.lower(),
                entry.name,
                entry.short_name,
                entry.database.lower(),
                jsonify(integrates.get(entry_acc), nullable=True),
                entry.integrated_in,
                jsonify(entry.go_terms, nullable=True),
                jsonify(entry.descriptions, nullable=True),
                jsonify(wiki.get(entry_acc), nullable=True),
                jsonify(pfam_details.get(entry.accession), nullable=True),
                jsonify(entry.literature, nullable=True),
                jsonify(hierarchy.get(entry_acc), nullable=True),
                jsonify(entry.cross_references, nullable=True),
                jsonify(entry.ppi, nullable=True),
                jsonify(pathways, nullable=True),
                jsonify(overlaps_with.get(entry_acc, []), nullable=True),
                jsonify(xrefs["taxa"]["tree"], nullable=True),
                0,
                1 if entry.deletion_date is None else 0,
                jsonify(history, nullable=True),
                entry.creation_date,
                entry.deletion_date,
                jsonify({
                    "domain_architectures": len(xrefs["dom_orgs"]),
                    "interactions": len(entry.ppi),
                    "matches": xrefs["matches"],
                    "pathways": sum([len(v) for v in pathways.values()]),
                    "proteins": len(xrefs["proteins"]),
                    "proteomes": len(xrefs["proteomes"]),
                    "sets": 1 if entry_acc in entries_in_clan else 0,
                    "structural_models": xrefs["struct_models"],
                    "structures": len(xrefs["structures"]),
                    "taxa": len(xrefs["taxa"]["all"]),
                }, nullable=False)
            )

            cur.execute(query, record)

    # Add entries without cross-references
    for entry in entries.values():
        if entry.old_names or entry.old_integrations:
            history = {
                "names": entry.old_names,
                "signatures": entry.old_integrations
            }
        else:
            history = {}

        record = (
            entry_acc,
            entry.type.lower(),
            entry.name,
            entry.short_name,
            entry.database.lower(),
            jsonify(integrates.get(entry_acc), nullable=True),
            entry.integrated_in,
            jsonify(entry.go_terms, nullable=True),
            jsonify(entry.descriptions, nullable=True),
            jsonify(wiki.get(entry_acc), nullable=True),
            jsonify(pfam_details.get(entry.accession), nullable=True),
            jsonify(entry.literature, nullable=True),
            jsonify(hierarchy.get(entry_acc), nullable=True),
            jsonify(entry.cross_references, nullable=True),
            jsonify(entry.ppi, nullable=True),
            jsonify(pathways, nullable=True),
            jsonify(overlaps_with.get(entry_acc, []), nullable=True),
            jsonify(xrefs["taxa"]["tree"], nullable=True),
            0,
            1 if entry.deletion_date is None else 0,
            jsonify(history, nullable=True),
            entry.creation_date,
            entry.deletion_date,
            jsonify({
                "domain_architectures": 0,
                "interactions": len(entry.ppi),
                "matches": 0,
                "pathways": 0,
                "proteins": 0,
                "proteomes": 0,
                "sets": 1 if entry_acc in entries_in_clan else 0,
                "structural_models": 0,
                "structures": 0,
                "taxa": 0,
            }, nullable=False)
        )

        cur.execute(query, record)

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
            for entry_acc, _, _ in clan["members"]:
                member2clan[entry_acc] = clan_acc

    logger.info("loading entries")
    with open(entries_file, "rb") as fh:
        entries = pickle.load(fh)

    logger.info("loading PDBe structures")
    with open(structures_file, "rb") as fh:
        protein2structures = pickle.load(fh)

    with open(structureinfo_file, "rb") as fh:
        structures = pickle.load(fh)["entries"]

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
            "total": 0,         # total number of proteins
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

        database["total"] += 1

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
    for key in ("total", "hit", "integrated"):
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

        if entry.database.lower() == "interpro":
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
                obj = member_databases[entry.database]
            except KeyError:
                database, version, _ = staging_databases[entry.database]

                is_new = is_updated = False
                if entry.database not in public_databases:
                    is_new = True
                elif version != public_databases[entry.database]:
                    is_updated = True

                obj = member_databases[entry.database] = {
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
            "version": seq_databases["uniprot"]["version"]
        },
        "taxonomy": {
            "total": len(taxa),
            "integrated": len(integrated_taxa),
            "version": seq_databases["uniprot"]["version"]
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
