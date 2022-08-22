import pickle

import MySQLdb

from interpro7dw import pfam
from interpro7dw.utils import logger
from interpro7dw.utils.mysql import uri2dict
from interpro7dw.interpro.oracle.entries import Entry
from interpro7dw.utils.store import BasicStore
from .utils import create_index, jsonify


def populate_annotations(uri: str, entries_file: str, hmms_file: str,
                         pfam_alignments: str):
    logger.info("loading entries")
    pfam2interpro = {}
    with open(entries_file, "rb") as fh:
        for e in pickle.load(fh).values():
            if e.database.lower() == "pfam" and e.integrated_in is not None:
                pfam2interpro[e.accession] = e.integrated_in

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

                # Pfam alignments: add for InterPro entry
                if (anno_type.startswith("alignment:")
                        and accession in pfam2interpro):
                    accession2 = pfam2interpro[accession]
                    cur.execute(
                        """
                        INSERT INTO webfront_entryannotation (
                            accession, type, value, mime_type, num_sequences
                        ) VALUES (%s, %s, %s, %s, %s)
                        """,
                        (accession2, anno_type, anno_value, mime_type, count)
                    )

    con.commit()
    cur.close()
    con.close()

    logger.info("done")


def index_annotations(uri: str):
    con = MySQLdb.connect(**uri2dict(uri), charset="utf8mb4")
    cur = con.cursor()
    create_index(
        cur,
        """
        CREATE INDEX i_entryannotation 
        ON webfront_entryannotation (accession)
        """
    )
    cur.close()
    con.close()



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
        if entry.database.lower() != "interpro":
            continue

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
    to_change, pfam2wiki = pfam.get_wiki(pfam_uri)
    for entry_acc, old_pages, new_pages in to_change:
        logger.warning(f"{entry_acc}: update following Wikipedia links:")
        for title in old_pages:
            logger.warning(f"\t- Remove: {title}")

        for title in new_pages:
            logger.warning(f"\t- Create: {title}")

    logger.info("loading Pfam curation/family details")
    pfam_details = pfam.get_details(pfam_uri)

    logger.info("loading clan members")
    entries_in_clan = set()
    with open(clans_file, "rb") as fh:
        for clan in pickle.load(fh).values():
            for entry_acc, _, _, _, _ in clan["members"]:
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

    # InterPro accession > member database > sig. accession > sig. name
    integrates = {}
    for entry_acc, entry in entries.items():
        if entry.integrated_in:
            try:
                mem_dbs = integrates[entry.integrated_in]
            except KeyError:
                mem_dbs = integrates[entry.integrated_in] = {}

            try:
                members = mem_dbs[entry.database.lower()]
            except KeyError:
                members = mem_dbs[entry.database.lower()] = {}

            members[entry_acc] = entry.name or entry.short_name or entry_acc

    logger.info("creating table")
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
          %s, %s, %s, %s, %s, %s, %s, %s, %s)
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

            r"""
            In several places, we convert the database names to lower case.
            This is because the API/client relies on lower cases. ¯\_(ツ)_/¯
            """
            if entry.old_names or entry.old_integrations:
                for key in list(entry.old_integrations.keys()):
                    value = entry.old_integrations.pop(key)
                    entry.old_integrations[key.lower()] = value

                history = {
                    "names": entry.old_names,
                    "signatures": entry.old_integrations
                }
            else:
                history = {}

            for key in list(entry.cross_references.keys()):
                value = entry.cross_references.pop(key)
                entry.cross_references[key.lower()] = value

            record = (
                None,
                entry.accession,
                entry.type.lower(),
                entry.name,
                entry.short_name,
                entry.database.lower(),
                jsonify(integrates.get(entry.accession), nullable=True),
                entry.integrated_in,
                jsonify(entry.go_terms, nullable=True),
                jsonify(entry.descriptions, nullable=True),
                # TODO: add support for multiple Wikipedia articles
                jsonify(pfam2wiki.get(entry.accession, [None])[0],
                        nullable=True),
                jsonify(pfam_details.get(entry.accession), nullable=True),
                jsonify(entry.literature, nullable=True),
                jsonify(hierarchy.get(entry.accession), nullable=True),
                jsonify(entry.cross_references, nullable=True),
                jsonify(entry.ppi, nullable=True),
                jsonify(pathways, nullable=True),
                jsonify(overlaps_with.get(entry.accession, []), nullable=True),
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
                    "sets": 1 if entry.accession in entries_in_clan else 0,
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
            None,
            entry.accession,
            entry.type.lower(),
            entry.name,
            entry.short_name,
            entry.database.lower(),
            jsonify(integrates.get(entry.accession), nullable=True),
            entry.integrated_in,
            jsonify(entry.go_terms, nullable=True),
            jsonify(entry.descriptions, nullable=True),
            # TODO: add support for multiple Wikipedia articles
            jsonify(pfam2wiki.get(entry.accession, [None])[0],
                    nullable=True),
            jsonify(pfam_details.get(entry.accession), nullable=True),
            jsonify(entry.literature, nullable=True),
            jsonify(hierarchy.get(entry.accession), nullable=True),
            jsonify(entry.cross_references, nullable=True),
            jsonify(entry.ppi, nullable=True),
            jsonify(pathways, nullable=True),
            jsonify(overlaps_with.get(entry.accession, []), nullable=True),
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
                "sets": 1 if entry.accession in entries_in_clan else 0,
                "structural_models": {
                    "alphafold": 0,
                    "rosettafold": 0
                },
                "structures": 0,
                "taxa": 0,
            }, nullable=False)
        )

        cur.execute(query, record)

    con.commit()
    cur.close()
    con.close()

    logger.info("done")


def index_entries(uri: str):
    con = MySQLdb.connect(**uri2dict(uri), charset="utf8mb4")
    cur = con.cursor()
    create_index(
        cur,
        """
        CREATE INDEX i_entry_database
        ON webfront_entry (source_database)
        """
    )
    create_index(
        cur,
        """
        CREATE INDEX i_entry_integrated
        ON webfront_entry (integrated_id)
        """
    )
    cur.close()
    con.close()


def populate_entry_taxa_distrib(uri: str, entries_file: str, xrefs_file: str):
    logger.info("loading entries")
    with open(entries_file, "rb") as fh:
        entries: dict[str, Entry] = pickle.load(fh)

    logger.info("creating table")
    con = MySQLdb.connect(**uri2dict(uri), charset="utf8mb4")
    cur = con.cursor()
    cur.execute("DROP TABLE IF EXISTS webfront_entrytaxa")
    cur.execute(
        """
        CREATE TABLE webfront_entrytaxa
        (
            accession VARCHAR(25) PRIMARY KEY NOT NULL,
            tree LONGTEXT
        ) CHARSET=utf8mb4 DEFAULT COLLATE=utf8mb4_unicode_ci
        """
    )

    logger.info("populating table")
    query = "INSERT INTO webfront_entrytaxa VALUES (%s, %s)"
    with BasicStore(xrefs_file, mode="r") as store:
        for accession, xrefs in store:
            tree = xrefs["taxa"]["tree"]
            cur.execute(query, (accession, jsonify(tree, nullable=True)))
            entries.pop(accession)

    for accession in entries:
        cur.execute(query, (accession, None))

    con.commit()
    cur.close()
    con.close()

    logger.info("done")
