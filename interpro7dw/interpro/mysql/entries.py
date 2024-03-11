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
            accession VARCHAR(30) NOT NULL,
            type VARCHAR(20) NOT NULL,
            value LONGBLOB NOT NULL,
            mime_type VARCHAR(32) NOT NULL,
            num_sequences INT
        ) CHARSET=utf8mb4 DEFAULT COLLATE=utf8mb4_unicode_ci
        """
    )

    ignore = {
        "alignment:rp15",
        "alignment:rp35",
        "alignment:rp55",
        "alignment:rp75",
        "alignment:full",
        "alignment:uniprot"
    }

    for file in [hmms_file, pfam_alignments]:
        with BasicStore(file, mode="r") as store:
            for accession, anno_type, anno_value, count in store:
                if anno_type in ignore:
                    continue
                elif anno_type == "logo":
                    mime_type = "application/json"
                else:
                    mime_type = "application/gzip"

                cur.execute(
                    """
                    INSERT INTO webfront_entryannotation (
                        accession, type, value, mime_type, num_sequences
                    ) VALUES (%s, %s, %s, %s, %s)
                    """,
                    [accession, anno_type, anno_value, mime_type, count]
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
                        [accession2, anno_type, anno_value, mime_type, count]
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
        if entry.deletion_date or not entry.public:
            continue

        # Find root
        accession = entry.accession
        parent_acc = child2parent.get(accession)

        while parent_acc:
            accession = parent_acc
            parent_acc = child2parent.get(accession)

        # Make hierarchy from root to entry
        hierarchy[entry.accession] = format_node(accession, entries,
                                                 parent2children)

    return hierarchy


def get_hierarchy(entry: Entry, hierarchy: dict[str, dict]) -> tuple:
    if entry.accession in hierarchy:
        entry_hierarchy = hierarchy[entry.accession]

        if entry.database.lower() == "interpro":
            # InterPro entry -> entry hierarchy
            return entry_hierarchy, 0
        elif entry.database.lower() in ("cathgene3d", "panther"):
            # Panther subfamilies, or CATH -> FunFams
            return None, len(entry_hierarchy["children"])

    return None, 0


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
                     xrefs_file: str, structures_file: str):
    logger.info("fetching Wikipedia data for Pfam entries")
    to_change, pfam2wiki = pfam.get_wiki(pfam_uri)
    # for entry_acc, old_pages, new_pages in to_change:
    #     logger.warning(f"{entry_acc}: update following Wikipedia links:")
    #     for title in old_pages:
    #         logger.warning(f"\t- Remove: {title}")
    #
    #     for title in new_pages:
    #         logger.warning(f"\t- Create: {title}")

    logger.info("loading Pfam curation/family details")
    pfam_details = pfam.get_details(pfam_uri)

    logger.info("loading clan members")
    entries_in_clan = {}
    with open(clans_file, "rb") as fh:
        for clan in pickle.load(fh).values():
            for entry_acc, _, _, _, _ in clan["members"]:
                entries_in_clan[entry_acc] = {
                    "accession": clan["accession"],
                    "name": clan["name"]
                }

    logger.info("loading structures")
    highres_structures = {}
    with open(structures_file, "rb") as fh:
        for s in pickle.load(fh).values():
            if s["resolution"] is not None and s["resolution"] <= 2:
                highres_structures[s["id"]] = s["name"]

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
        if entry.integrated_in is None:
            continue

        parent = entries[entry.integrated_in]
        if parent.database.lower() != "interpro":
            # Ignore PANTHER, FunFam hierarchies
            continue

        try:
            mem_dbs = integrates[parent.accession]
        except KeyError:
            mem_dbs = integrates[parent.accession] = {}

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
            accession VARCHAR(30) PRIMARY KEY NOT NULL,
            type VARCHAR(50) NOT NULL,
            name VARCHAR(400),
            short_name VARCHAR(100),
            source_database VARCHAR(10) NOT NULL,
            member_databases LONGTEXT,
            integrated_id VARCHAR(30),
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
            is_llm TINYINT NOT NULL,
            is_reviewed_llm TINYINT NOT NULL,
            is_public TINYINT NOT NULL,
            history LONGTEXT,
            entry_date DATETIME NOT NULL,
            deletion_date DATETIME,
            set_info TEXT,
            representative_structure LONGTEXT,
            counts LONGTEXT NOT NULL
        ) CHARSET=utf8mb4 DEFAULT COLLATE=utf8mb4_unicode_ci
        """
    )

    query = """
        INSERT INTO webfront_entry
        VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s,
          %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
    """

    inserted_entries = set()
    records = []
    with BasicStore(xrefs_file, mode="r") as store:
        for entry_acc, xrefs in store:
            entry = entries.pop(entry_acc)
            inserted_entries.add(entry_acc.lower())

            if xrefs["enzymes"]:
                entry.cross_references["ec"] = sorted(xrefs["enzymes"])

            pathways = {}
            for key in ["metacyc", "reactome"]:
                if xrefs[key]:
                    pathways[key] = [dict(zip(("id", "name"), item))
                                     for item in xrefs[key]]

            history = {}
            if entry.old_names:
                history["names"] = entry.old_names

            if entry.old_short_names:
                history["short_names"] = entry.old_short_names
            
            if entry.old_integrations:
                # Convert DB name to lower cases (API/client relies on LC)
                history["signatures"] = {
                    k.lower(): v 
                    for k, v in entry.old_integrations.items()
                }

            # Force keys of cross-references to lower case
            for key in list(entry.cross_references.keys()):
                value = entry.cross_references.pop(key)
                entry.cross_references[key.lower()] = value

            best_coverage = 0
            best_structure = None
            for pdb_id, coverage in xrefs["structures"]:
                if pdb_id not in highres_structures or coverage < 0.5:
                    continue
                elif coverage > best_coverage:
                    best_coverage = coverage
                    best_structure = {
                        "accession": pdb_id,
                        "name": highres_structures[pdb_id]
                    }

            entry_hierarchy, num_subfamilies = get_hierarchy(entry, hierarchy)
            entry_clan = entries_in_clan.get(entry.accession)
            records.append((
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
                jsonify(entry_hierarchy, nullable=True),
                jsonify(entry.cross_references, nullable=True),
                jsonify(entry.ppi, nullable=True),
                jsonify(pathways, nullable=True),
                jsonify(overlaps_with.get(entry.accession, []), nullable=True),
                1 if entry.llm else 0,
                1 if entry.llm_reviewed else 0,
                1 if entry.public else 0,
                jsonify(history, nullable=True),
                entry.creation_date,
                entry.deletion_date,
                jsonify(entry_clan, nullable=True),
                jsonify(best_structure),
                jsonify({
                    "subfamilies": num_subfamilies,
                    "domain_architectures": len(xrefs["dom_orgs"]),
                    "interactions": len(entry.ppi),
                    "matches": xrefs["matches"],
                    "pathways": sum([len(v) for v in pathways.values()]),
                    "proteins": len(xrefs["proteins"]),
                    "proteomes": len(xrefs["proteomes"]),
                    "sets": 1 if entry_clan else 0,
                    "structural_models": xrefs["struct_models"],
                    "structures": len(xrefs["structures"]),
                    "taxa": len(xrefs["taxa"]["all"]),
                }, nullable=False)
            ))

            if len(records) == 1000:
                cur.executemany(query, records)
                records.clear()

    # Add entries without cross-references
    for entry in entries.values():
        if entry.accession.lower() in inserted_entries:
            continue

        history = {}
        if entry.old_names:
            history["names"] = entry.old_names

        if entry.old_short_names:
            history["short_names"] = entry.old_short_names
        
        if entry.old_integrations:
            # Convert DB name to lower cases (API/client relies on LC)
            history["signatures"] = {
                k.lower(): v 
                for k, v in entry.old_integrations.items()
            }

        entry_clan = entries_in_clan.get(entry.accession)
        entry_hierarchy, num_subfamilies = get_hierarchy(entry, hierarchy)
        records.append((
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
            jsonify(entry_hierarchy, nullable=True),
            jsonify(entry.cross_references, nullable=True),
            jsonify(entry.ppi, nullable=True),
            jsonify(pathways, nullable=True),
            jsonify(overlaps_with.get(entry.accession, []), nullable=True),
            1 if entry.llm else 0,
            1 if entry.llm_reviewed else 0,
            1 if entry.public else 0,
            jsonify(history, nullable=True),
            entry.creation_date,
            entry.deletion_date,
            jsonify(entry_clan, nullable=True),
            None,
            jsonify({
                "subfamilies": num_subfamilies,
                "domain_architectures": 0,
                "interactions": len(entry.ppi),
                "matches": 0,
                "pathways": 0,
                "proteins": 0,
                "proteomes": 0,
                "sets": 1 if entry_clan else 0,
                "structural_models": {
                    "alphafold": 0,
                },
                "structures": 0,
                "taxa": 0,
            }, nullable=False)
        ))

    for i in range(0, len(records), 1000):
        cur.executemany(query, records[i:i+1000])

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
    create_index(
        cur,
        """
        CREATE INDEX i_entry_name
        ON webfront_entry (name)
        """
    )
    create_index(
        cur,
        """
        CREATE INDEX i_entry_short_name
        ON webfront_entry (short_name)
        """
    )
    create_index(
        cur,
        """
        CREATE INDEX i_entry_deletion_date
        ON webfront_entry (deletion_date)
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
            accession VARCHAR(30) PRIMARY KEY NOT NULL,
            tree LONGTEXT
        ) CHARSET=utf8mb4 DEFAULT COLLATE=utf8mb4_unicode_ci
        """
    )

    logger.info("populating table")
    query = "INSERT INTO webfront_entrytaxa VALUES (%s, %s)"
    with BasicStore(xrefs_file, mode="r") as store:
        for accession, xrefs in store:
            entry = entries.pop(accession)

            if entry.deletion_date is None and entry.public:
                tree = xrefs["taxa"]["tree"]
                cur.execute(query, (accession, jsonify(tree, nullable=True)))
                con.commit()

    for entry in entries.values():
        if entry.deletion_date is None and entry.public:
            cur.execute(query, (entry.accession, None))

    con.commit()
    cur.close()
    con.close()

    logger.info("done")
