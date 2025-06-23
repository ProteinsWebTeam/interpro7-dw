import pickle

from interpro7dw import intact
from interpro7dw.utils import logger
from interpro7dw.interpro.oracle.entries import Entry
from interpro7dw.utils.store import BasicStore
from .utils import connect, jsonify


_REPR_STRUCT_MIN_COVERAGE = 0.5
_REPR_STRUCT_MAX_RESOLUTION = 2


def populate_annotations(
    uri: str, entries_file: str, hmms_file: str, pfam_alignments: str
):
    logger.info("loading entries")
    pfam2interpro = {}
    with open(entries_file, "rb") as fh:
        for e in pickle.load(fh).values():
            if e.database.lower() == "pfam" and e.integrated_in is not None:
                pfam2interpro[e.accession] = e.integrated_in

    logger.info("creating webfront_entryannotation")
    con = connect(uri)
    cur = con.cursor()
    cur.execute("DROP TABLE IF EXISTS webfront_entryannotation")
    cur.execute(
        """
        CREATE TABLE webfront_entryannotation
        (
            annotation_id SERIAL PRIMARY KEY,
            accession VARCHAR(30) NOT NULL,
            type VARCHAR(20) NOT NULL,
            value BYTEA NOT NULL,
            mime_type VARCHAR(32) NOT NULL,
            num_sequences INTEGER
        )
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
                    """,
                    [accession, anno_type, anno_value, mime_type, count],
                )

                # Pfam alignments: add for InterPro entry
                if anno_type.startswith("alignment:") and accession in pfam2interpro:
                    accession2 = pfam2interpro[accession]
                    cur.execute(
                        """
                        INSERT INTO webfront_entryannotation (
                            accession, type, value, mime_type, num_sequences
                        ) VALUES (%s, %s, %s, %s, %s)
                        """,
                        [accession2, anno_type, anno_value, mime_type, count],
                    )

                con.commit()

    cur.close()
    con.close()

    logger.info("done")


def index_annotations(uri: str):
    con = connect(uri)
    cur = con.cursor()
    cur.execute(
        """
        CREATE INDEX IF NOT EXISTS i_entryannotation 
        ON webfront_entryannotation (accession)
        """
    )
    con.commit()
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
        hierarchy[entry.accession] = format_node(accession, entries, parent2children)

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


def format_node(
    accession: str, entries: dict[str, Entry], parent2children: dict[str, list[str]]
) -> dict:
    children = []
    for child_acc in sorted(parent2children.get(accession, [])):
        children.append(format_node(child_acc, entries, parent2children))

    entry = entries[accession]
    return {
        "accession": accession,
        "name": entry.name,
        "type": entry.type,
        "children": children,
    }


def populate_entries(
    ipr_stg_uri: str,
    clans_file: str,
    entries_file: str,
    overlapping_file: str,
    xrefs_file: str,
    structures_file: str,
    pfam_file: str,
    intact_file: str,
):
    logger.info("loading clan members")
    entries_in_clan = {}
    with open(clans_file, "rb") as fh:
        for clan in pickle.load(fh).values():
            for entry_acc, _, _, _, _ in clan["members"]:
                entries_in_clan[entry_acc] = {
                    "accession": clan["accession"],
                    "name": clan["name"],
                }

    logger.info("loading structures")
    highres_structures = {}
    with open(structures_file, "rb") as fh:
        for s in pickle.load(fh).values():
            if (
                s["resolution"] is not None
                and s["resolution"] <= _REPR_STRUCT_MAX_RESOLUTION
            ):
                highres_structures[s["id"]] = (s["name"], s["resolution"])

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
                others.append(
                    {
                        "accession": other.accession,
                        "name": other.name,
                        "type": other.type.lower(),
                    }
                )

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

    with open(pfam_file, "rb") as fh:
        pfam_families = pickle.load(fh)

    logger.info("parsing IntAct data")
    intact_data = intact.get_interpro_interactions(intact_file)

    logger.info("creating table")
    con = connect(ipr_stg_uri)
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
            member_databases JSONB,
            integrated_id VARCHAR(30),
            go_terms JSONB,
            description JSONB,
            wikipedia JSONB,
            details JSONB,
            literature JSONB,
            hierarchy JSONB,
            cross_references JSONB,
            interactions JSONB,
            pathways JSONB,
            overlaps_with JSONB,
            is_llm BOOLEAN NOT NULL,
            is_reviewed_llm BOOLEAN NOT NULL,
            is_updated_llm BOOLEAN NOT NULL,
            is_public BOOLEAN NOT NULL,
            history JSONB,
            entry_date TIMESTAMP NOT NULL,
            deletion_date TIMESTAMP,
            set_info JSONB,
            representative_structure JSONB,
            counts JSONB NOT NULL
        )
        """
    )

    query = """
        INSERT INTO webfront_entry
        VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s,
          %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
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
                    pathways[key] = [
                        dict(zip(("id", "name"), item)) for item in xrefs[key]
                    ]

            history = {}
            if entry.old_names:
                history["names"] = entry.old_names

            if entry.old_short_names:
                history["short_names"] = entry.old_short_names

            if entry.old_integrations:
                # Convert DB name to lower cases (API/client relies on LC)
                history["signatures"] = {
                    k.lower(): v for k, v in entry.old_integrations.items()
                }

            # Force keys of cross-references to lower case
            for key in list(entry.cross_references.keys()):
                value = entry.cross_references.pop(key)
                entry.cross_references[key.lower()] = value

            best_coverage = _REPR_STRUCT_MIN_COVERAGE
            best_resolution = _REPR_STRUCT_MAX_RESOLUTION
            best_structure = None
            for pdb_id, coverage in xrefs["structures"]:
                try:
                    s_name, s_reso = highres_structures[pdb_id]
                except KeyError:
                    continue

                if coverage < best_coverage:
                    continue
                elif coverage > best_coverage or s_reso < best_resolution:
                    best_coverage = coverage
                    best_resolution = s_reso
                    best_structure = {"accession": pdb_id, "name": s_name}

            entry_hierarchy, num_subfamilies = get_hierarchy(entry, hierarchy)
            entry_clan = entries_in_clan.get(entry.accession)
            ppi = intact_data.get(entry.accession, [])

            pfam_family = pfam_families.get(entry.accession, {})
            pfam_details = pfam_family.get("details")
            pfam_wiki = pfam_family.get("wikipedia", [])

            records.append(
                (
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
                    jsonify(pfam_wiki, nullable=True),
                    jsonify(pfam_details, nullable=True),
                    jsonify(entry.literature, nullable=True),
                    jsonify(entry_hierarchy, nullable=True),
                    jsonify(entry.cross_references, nullable=True),
                    jsonify(ppi, nullable=True),
                    jsonify(pathways, nullable=True),
                    jsonify(overlaps_with.get(entry.accession, []), nullable=True),
                    entry.llm,
                    entry.llm_reviewed,
                    entry.llm_updated,
                    entry.public,
                    jsonify(history, nullable=True),
                    entry.creation_date,
                    entry.deletion_date,
                    jsonify(entry_clan, nullable=True),
                    jsonify(best_structure),
                    jsonify(
                        {
                            "subfamilies": num_subfamilies,
                            "domain_architectures": len(xrefs["dom_orgs"]),
                            "interactions": len(ppi),
                            "matches": xrefs["matches"],
                            "pathways": sum([len(v) for v in pathways.values()]),
                            "proteins": len(xrefs["proteins"]),
                            "proteomes": len(xrefs["proteomes"]),
                            "sets": 1 if entry_clan else 0,
                            "structural_models": xrefs["struct_models"],
                            "structures": len(xrefs["structures"]),
                            "taxa": len(xrefs["taxa"]["all"]),
                        },
                        nullable=False,
                    ),
                )
            )

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
                k.lower(): v for k, v in entry.old_integrations.items()
            }

        pfam_family = pfam_families.get(entry.accession, {})
        pfam_details = pfam_family.get("details")
        pfam_wiki = pfam_family.get("wikipedia", [])

        entry_clan = entries_in_clan.get(entry.accession)
        entry_hierarchy, num_subfamilies = get_hierarchy(entry, hierarchy)
        ppi = intact_data.get(entry.accession, [])
        records.append(
            (
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
                jsonify(pfam_wiki, nullable=True),
                jsonify(pfam_details, nullable=True),
                jsonify(entry.literature, nullable=True),
                jsonify(entry_hierarchy, nullable=True),
                jsonify(entry.cross_references, nullable=True),
                jsonify(ppi, nullable=True),
                jsonify(pathways, nullable=True),
                jsonify(overlaps_with.get(entry.accession, []), nullable=True),
                entry.llm,
                entry.llm_reviewed,
                entry.llm_updated,
                entry.public,
                jsonify(history, nullable=True),
                entry.creation_date,
                entry.deletion_date,
                jsonify(entry_clan, nullable=True),
                None,
                jsonify(
                    {
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
                    },
                    nullable=False,
                ),
            )
        )

    for i in range(0, len(records), 1000):
        cur.executemany(query, records[i : i + 1000])

    con.commit()
    cur.close()
    con.close()

    logger.info("done")


def index_entries(uri: str):
    con = connect(uri)
    cur = con.cursor()
    cur.execute(
        """
        CREATE INDEX IF NOT EXISTS i_entry_database
        ON webfront_entry (source_database)
        """
    )
    cur.execute(
        """
        CREATE INDEX IF NOT EXISTS i_entry_integrated
        ON webfront_entry (integrated_id)
        """
    )
    cur.execute(
        """
        CREATE INDEX IF NOT EXISTS i_entry_name
        ON webfront_entry (name)
        """
    )
    cur.execute(
        """
        CREATE INDEX IF NOT EXISTS i_entry_short_name
        ON webfront_entry (short_name)
        """
    )
    cur.execute(
        """
        CREATE INDEX IF NOT EXISTS i_entry_deletion_date
        ON webfront_entry (deletion_date)
        """
    )
    con.commit()
    cur.close()
    con.close()


def populate_entry_taxa_distrib(uri: str, entries_file: str, xrefs_file: str):
    logger.info("loading entries")
    with open(entries_file, "rb") as fh:
        entries: dict[str, Entry] = pickle.load(fh)

    logger.info("creating table")
    con = connect(uri)
    cur = con.cursor()
    cur.execute("DROP TABLE IF EXISTS webfront_entrytaxa")
    cur.execute(
        """
        CREATE TABLE webfront_entrytaxa
        (
            accession VARCHAR(30) PRIMARY KEY NOT NULL,
            tree TEXT
        )
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
