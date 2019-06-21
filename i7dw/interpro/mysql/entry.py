import json

from . import reduce
from .. import oracle
from ... import dbms, cdd, logger, pfam
from ...io import KVdb, Store, TempFile


def jsonify(x):
    return json.dumps(x) if x else None


def parse(value, default=None):
    if value is None:
        return default
    else:
        return json.loads(value)


def find_node(node, accession, relations=[]):
    """
    Find a entry (node) in its hierarchy tree
    and store its relations (ancestors + direct children)
    """
    if node['accession'] == accession:
        relations += [child['accession'] for child in node['children']]
        return node

    for child in node['children']:
        child = find_node(child, accession, relations)

        if child:
            relations.append(node['accession'])
            return child

    return None


def insert_entries(ora_uri, pfam_uri, my_uri, chunk_size=100000):
    wiki = pfam.get_wiki(pfam_uri)

    data = []
    for e in oracle.get_entries(ora_uri):
        data.append((
            e["accession"],
            e["type"],
            e["name"],
            e["short_name"],
            e["database"],
            jsonify(e["member_databases"]),
            e["integrated"],
            jsonify(e["go_terms"]),
            jsonify(e["descriptions"]),
            jsonify(wiki.get(e["accession"])),
            jsonify(e["citations"]),
            jsonify(e["hierarchy"]),
            jsonify(e["cross_references"]),
            1,      # is alive
            e["date"],
            None    # deletion date
        ))

    for e in oracle.get_deleted_entries(ora_uri):
        data.append((
            e["accession"],
            e["type"],
            e["name"],
            e["short_name"],
            e["database"],
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            0,
            e["creation_date"],
            e["deletion_date"]
        ))

    con, cur = dbms.connect(my_uri)

    for i in range(0, len(data), chunk_size):
        cur.executemany(
            """
              INSERT INTO webfront_entry (
                accession,
                type,
                name,
                short_name,
                source_database,
                member_databases,
                integrated_id,
                go_terms,
                description,
                wikipedia,
                literature,
                hierarchy,
                cross_references,
                is_alive,
                entry_date,
                deletion_date
              )
              VALUES (
                %s, %s, %s, %s, %s, %s, %s, %s, %s, %s,
                %s, %s, %s, %s, %s, %s
              )
            """,
            data[i:i+chunk_size]
        )

    cur.close()
    con.commit()
    con.close()


def get_entries(uri: str, alive_only: bool=True) -> dict:
    query = """
        SELECT
            accession, source_database, entry_date, description,
            integrated_id, name, type, short_name, member_databases,
            go_terms, literature, cross_references, hierarchy
        FROM webfront_entry
    """

    if alive_only:
        query += "WHERE is_alive = 1"

    con, cur = dbms.connect(uri)
    cur.execute(query)

    entries = {}
    for row in cur:
        accession = row[0]
        relations = []
        hierarchy = parse(row[12])
        if hierarchy:
            find_node(hierarchy, accession, relations)
            _hierarchy = hierarchy.get("accession")
        else:
            _hierarchy = None

        entries[accession] = {
            "accession": accession,
            "database": row[1],
            "date": row[2],
            "descriptions": parse(row[3], []),
            "integrated": row[4],
            "name": row[5],
            "type": row[6],
            "short_name": row[7],
            "member_databases": parse(row[8], {}),
            "go_terms": parse(row[9], []),
            "citations": parse(row[10], {}),
            "cross_references": parse(row[11], {}),
            "root": _hierarchy,
            "relations": relations
        }

    cur.close()
    con.close()

    return entries


def insert_annotations(pfam_uri, my_uri, chunk_size=10000):
    data = pfam.get_annotations(pfam_uri)

    con, cur = dbms.connect(my_uri)
    for i in range(0, len(data), chunk_size):
        cur.executemany(
            """
            INSERT INTO webfront_entryannotation (
                annotation_id,
                accession_id,
                type,
                value,
                mime_type
            ) VALUES (%s, %s, %s, %s, %s)
            """,
            data[i:i+chunk_size]
        )

    cur.close()
    con.commit()
    con.close()


def insert_sets_new(ora_uri, pfam_uri, my_uri):
    sets = pfam.get_clans(pfam_uri)
    sets.update(cdd.get_superfamilies())

    con, cur = dbms.connect(my_uri)
    cur.close()
    table1 = dbms.Populator(
        uri=my_uri,
        query="INSERT INTO webfront_set (accession, name, description, "
              "source_database, is_set, relationships) "
              "VALUES (%s, %s, %s, %s, %s, %s)"
    )
    table2 = dbms.Populator(
        uri=my_uri,
        query="INSERT INTO webfront_alignment (set_acc, entry_acc, "
              "target_acc, target_set_acc, score, seq_length, domains) "
              "VALUES (%s, %s, %s, %s, %s, %s, %s)"
    )
    for set_ac, db, rels, alns in oracle.get_profile_alignments_new(ora_uri):
        if db in ("cdd", "pfam"):
            try:
                s = sets[set_ac]
            except KeyError:
                logger.warning("unknown CDD/Pfam set: {}".format(set_ac))
                continue

            if db == "pfam":
                # Use nodes from Pfam DB (they have a 'real' score)
                rels["nodes"] = s["relationships"]["nodes"]

            name = s["name"] or set_ac
            desc = s["description"]
        else:
            name = set_ac
            desc = None

        table1.insert((set_ac, name, desc, db, 1, json.dumps(rels)))
        for entry_acc, targets in alns.items():
            for target in targets:
                table2.insert((set_ac, entry_acc, *target))

    table1.close()
    table2.close()


def insert_sets(ora_uri, pfam_uri, my_uri, tmpdir=None):
    with TempFile(dir=tmpdir) as f:
        logger.info("Pfam clans")
        sets = pfam.get_clans(pfam_uri)
        gen = oracle.get_profile_alignments(ora_uri, "pfam")
        for set_ac, relationships, alignments in gen:
            try:
                clan = sets[set_ac]
            except KeyError:
                logger.warning("unknown Pfam clan: {}".format(set_ac))
                continue

            # Use nodes from Pfam DB for the score
            relationships["nodes"] = clan["relationships"]["nodes"]

            f.write((
                set_ac,
                clan["name"] or set_ac,
                clan["description"],
                "pfam",
                1,
                relationships,
                alignments
            ))

        logger.info("CDD superfamilies")
        sets = cdd.get_superfamilies()
        gen = oracle.get_profile_alignments(ora_uri, "cdd")
        for set_ac, relationships, alignments in gen:
            try:
                supfam = sets[set_ac]
            except KeyError:
                logger.warning("unknown CDD superfamily: {}".format(set_ac))
                continue

            f.write((
                set_ac,
                supfam["name"] or set_ac,
                supfam["description"],
                "cdd",
                1,
                relationships,
                alignments
            ))

        logger.info("PANTHER superfamilies")
        gen = oracle.get_profile_alignments(ora_uri, "panther")
        for set_ac, relationships, alignments in gen:
            f.write((
                set_ac,
                set_ac,
                None,
                "panther",
                1,
                relationships,
                alignments
            ))

        logger.info("PIRSF superfamilies")
        gen = oracle.get_profile_alignments(ora_uri, "pirsf")
        for set_ac, relationships, alignments in gen:
            f.write((
                set_ac,
                set_ac,
                None,
                "pirsf",
                1,
                relationships,
                alignments
            ))

        logger.info("{:,} bytes".format(f.size))

        con, cur = dbms.connect(my_uri)
        for acc, name, descr, db, is_set, relationships, alignments in f:
            cur.execute(
                """
                INSERT INTO webfront_set (
                  accession, name, description, source_database, is_set,
                  relationships
                ) VALUES (%s, %s, %s, %s, %s, %s)
                """,
                (acc, name, descr, db, is_set, json.dumps(relationships))
            )

            for entry_acc, targets in alignments.items():
                for target_acc, t in targets.items():
                    cur.execute(
                        """
                        INSERT INTO webfront_alignment (
                            set_acc, entry_acc, target_acc, target_set_acc, 
                            score, seq_length, domains
                        ) VALUES (%s, %s, %s, %s, %s, %s, %s)
                        """,
                        (acc, entry_acc, target_acc, t["set_acc"], t["score"],
                         t["length"], json.dumps(t["domains"]))
                    )

        con.commit()
        cur.execute(
            """
            CREATE INDEX i_webfront_alignment
            ON webfront_alignment (set_acc, entry_acc)
            """
        )
        cur.execute("ANALYZE TABLE webfront_set")
        cur.execute("ANALYZE TABLE webfront_alignment")
        cur.close()
        con.close()


def get_sets(uri: str) -> dict:
    con, cur = dbms.connect(uri, sscursor=True)
    cur.execute(
        """
        SELECT accession, name, description, source_database, relationships
        FROM webfront_set
        """
    )

    sets = {}
    for acc, name, description, database, relationships in cur:
        sets[acc] = {
            "name": name,
            "description": description,
            "database": database,
            "members": [
                n["accession"]
                for n in json.loads(relationships)["nodes"]
            ]
        }

    cur.close()
    con.close()

    return sets


def update_counts(uri: str, src_entries: str, tmpdir: str=None):
    logger.info("updating webfront_entry")
    sets = get_sets(uri)
    entry2set = {}
    for set_ac, s in sets.items():
        for entry_ac in s["members"]:
            entry2set[entry_ac] = set_ac

    all_entries = set(get_entries(uri, alive_only=False))
    cnt = 0
    with KVdb(cache_size=10, dir=tmpdir) as entries:
        con, cur = dbms.connect(uri)

        with Store(src_entries) as store:
            for entry_ac, xrefs in store:
                all_entries.remove(entry_ac)

                counts = reduce(xrefs)
                counts["matches"] = xrefs["matches"].pop()  # set of one item
                if entry2set.get(entry_ac):
                    entries[entry_ac] = xrefs
                    counts["sets"] = 1
                else:
                    counts["sets"] = 0

                cur.execute(
                    """
                    UPDATE webfront_entry
                    SET counts = %s
                    WHERE accession = %s
                    """,
                    (json.dumps(counts), entry_ac)
                )

                cnt += 1
                if not cnt % 10000:
                    logger.info("{:>9,}".format(cnt))

        for entry_ac in all_entries:
            cur.execute(
                """
                UPDATE webfront_entry
                SET counts = %s
                WHERE accession = %s
                """,
                (json.dumps({
                    "matches": 0,
                    "proteins": 0,
                    "proteomes": 0,
                    "sets": 1 if entry_ac in entry2set else 0,
                    "structures": 0,
                    "taxa": 0
                }), entry_ac)
            )

            cnt += 1
            if not cnt % 10000:
                logger.info("{:>9,}".format(cnt))

        logger.info("{:>9,}".format(cnt))

        logger.info("updating webfront_set")
        cnt = 0
        for set_ac, s in sets.items():
            xrefs = {
                "domain_architectures": set(),
                "entries": {
                    s["database"]: len(s["members"]),
                    "total": len(s["members"])
                },
                "proteins": set(),
                "proteomes": set(),
                "structures": set(),
                "taxa": set()
            }

            for entry_ac in s["members"]:
                try:
                    entry = entries[entry_ac]
                except KeyError:
                    continue

                for _type in ("domain_architectures", "proteins", "proteomes", "structures", "taxa"):
                    try:
                        accessions = entry[_type]
                    except KeyError:
                        # Type absent
                        continue
                    else:
                        xrefs[_type] |= accessions

            cur.execute(
                """
                UPDATE webfront_set
                SET counts = %s
                WHERE accession = %s
                """,
                (json.dumps(reduce(xrefs)), set_ac)
            )

            cnt += 1
            if not cnt % 1000:
                logger.info("{:>6,}".format(cnt))

        con.commit()
        cur.close()
        con.close()
        logger.info("{:>6,}".format(cnt))
        logger.info("database size: {:,}".format(entries.getsize()))
