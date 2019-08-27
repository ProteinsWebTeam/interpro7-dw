import json
import os
from typing import Dict, Optional

from . import reduce, structure
from .. import oracle
from ... import dbms, cdd, logger, pfam
from ...io import KVdb, Store


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
            "descriptions": parse(row[3], default=[]),
            "integrated": row[4],
            "name": row[5],
            "type": row[6],
            "short_name": row[7],
            "member_databases": parse(row[8], default={}),
            "go_terms": parse(row[9], default=[]),
            "citations": parse(row[10], default={}),
            "cross_references": parse(row[11], default={}),
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


def insert_sets(ora_uri, pfam_uri, my_uri):
    sets = pfam.get_clans(pfam_uri)
    sets.update(cdd.get_superfamilies())

    con, cur = dbms.connect(my_uri)
    query = """
        INSERT INTO webfront_set (accession, name, description,
          source_database, is_set, relationships
        ) VALUES (%s, %s, %s, %s, %s, %s)
    """
    table = dbms.Populator(
        con=con,
        query="INSERT INTO webfront_alignment (set_acc, entry_acc, "
              "target_acc, target_set_acc, score, seq_length, domains) "
              "VALUES (%s, %s, %s, %s, %s, %s, %s)"
    )
    for set_ac, db, rels, alns in oracle.get_profile_alignments(ora_uri):
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

        cur.execute(query, (set_ac, name, desc, db, 1, json.dumps(rels)))
        # con.commit()  # because of FK constraint
        for entry_acc, targets in alns.items():
            for target in targets:
                table.insert((set_ac, entry_acc, *target))

    table.close()
    con.commit()
    cur.close()
    con.close()


# def insert_sets(ora_uri, pfam_uri, my_uri, tmpdir=None):
#     with TempFile(dir=tmpdir) as f:
#         logger.info("Pfam clans")
#         sets = pfam.get_clans(pfam_uri)
#         gen = oracle.get_profile_alignments(ora_uri, "pfam")
#         for set_ac, relationships, alignments in gen:
#             try:
#                 clan = sets[set_ac]
#             except KeyError:
#                 logger.warning("unknown Pfam clan: {}".format(set_ac))
#                 continue
#
#             # Use nodes from Pfam DB for the score
#             relationships["nodes"] = clan["relationships"]["nodes"]
#
#             f.write((
#                 set_ac,
#                 clan["name"] or set_ac,
#                 clan["description"],
#                 "pfam",
#                 1,
#                 relationships,
#                 alignments
#             ))
#
#         logger.info("CDD superfamilies")
#         sets = cdd.get_superfamilies()
#         gen = oracle.get_profile_alignments(ora_uri, "cdd")
#         for set_ac, relationships, alignments in gen:
#             try:
#                 supfam = sets[set_ac]
#             except KeyError:
#                 logger.warning("unknown CDD superfamily: {}".format(set_ac))
#                 continue
#
#             f.write((
#                 set_ac,
#                 supfam["name"] or set_ac,
#                 supfam["description"],
#                 "cdd",
#                 1,
#                 relationships,
#                 alignments
#             ))
#
#         logger.info("PANTHER superfamilies")
#         gen = oracle.get_profile_alignments(ora_uri, "panther")
#         for set_ac, relationships, alignments in gen:
#             f.write((
#                 set_ac,
#                 set_ac,
#                 None,
#                 "panther",
#                 1,
#                 relationships,
#                 alignments
#             ))
#
#         logger.info("PIRSF superfamilies")
#         gen = oracle.get_profile_alignments(ora_uri, "pirsf")
#         for set_ac, relationships, alignments in gen:
#             f.write((
#                 set_ac,
#                 set_ac,
#                 None,
#                 "pirsf",
#                 1,
#                 relationships,
#                 alignments
#             ))
#
#         logger.info("{:,} bytes".format(f.size))
#
#         con, cur = dbms.connect(my_uri)
#         for acc, name, descr, db, is_set, relationships, alignments in f:
#             cur.execute(
#                 """
#                 INSERT INTO webfront_set (
#                   accession, name, description, source_database, is_set,
#                   relationships
#                 ) VALUES (%s, %s, %s, %s, %s, %s)
#                 """,
#                 (acc, name, descr, db, is_set, json.dumps(relationships))
#             )
#
#             for entry_acc, targets in alignments.items():
#                 for target_acc, t in targets.items():
#                     cur.execute(
#                         """
#                         INSERT INTO webfront_alignment (
#                             set_acc, entry_acc, target_acc, target_set_acc,
#                             score, seq_length, domains
#                         ) VALUES (%s, %s, %s, %s, %s, %s, %s)
#                         """,
#                         (acc, entry_acc, target_acc, t["set_acc"], t["score"],
#                          t["length"], json.dumps(t["domains"]))
#                     )
#
#         con.commit()
#         cur.execute(
#             """
#             CREATE INDEX i_webfront_alignment
#             ON webfront_alignment (set_acc, entry_acc)
#             """
#         )
#         cur.execute("ANALYZE TABLE webfront_set")
#         cur.execute("ANALYZE TABLE webfront_alignment")
#         cur.close()
#         con.close()


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


def _export(my_uri: str, src_proteins: str, src_proteomes:str,
            src_matches: str, src_ida: str, dst_xrefs: str,
            sync_frequency: int, tmpdir: Optional[str]) -> Store:
    # Get required MySQL data
    entries = get_entries(my_uri)
    accessions = set(entries)
    protein2structures = structure.get_protein2structures(my_uri)
    entry2set = get_entry2set(my_uri)

    # Open existing stores containing protein-related info
    proteins = Store(src_proteins)
    protein2proteome = Store(src_proteomes)
    protein2matches = Store(src_matches)
    protein2ida = Store(src_ida)

    xrefs = Store(dst_xrefs, Store.chunk_keys(accessions, 10), tmpdir)
    entry_match_counts = {
        "mobidb-lite": 0
    }
    cnt_proteins = 0
    mobidb_proteins = 0
    for protein_acc, p in proteins:
        cnt_proteins += 1
        if not cnt_proteins % 10000000:
            logger.info(f"{cnt_proteins:>12}")

        tax_id = p["taxon"]
        matches = {}
        for match in protein2matches.get(protein_acc, []):
            method_acc = match["method_ac"]
            entry_acc = entries[method_acc]["integrated"]

            if method_acc in matches:
                matches[method_acc] += 1
            else:
                matches[method_acc] = 1

            if entry_acc:
                if entry_acc in matches:
                    matches[entry_acc] += 1
                else:
                    matches[entry_acc] = 1

        _xrefs = {
            "domain_architectures": set(),
            "proteomes": set(),
            "structures": set(),
            "taxa": {tax_id}
        }

        try:
            ida, ida_id = protein2ida[protein_acc]
        except KeyError:
            pass
        else:
            _xrefs["domain_architectures"].add(ida)

        try:
            upid = protein2proteome[protein_acc]
        except KeyError:
            pass
        else:
            _xrefs["proteomes"].add(upid)

        try:
            pdbe_ids = protein2structures[protein_acc]
        except KeyError:
            pass
        else:
            _xrefs["structures"] = pdbe_ids

        """
        As MobiDB-lite is not included in EBISearch,
        we do not need to keep track of matched proteins.
        We only need the number of proteins matched.
        """
        try:
            cnt = matches.pop("mobidb-lite")
        except KeyError:
            pass
        else:
            mobidb_proteins += 1
            entry_match_counts["mobidb-lite"] += cnt

        # Add proteins for other entries
        _xrefs["proteins"] = {(protein_acc, p["identifier"])}

        for entry_acc, cnt in matches.items():
            xrefs.update(entry_acc, _xrefs)
            if entry_acc in entry_match_counts:
                entry_match_counts[entry_acc] += cnt
            else:
                entry_match_counts[entry_acc] = cnt

        if not cnt_proteins % sync_frequency:
            xrefs.sync()

    proteins.close()
    protein2proteome.close()
    protein2matches.close()
    protein2ida.close()
    logger.info(f"{cnt_proteins:>12}")

    # Add protein count for MobiDB-lite
    xrefs.update("mobidb-lite", {"proteins": mobidb_proteins})

    # Add match counts and set for entries *with* protein matches
    for entry_acc, cnt in entry_match_counts.items():
        _xrefs = {"matches": cnt}
        try:
            set_acc = entry2set[entry_acc]
        except KeyError:
            _xrefs["sets"] = set()
        else:
            _xrefs["sets"] = {set_acc}
        finally:
            xrefs.update(entry_acc, _xrefs)
            accessions.remove(entry_acc)

    # Remaining entries without protein matches
    for entry_acc in accessions:
        _xrefs = {
            "domain_architectures": set(),
            "matches": 0,
            "proteins": set(),
            "proteomes": set(),
            "structures": set(),
            "taxa": set()
        }

        try:
            set_acc = entry2set[entry_acc]
        except KeyError:
            _xrefs["sets"] = set()
        else:
            _xrefs["sets"] = {set_acc}
        finally:
            xrefs.update(entry_acc, _xrefs)

    return xrefs


def export_xrefs(my_uri: str, src_proteins: str, src_proteomes:str,
                 src_matches: str, src_ida: str, dst: str, processes: int=1,
                 sync_frequency: int=100000, tmpdir: Optional[str]=None):
    logger.info("starting")
    if tmpdir:
        os.makedirs(tmpdir, exist_ok=True)

    xrefs = _export(my_uri, src_proteins, src_proteomes, src_matches, src_ida,
                    dst, sync_frequency, tmpdir)
    size = xrefs.merge(processes=processes)
    logger.info("disk usage: {:.0f}MB".format(size/1024**2))


def update_counts(my_uri: str, src_entries: str, tmpdir: Optional[str]=None):
    logger.info("starting")
    if tmpdir:
        os.makedirs(tmpdir, exist_ok=True)

    entry2set = get_entry2set(my_uri)

    con, cur = dbms.connect(my_uri)
    cur.close()
    with KVdb(dir=tmpdir) as kvdb:
        logger.info("updating webfront_entry")
        query = "UPDATE webfront_entry SET counts = %s WHERE accession = %s"
        with dbms.Populator(con, query) as table, Store(src_entries) as xrefs:
            for entry_acc, _xrefs in xrefs:
                table.update((json.dumps(reduce(_xrefs)), entry_acc))

                if entry_acc in entry2set:
                    kvdb[entry_acc] = _xrefs

        logger.info("updating webfront_set")
        query = "UPDATE webfront_set SET counts = %s WHERE accession = %s"
        with dbms.Populator(con, query) as table:
            for set_acc, s in get_sets(my_uri).items():
                _xrefs = {
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
                for entry_acc in s["members"]:
                    try:
                        xrefs = kvdb[entry_acc]
                    except KeyError:
                        continue
                    else:
                        for key in xrefs:
                            try:
                                _xrefs[key] |= xrefs[key]
                            except KeyError:
                                pass

                table.update((json.dumps(reduce(_xrefs)), set_acc))

        logger.info("disk usage: {:.0f}MB".format(kvdb.size/1024**2))

    con.commit()
    con.close()
    logger.info("complete")


def get_entry2set(my_uri: str) -> Dict[str, str]:
    entry2set = {}
    for set_acc, s in get_sets(my_uri).items():
        for entry_acc in s["members"]:
            entry2set[entry_acc] = set_acc

    return entry2set
