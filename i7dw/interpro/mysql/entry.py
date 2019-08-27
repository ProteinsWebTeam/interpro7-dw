import json
import os
from typing import Dict, List, Optional

from . import reduce, structure
from .. import is_overlapping, oracle, repr_frag
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


def _is_structure_overlapping(start: int, end: int, chains: Dict[str, List[dict]]) -> bool:
    for locations in chains.values():
        for loc in locations:
            if is_overlapping(start, end, loc["protein_start"], loc["protein_end"]):
                return True
    return False


def _export(my_uri: str, src_proteins: str, src_proteomes:str,
            src_matches: str, src_ida: str, dst_xrefs: str,
            sync_frequency: int, tmpdir: Optional[str]) -> Store:
    # Get required MySQL data
    entries = get_entries(my_uri)
    accessions = set(entries)
    entry2set = get_entry2set(my_uri)
    protein2structures = {}
    for pdb_id, s in structure.get_structures(my_uri).items():
        for protein_acc, chains in s["proteins"].items():
            try:
                protein = protein2structures[protein_acc]
            except KeyError:
                protein = protein2structures[protein_acc] = {}
            finally:
                protein[pdb_id] = chains

    # Open existing stores containing protein-related info
    proteins = Store(src_proteins)
    protein2proteome = Store(src_proteomes)
    protein2matches = Store(src_matches)
    protein2ida = Store(src_ida)

    store = Store(dst_xrefs, Store.chunk_keys(accessions, 10), tmpdir)
    entry_match_counts = {
        "mobidb-lite": 0
    }
    mobidb_proteins = 0
    cnt_proteins = 0
    for protein_acc, p in proteins:
        cnt_proteins += 1
        if not cnt_proteins % 10000000:
            logger.info(f"{cnt_proteins:>12}")

        xrefs = {
            "domain_architectures": set(),
            "proteins": {(protein_acc, p["identifier"])},
            "proteomes": set(),
            "structures": set(),
            "taxa": {p["taxon"]}
        }

        try:
            ida, ida_id = protein2ida[protein_acc]
        except KeyError:
            pass
        else:
            xrefs["domain_architectures"].add(ida)

        try:
            upid = protein2proteome[protein_acc]
        except KeyError:
            pass
        else:
            xrefs["proteomes"].add(upid)

        structures = protein2structures.get(protein_acc, {})

        matches = {}
        match_counts = {}
        for match in protein2matches.get(protein_acc, []):
            method_acc = match["method_ac"]
            try:
                matches[method_acc] += match["fragments"]
            except KeyError:
                matches[method_acc] = list(match["fragments"])
                match_counts[method_acc] = 0
            finally:
                match_counts[method_acc] += 1

            entry_acc = entries[method_acc]["integrated"]
            if entry_acc:
                try:
                    match_counts[entry_acc] += 1
                except KeyError:
                    match_counts[entry_acc] = 1

        """
        As MobiDB-lite is not included in EBISearch,
        we do not need to keep track of matched proteins.
        We only need the number of proteins matched.
        """
        try:
            cnt = match_counts.pop("mobidb-lite")
        except KeyError:
            pass
        else:
            mobidb_proteins += 1
            entry_match_counts["mobidb-lite"] += cnt
            del matches["mobidb-lite"]

        for entry_acc, cnt in match_counts.items():
            try:
                entry_match_counts[entry_acc] += cnt
            except KeyError:
                entry_match_counts[entry_acc] = cnt

        for method_acc, fragments in matches.items():
            _xrefs = xrefs.copy()

            fragments.sort(key=repr_frag)
            start = fragments[0]["start"]
            end = fragments[-1]["end"]

            for pdb_id, chains in structures.items():
                if _is_structure_overlapping(start, end, chains):
                    _xrefs["structures"].add(pdb_id)

            store.update(method_acc, _xrefs)

            entry_acc = entries[method_acc]["integrated"]
            if entry_acc:
                store.update(entry_acc, _xrefs)

        if not cnt_proteins % sync_frequency:
            store.sync()

    proteins.close()
    protein2proteome.close()
    protein2matches.close()
    protein2ida.close()
    logger.info(f"{cnt_proteins:>12}")

    # Add protein count for MobiDB-lite
    store.update("mobidb-lite", {"proteins": mobidb_proteins})

    # Add match counts and set for entries *with* protein matches
    for entry_acc, cnt in entry_match_counts.items():
        xrefs = {"matches": cnt}
        try:
            set_acc = entry2set[entry_acc]
        except KeyError:
            xrefs["sets"] = set()
        else:
            xrefs["sets"] = {set_acc}
        finally:
            store.update(entry_acc, xrefs)
            accessions.remove(entry_acc)

    # Remaining entries without protein matches
    for entry_acc in accessions:
        xrefs = {
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
            xrefs["sets"] = set()
        else:
            xrefs["sets"] = {set_acc}
        finally:
            store.update(entry_acc, xrefs)

    return store


def export_xrefs(my_uri: str, src_proteins: str, src_proteomes:str,
                 src_matches: str, src_ida: str, dst: str, processes: int=1,
                 sync_frequency: int=100000, tmpdir: Optional[str]=None):
    logger.info("starting")
    if tmpdir:
        os.makedirs(tmpdir, exist_ok=True)

    store = _export(my_uri, src_proteins, src_proteomes, src_matches, src_ida,
                    dst, sync_frequency, tmpdir)
    size = store.merge(processes=processes)
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
        with dbms.Populator(con, query) as table, Store(src_entries) as store:
            for entry_acc, xrefs in store:
                table.update((json.dumps(reduce(xrefs)), entry_acc))

                if entry_acc in entry2set:
                    kvdb[entry_acc] = xrefs

        logger.info("updating webfront_set")
        query = "UPDATE webfront_set SET counts = %s WHERE accession = %s"
        with dbms.Populator(con, query) as table:
            for set_acc, s in get_sets(my_uri).items():
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
                for entry_acc in s["members"]:
                    try:
                        xrefs = kvdb[entry_acc]
                    except KeyError:
                        continue
                    else:
                        for key in xrefs:
                            try:
                                xrefs[key] |= xrefs[key]
                            except KeyError:
                                pass

                table.update((json.dumps(reduce(xrefs)), set_acc))

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
