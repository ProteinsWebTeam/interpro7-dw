import json
import os
from typing import Dict, List, Optional

from . import entry, reduce
from .. import is_overlapping, repr_frag
from ... import dbms, io, logger, pdbe


def insert_structures(ora_uri, uri):
    structures = pdbe.get_structures(ora_uri)
    sec_structures = pdbe.get_secondary_structures(ora_uri)

    con, cur = dbms.connect(uri)
    cur.close()
    table = dbms.Populator(
        con=con,
        query="""
            INSERT INTO webfront_structure (
                  accession,
                  name,
                  source_database,
                  experiment_type,
                  release_date,
                  resolution,
                  literature,
                  chains,
                  proteins,
                  secondary_structures
              ) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
        """
    )

    for pdb_id, s in structures.items():
        table.insert((
            pdb_id,
            s["name"],
            "pdb",
            s["evidence"],
            s["date"],
            s["resolution"],
            json.dumps(s["citations"]),
            json.dumps(sorted({chain_id
                               for chains in s["proteins"].values()
                               for chain_id in chains})),
            json.dumps(s["proteins"]),
            json.dumps(sec_structures.get(pdb_id, []))
        ))
    table.close()
    con.commit()
    con.close()


def _is_entry_overlapping(fragments: List[dict],
                          chains: Dict[str, List[dict]]) -> bool:
    fragments.sort(key=repr_frag)
    start = fragments[0]["start"]
    end = fragments[-1]["end"]

    for locations in chains.values():
        for loc in locations:
            if is_overlapping(start, end, loc["protein_start"], loc["protein_end"]):
                return True

    return False


def get_structures(uri: str) -> dict:
    structures = {}
    con, cur = dbms.connect(uri)
    cur.execute(
        """
        SELECT accession, name, experiment_type, release_date, 
               resolution, literature, proteins
        FROM webfront_structure
        """
    )

    for row in cur:
        structures[row[0]] = {
            "accession": row[0],
            "name": row[1],
            "evidence": row[2],
            "date": row[3],
            "resolution": row[4],
            "citations": json.loads(row[5]),
            "proteins": json.loads(row[6])
        }

    cur.close()
    con.close()

    return structures


def update_counts(my_uri: str, src_proteins: str, src_proteomes:str,
                  src_matches: str, src_ida: str, processes: int=1,
                  sync_frequency: int=100000, tmpdir: Optional[str]=None):
    logger.info("starting")
    if tmpdir:
        os.makedirs(tmpdir, exist_ok=True)

    # Get required MySQL data
    entries = entry.get_entries(my_uri)
    entry2set = entry.get_entry2set(my_uri)
    pdb_ids = set()
    protein2structures = {}
    for pdb_id, s in get_structures(my_uri).items():
        pdb_ids.add(pdb_id)
        for protein_acc, chains in s["proteins"].items():
            try:
                protein = protein2structures[protein_acc]
            except KeyError:
                protein = protein2structures[protein_acc] = {}
            finally:
                protein[pdb_id] = chains

    # Open existing stores containing protein-related info
    proteins = io.Store(src_proteins)
    protein2proteome = io.Store(src_proteomes)
    protein2matches = io.Store(src_matches)
    protein2ida = io.Store(src_ida)

    protein_counts = {}
    cnt_proteins = 0
    cnt_updates = 0
    with io.Store(keys=io.Store.chunk_keys(pdb_ids, 100), tmpdir=tmpdir) as store:
        for protein_acc, p in proteins:
            cnt_proteins += 1
            if not cnt_proteins % 10000000:
                logger.info(f"{cnt_proteins:>12}")

            try:
                structures = protein2structures[protein_acc]
            except KeyError:
                continue

            xrefs = {
                "domain_architectures": set(),
                "entries": {},
                "proteomes": set(),
                "sets": set(),
                "taxa": {p["taxon"]}
            }

            try:
                upid = protein2proteome[protein_acc]
            except KeyError:
                pass
            else:
                xrefs["proteomes"].add(upid)

            try:
                ida, ida_id = protein2ida[protein_acc]
            except KeyError:
                pass
            else:
                xrefs["domain_architectures"].add(ida)

            matches = {}
            for match in protein2matches.get(protein_acc, []):
                method_acc = match["method_ac"]
                try:
                    matches[method_acc] += match["fragments"]
                except KeyError:
                    matches[method_acc] = list(match["fragments"])

            for pdb_id, chains in structures.items():
                _xrefs = xrefs.copy()

                # Only count entries that are overlapping
                for method_acc, fragments in matches.items():
                    if _is_entry_overlapping(fragments, chains):
                        database = entries[method_acc]["database"]
                        try:
                            _xrefs["entries"][database].add(method_acc)
                        except KeyError:
                            _xrefs["entries"][database] = {method_acc}

                        entry_acc = entries[method_acc]["integrated"]
                        if entry_acc:
                            database = "interpro"
                            try:
                                _xrefs["entries"][database].add(entry_acc)
                            except KeyError:
                                _xrefs["entries"][database] = {entry_acc}

                # Adding sets for overlapping entries
                for accessions in _xrefs["entries"].values():
                    for entry_acc in accessions:
                        try:
                            set_acc = entry2set[entry_acc]
                        except KeyError:
                            pass
                        else:
                            _xrefs["sets"].add(set_acc)

                store.update(pdb_id, _xrefs)

                if pdb_id in protein_counts:
                    protein_counts[pdb_id] += 1
                else:
                    protein_counts[pdb_id] = 1

            cnt_updates += 1
            if not cnt_updates % sync_frequency:
                store.sync()

        proteins.close()
        protein2proteome.close()
        protein2matches.close()
        protein2ida.close()
        logger.info(f"{cnt_proteins:>12}")

        for pdb_id, cnt in protein_counts.items():
            store.update(pdb_id, {"proteins": cnt})
            pdb_ids.remove(pdb_id)

        # Remaining structures
        for pdb_id in pdb_ids:
            store.update(pdb_id, {
                "domain_architectures": set(),
                "entries": {},
                "proteomes": set(),
                "proteins": set(),
                "sets": set(),
                "taxa": set(),
            })

        size = store.merge(processes=processes)
        logger.info("Disk usage: {:.0f}MB".format(size / 1024 ** 2))

        con, cur = dbms.connect(my_uri)
        cur.close()
        query = "UPDATE webfront_structure SET counts = %s WHERE accession = %s"
        with dbms.Populator(con, query) as table:
            for upid, _xrefs in store:
                counts = reduce(_xrefs)
                counts["entries"]["total"] = sum(
                    counts["entries"].values())
                table.update((json.dumps(counts), upid))

        con.commit()
        con.close()

    logger.info("complete")
