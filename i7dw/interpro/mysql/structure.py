import json
import os
from typing import Dict, Optional, Set

from . import entry, reduce
from ... import dbms, logger, pdbe
from ...io import Store


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

    for pdbe_id, s in structures.items():
        proteins = {}
        all_chains = set()

        for acc, chains in s["proteins"].items():
            proteins[acc] = []
            for chain_id in chains:
                all_chains.add(chain_id)
                proteins[acc].append(chain_id)

        table.insert((
            pdbe_id,
            s["name"],
            "pdb",
            s["evidence"],
            s["date"],
            s["resolution"],
            json.dumps(s["citations"]),
            json.dumps(sorted(all_chains)),
            json.dumps(proteins),
            json.dumps(sec_structures.get(pdbe_id, []))
        ))
    table.close()
    con.commit()
    con.close()


def get_structures(uri: str) -> dict:
    structures = {}
    con, cur = dbms.connect(uri)
    cur.execute(
        """
        SELECT accession, name, experiment_type, resolution, proteins
        FROM webfront_structure
        """
    )

    for acc, name, _type, resolution, proteins in cur:
        structures[acc] = {
            "accession": acc,
            "name": name,
            "type": _type,
            "resolution": resolution,
            "proteins": json.loads(proteins)
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

    # Open existing stores containing protein-related info
    proteins = Store(src_proteins)
    protein2proteome = Store(src_proteomes)
    protein2matches = Store(src_matches)
    protein2ida = Store(src_ida)

    # Get required MySQL data
    entries = entry.get_entries(my_uri)
    protein2structures = get_protein2structures(my_uri)
    entry2set = entry.get_entry2set(my_uri)
    structures = get_structures(my_uri)

    protein_counts = {}
    cnt_proteins = 0
    cnt_updates = 0
    with Store(keys=Store.chunk_keys(structures, 100), tmpdir=tmpdir) as xrefs:
        for protein_acc, p in proteins:
            cnt_proteins += 1
            if not cnt_proteins % 10000000:
                logger.info(f"{cnt_proteins:>12}")

            try:
                pdbe_ids = protein2structures[protein_acc]
            except KeyError:
                continue

            prot_entries = set()
            for match in protein2matches.get(protein_acc, []):
                method_acc = match["method_ac"]
                prot_entries.add(method_acc)

                entry_acc = entries[method_acc]["integrated"]
                if entry_acc:
                    prot_entries.add(entry_acc)

            entry_databases = {}
            entry_sets = set()
            for entry_acc in prot_entries:
                database = entries[entry_acc]["database"]
                if database in entry_databases:
                    entry_databases[database].add(entry_acc)
                else:
                    entry_databases[database] = {entry_acc}

                try:
                    set_acc = entry2set[entry_acc]
                except KeyError:
                    pass
                else:
                    entry_sets.add(set_acc)

            _xrefs = {
                "domain_architectures": set(),
                "entries": entry_databases,
                "proteomes": set(),
                "sets": entry_sets,
                "taxa": {p["taxon"]}
            }

            try:
                upid = protein2proteome[protein_acc]
            except KeyError:
                pass
            else:
                _xrefs["proteomes"].add(upid)

            try:
                ida, ida_id = protein2ida[protein_acc]
            except KeyError:
                pass
            else:
                _xrefs["domain_architectures"].add(ida)

            for pdbe_id in pdbe_ids:
                xrefs.update(pdbe_id, _xrefs)

                if pdbe_id in protein_counts:
                    protein_counts[pdbe_id] += 1
                else:
                    protein_counts[pdbe_id] = 1

            cnt_updates += 1
            if not cnt_updates % sync_frequency:
                xrefs.sync()

        proteins.close()
        protein2proteome.close()
        protein2matches.close()
        protein2ida.close()
        logger.info(f"{cnt_proteins:>12}")

        for pdbe_id, cnt in protein_counts.items():
            xrefs.update(pdbe_id, {"proteins": cnt})
            structures.pop(pdbe_id)

        # Remaining structures
        for pdbe_id in structures:
            xrefs.update(pdbe_id, {
                "domain_architectures": set(),
                "entries": {},
                "proteomes": set(),
                "proteins": set(),
                "sets": set(),
                "taxa": set(),
            })

        size = xrefs.merge(processes=processes)
        logger.info("Disk usage: {:.0f}MB".format(size / 1024 ** 2))

        con, cur = dbms.connect(my_uri)
        cur.close()
        query = "UPDATE webfront_structure SET counts = %s WHERE accession = %s"
        with dbms.Populator(con, query) as table:
            for upid, _xrefs in xrefs:
                counts = reduce(_xrefs)
                counts["entries"]["total"] = sum(
                    counts["entries"].values())
                table.update((json.dumps(counts), upid))

        con.commit()
        con.close()

    logger.info("complete")


def get_protein2structures(my_uri: str) -> Dict[str, Set[str]]:
    protein2pdb = {}
    for pdb_id, s in get_structures(my_uri).items():
        for protein_ac in s["proteins"]:
            if protein_ac in protein2pdb:
                protein2pdb[protein_ac].add(pdb_id)
            else:
                protein2pdb[protein_ac] = {pdb_id}

    return protein2pdb
