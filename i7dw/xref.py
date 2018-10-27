#!/usr/bin/env python
# -*- coding: utf-8 -*-

import json
import logging
import os
import time
from multiprocessing import Process, Queue

from . import dbms, disk, mysql

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s: %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)


def feed_store(filepath: str, queue: Queue, **kwargs: dict):
    with disk.Store(filepath, **kwargs) as store:
        while True:
            chunk = queue.get()
            if chunk is None:
                break

            for args in chunk:
                store.update_from_seq(*args)

            store.flush()

        logging.info("temporary files ({}): {:,} bytes)".format(
            os.path.basename(filepath), store.getsize()
        ))

        store.save()


def chunk_keys(keys: list, chunk_size: int=1000) -> list:
    return [keys[i] for i in range(0, len(keys), chunk_size)]


def count_xrefs(my_uri, src_proteins, src_matches, src_proteomes,
                dst_entries, dst_taxa, dst_proteomes, dst_sets,
                dst_structures, chunk_size=100000, tmpdir=None):

    logging.info("starting")

    if tmpdir is not None:
        os.makedirs(tmpdir, exist_ok=True)

    entries = {}
    integrated = {}
    for acc, e in mysql.get_entries(my_uri).items():
        entries[acc] = {
            "database": e["database"],
            "matches": 0
        }

        if e["integrated"]:
            integrated[acc] = e["integrated"]

    entry_keys = chunk_keys(sorted(entries), chunk_size=100)
    taxa = chunk_keys(sorted(mysql.get_taxa(my_uri, method="basic")))
    proteomes = chunk_keys(sorted(mysql.get_proteomes(my_uri)))

    sets = mysql.get_sets(my_uri)
    entry2set = {}
    for acc, s in sets.items():
        for entry_ac in s["members"]:
            entry2set[entry_ac] = acc
    sets = chunk_keys(sorted(sets))

    structures = mysql.get_structures(my_uri)
    protein2pdb = {}
    for pdb_id, s in structures.items():
        for acc in s["proteins"]:
            if acc in protein2pdb:
                protein2pdb[acc].add(pdb_id)
            else:
                protein2pdb[acc] = {pdb_id}
    structures = chunk_keys(sorted(structures))

    entries_chunk = []
    entries_queue = Queue(maxsize=1)
    entries_proc = Process(target=feed_store,
                           args=(dst_entries, entries_queue),
                           kwargs={"keys": entry_keys,
                                   "tmpdir": tmpdir})

    taxa_chunk = []
    taxa_queue = Queue(maxsize=1)
    taxa_proc = Process(target=feed_store,
                        args=(dst_taxa, taxa_queue),
                        kwargs={"keys": taxa,
                                "tmpdir": tmpdir})

    proteomes_chunk = []
    proteomes_queue = Queue(maxsize=1)
    proteomes_proc = Process(target=feed_store,
                             args=(dst_proteomes, proteomes_queue),
                             kwargs={"keys": proteomes,
                                     "tmpdir": tmpdir})

    sets_chunk = []
    sets_queue = Queue(maxsize=1)
    sets_proc = Process(target=feed_store,
                        args=(dst_sets, sets_queue),
                        kwargs={"keys": sets,
                                "tmpdir": tmpdir})

    structures_chunk = []
    structures_queue = Queue(maxsize=1)
    structures_proc = Process(target=feed_store,
                              args=(dst_structures, structures_queue),
                              kwargs={"keys": structures,
                                      "tmpdir": tmpdir})

    for p in (entries_proc, taxa_proc, proteomes_proc, sets_proc,
              structures_proc):
        p.start()

    proteins = disk.Store(src_proteins)
    protein2matches = disk.Store(src_matches)
    protein2proteome = disk.Store(src_proteomes)

    n_proteins = 0
    ts = time.time()
    for acc, protein in proteins:
        tax_id = protein["taxon"]
        # Taxon ---> protein
        taxa_chunk.append((tax_id, "proteins", acc))

        _entries = {}
        for m in protein2matches.get(acc, []):
            method_ac = m["method_ac"]
            if method_ac in _entries:
                _entries[method_ac] += 1
            else:
                _entries[method_ac] = 1

            if method_ac in integrated:
                entry_ac = integrated[method_ac]
                if entry_ac in _entries:
                    _entries[entry_ac] += 1
                else:
                    _entries[entry_ac] = 1

        for entry_ac in _entries:
            try:
                e = entries[entry_ac]
            except KeyError:
                # TODO: remove try/except after debug
                db = None
            else:
                e["matches"] += _entries[entry_ac]
                db = e["database"]
            finally:
                _entries[entry_ac] = db

        upid = protein2proteome.get(acc)
        pdbe_ids = protein2pdb.get(acc, [])

        if upid:
            # Proteome <---> protein
            proteomes_chunk.append((upid, "proteins", acc))

            # Proteome <---> taxon
            proteomes_chunk.append((upid, "taxa", tax_id))
            taxa_chunk.append((tax_id, "proteomes", upid))

            for pdbe_id in pdbe_ids:
                # Structure ---> protein
                structures_chunk.append((pdbe_id, "proteins", acc))

                # Structure <---> taxon
                structures_chunk.append((pdbe_id, "taxa", tax_id))
                taxa_chunk.append((tax_id, "structures", pdbe_id))

                # Structure <---> proteome
                structures_chunk.append((pdbe_id, "proteomes", upid))
                proteomes_chunk.append((upid, "structures", pdbe_id))

                _sets = set()
                for entry_ac, entry_db in _entries.items():
                    if not entry_db:
                        continue  # TODO: remove after debug

                    # Structure <---> entry
                    structures_chunk.append(
                        (pdbe_id, "entries", entry_db, entry_ac)
                    )
                    entries_chunk.append((entry_ac, "structures", pdbe_id))

                    if entry_ac in entry2set:
                        _sets.add(entry2set[entry_ac])

                for set_ac in _sets:
                    # Structure <---> set
                    structures_chunk.append((pdbe_id, "sets", set_ac))
                    sets_chunk.append((set_ac, "structures", pdbe_id))

            _sets = set()
            for entry_ac, entry_db in _entries.items():
                if not entry_db:
                    continue  # TODO: remove after debug

                # Entry ---> Protein
                entries_chunk.append(
                    (entry_ac, "proteins", (acc, protein["identifier"]))
                )

                # Entry <---> taxon
                entries_chunk.append((entry_ac, "taxa", tax_id))
                taxa_chunk.append((tax_id, "entries", entry_db, entry_ac))

                # Proteome <---> entry
                proteomes_chunk.append((upid, "entries", entry_db, entry_ac))
                entries_chunk.append((entry_ac, "proteomes", upid))

                if entry_ac in entry2set:
                    _sets.add(entry2set[entry_ac])

                for set_ac in _sets:
                    # Set ---> protein
                    sets_chunk.append((set_ac, "proteins", acc))

                    # Proteome <---> set
                    proteomes_chunk.append((upid, "sets", set_ac))
                    sets_chunk.append((set_ac, "proteomes", upid))

                    # Taxon <---> set
                    taxa_chunk.append((tax_id, "sets", set_ac))
                    sets_chunk.append((set_ac, "taxa", tax_id))
        else:
            for pdbe_id in pdbe_ids:
                # Structure ---> protein
                structures_chunk.append((pdbe_id, "proteins", acc))

                # Structure <---> taxon
                structures_chunk.append((pdbe_id, "taxa", tax_id))
                taxa_chunk.append((tax_id, "structures", pdbe_id))

                _sets = set()
                for entry_ac, entry_db in _entries.items():
                    if not entry_db:
                        continue  # TODO: remove after debug

                    # Structure <---> entry
                    structures_chunk.append(
                        (pdbe_id, "entries", entry_db, entry_ac)
                    )
                    entries_chunk.append((entry_ac, "structures", pdbe_id))

                    if entry_ac in entry2set:
                        _sets.add(entry2set[entry_ac])

                for set_ac in _sets:
                    # Structure <---> set
                    structures_chunk.append((pdbe_id, "sets", set_ac))
                    sets_chunk.append((set_ac, "structures", pdbe_id))

            _sets = set()
            for entry_ac, entry_db in _entries.items():
                if not entry_db:
                    continue  # TODO: remove after debug

                # Entry ---> Protein
                entries_chunk.append(
                    (entry_ac, "proteins", (acc, protein["identifier"]))
                )

                # Entry <---> taxon
                entries_chunk.append((entry_ac, "taxa", tax_id))
                taxa_chunk.append((tax_id, "entries", entry_db, entry_ac))

                if entry_ac in entry2set:
                    _sets.add(entry2set[entry_ac])

                for set_ac in _sets:
                    # Set ---> protein
                    sets_chunk.append((set_ac, "proteins", acc))

                    # Taxon <---> set
                    taxa_chunk.append((tax_id, "sets", set_ac))
                    sets_chunk.append((set_ac, "taxa", tax_id))

        n_proteins += 1
        if not n_proteins % chunk_size:
            entries_queue.put(entries_chunk)
            entries_chunk = []
            taxa_queue.put(taxa_chunk)
            taxa_chunk = []
            proteomes_queue.put(proteomes_chunk)
            proteomes_chunk = []
            sets_queue.put(sets_chunk)
            sets_chunk = []
            structures_queue.put(structures_chunk)
            structures_chunk = []

        if not n_proteins % 1000000:
            logging.info('{:>12,} ({:.0f} proteins/sec)'.format(
                n_proteins, n_proteins / (time.time() - ts)
            ))

    if entries_chunk:
        entries_queue.put(entries_chunk)
    if taxa_chunk:
        taxa_queue.put(taxa_chunk)
    if proteomes_chunk:
        proteomes_queue.put(proteomes_chunk)
    if sets_chunk:
        sets_queue.put(sets_chunk)
    if structures_chunk:
        structures_queue.put(structures_chunk)

    logging.info('{:>12,} ({:.0f} proteins/sec)'.format(
        n_proteins, n_proteins / (time.time() - ts)
    ))

    proteins.close()
    protein2matches.close()
    protein2proteome.close()

    for q in (entries_queue, taxa_queue, proteomes_queue, sets_queue,
              structures_queue):
        q.put(None)

    for p in (entries_proc, taxa_proc, proteomes_proc, sets_proc,
              structures_proc):
        p.join()

    for f in (dst_entries, dst_taxa, dst_proteomes, dst_sets, dst_structures):
        logging.info("{}: merging".format(os.path.basename(f)))

        with disk.Store(f) as store:
            store.reload()
            store.merge(processes=6)

        logging.info("{}: merged".format(os.path.basename(f)))

    logging.info("updating tables")

    con, cur = dbms.connect(my_uri)
    with disk.Store(dst_entries) as store:
        for entry_ac, data in store:
            counts = aggregate(data)
            counts["sets"] = 1 if entry_ac in entry2set else 0
            counts["matches"] = entries[entry_ac]["matches"]

            cur.execute(
                """
                UPDATE webfront_entry
                SET counts = %s
                WHERE accession = %s
                """,
                (json.dumps(counts), entry_ac)
            )

    with disk.Store(dst_taxa) as store:
        for tax_id, data in store:
            cur.execute(
                """
                UPDATE webfront_taxonomy
                SET counts = %s
                WHERE accession = %s
                """,
                (json.dumps(aggregate(data)), tax_id)
            )

    with disk.Store(dst_proteomes) as store:
        for upid, data in store:
            cur.execute(
                """
                UPDATE webfront_proteome
                SET counts = %s
                WHERE accession = %s
                """,
                (json.dumps(aggregate(data)), upid)
            )

    with disk.Store(dst_sets) as store:
        sets = mysql.get_sets(my_uri)
        for set_ac, data in store:
            counts = aggregate(data)
            counts["entries"] = len(sets[set_ac]["members"])
            cur.execute(
                """
                UPDATE webfront_set
                SET counts = %s
                WHERE accession = %s
                """,
                (json.dumps(counts), set_ac)
            )

    with disk.Store(dst_structures) as store:
        for pdbe_id, data in store:
            cur.execute(
                """
                UPDATE webfront_structure
                SET counts = %s
                WHERE accession = %s
                """,
                (json.dumps(aggregate(data)), pdbe_id)
            )

    con.commit()
    cur.close()
    con.close()

    logging.info("complete")


def aggregate(src: dict):
    dst = {}
    for k, v in src.items():
        if isinstance(v, dict):
            dst[k] = aggregate(v)
        else:
            dst[k] = len(v)

    return dst
