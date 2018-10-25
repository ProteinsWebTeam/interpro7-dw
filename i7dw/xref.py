#!/usr/bin/env python
# -*- coding: utf-8 -*-

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

        logging.info("{} filled".format(filepath))
        size = store.merge()
        logging.info("{}: temporary files: {} bytes)".format(
            filepath, size
        ))


def chunk_keys(keys: list, chunk_size: int=1000) -> list:
    return [keys[i] for i in range(0, len(keys), chunk_size)]


def count_xrefs(my_uri, src_proteins, src_matches, src_proteomes,
                dst_entries, dst_taxa, dst_proteomes, dst_sets,
                dst_structures, chunk_size=100000, tmpdir=None):

    logging.info("starting")

    if tmpdir is not None:
        os.makedirs(tmpdir, exist_ok=True)

    entry2db = {
        acc: e["database"]
        for acc, e
        in mysql.get_entries(my_uri).items()
    }
    entries = chunk_keys(sorted(entry2db))
    taxa = chunk_keys(sorted(mysql.get_taxa(my_uri, method="basic")))
    proteomes = chunk_keys(sorted(mysql.get_proteomes(my_uri)))

    sets = mysql.get_sets(my_uri)
    entry2set = {}
    for acc, s in sets.items():
        for entry_ac in s["members"]:
            entry2set[entry_ac] = acc
    sets = sorted(sets)

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
                           kwargs={"keys": entries,
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

        _entries = set()
        for m in protein2matches.get(acc, []):
            _entries.add(m["method_ac"])

            if m["entry_ac"]:
                _entries.add(m["entry_ac"])

        # Add source databases
        # TODO: remove `if in` check after debug
        _entries = [
            (entry_ac, entry2db[entry_ac])
            for entry_ac in _entries if entry_ac in entry2db
        ]

        upid = protein2proteome.get(acc)
        pdbe_ids = protein2pdb.get(acc, [])

        if upid:
            # Proteome ---> protein
            proteomes_chunk.append((upid, "proteins", acc))

            # Proteome <---> taxon
            proteomes_chunk.append((upid, "taxa", tax_id))
            taxa_chunk.append((tax_id, "proteomes", upid))

            for pdbe_id in pdbe_ids:
                # # Structure ---> protein
                # structures_chunk.append((pdbe_id, {"proteins": {acc}}))

                # Structure <---> taxon
                structures_chunk.append((pdbe_id, "taxa", tax_id))
                taxa_chunk.append((tax_id, "structures", pdbe_id))

                # Structure <---> proteome
                structures_chunk.append((pdbe_id, "proteomes", upid))
                # proteomes_chunk.append((upid, {"structures": {pdbe_id}}))

                _sets = set()
                for entry_ac, entry_db in _entries:
                    # Structure <---> entry
                    structures_chunk.append(
                        (pdbe_id, "entries", entry_db, entry_ac)
                    )
                    # entries_chunk.append(
                    #     (entry_ac, {"structures": {pdbe_id}})
                    # )

                    # if entry_ac in entry2set:
                    #     _sets.add(entry2set[entry_ac])

                # for set_ac in _sets:
                #     # Structure <---> set
                #     structures_chunk.append((pdbe_id, {"sets": {set_ac}}))
                #     sets_chunk.append((set_ac, {"structures": {pdbe_id}}))

            _sets = set()
            for entry_ac, entry_db in _entries:
                # Entry ---> Protein
                entries_chunk.append((entry_ac, "proteins", acc))

                # Entry <---> taxon
                entries_chunk.append((entry_ac, "taxa", tax_id))
                taxa_chunk.append(
                    (tax_id, "entries", entry_db, entry_ac)
                )

                # Proteome <---> entry
                proteomes_chunk.append(
                    (upid, "entries", entry_db, entry_ac)
                )
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
                # # Structure ---> protein
                # structures_chunk.append((pdbe_id, {"proteins": {acc}}))

                # Structure <---> taxon
                structures_chunk.append((pdbe_id, "taxa", tax_id))
                taxa_chunk.append((tax_id, "structures", pdbe_id))

                _sets = set()
                for entry_ac, entry_db in _entries:
                    # Structure <---> entry
                    structures_chunk.append(
                        (pdbe_id, "entries", entry_db, entry_ac)
                    )
                    # entries_chunk.append(
                    #     (entry_ac, {"structures": {pdbe_id}})
                    # )

                    # if entry_ac in entry2set:
                    #     _sets.add(entry2set[entry_ac])

                # for set_ac in _sets:
                #     # Structure <---> set
                #     structures_chunk.append((pdbe_id, {"sets": {set_ac}}))
                #     sets_chunk.append((set_ac, {"structures": {pdbe_id}}))

            _sets = set()
            for entry_ac, entry_db in _entries:
                # Entry ---> Protein
                entries_chunk.append((entry_ac, "proteins", acc))

                # Entry <---> taxon
                entries_chunk.append((entry_ac, "taxa", tax_id))
                taxa_chunk.append(
                    (tax_id, "entries", entry_db, entry_ac)
                )

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

    # #con, cur = dbms.connect(my_uri)
    # with disk.Store(dst_entries) as store:
    #     for acc, data in store:
    #         print(acc, aggregate(data))
    #         break
    #
    #         if acc in entry2set:
    #             pass  # has a set! so count = 1

    logging.info("complete")


def aggregate(src: dict):
    dst = {}
    for k, v in src.items():
        if isinstance(v, dict):
            dst[k] = aggregate(v)
        else:
            if k == "taxonomy":
                k = "taxa"
            elif k == "entry":
                k = "entries"
            else:
                k += "s"
            dst[k] = len(v)

    return dst
