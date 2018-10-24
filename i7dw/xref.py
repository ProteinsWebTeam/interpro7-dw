#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
import os
import time
from multiprocessing import Process, Queue

from . import disk, mysql

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

            for k, v in chunk:
                store.update(k, v)

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

    os.makedirs(tmpdir, exist_ok=True)

    entry2db = {
        acc: e["database"]
        for acc, e
        in mysql.get_entries(my_uri).items()
    }
    entries = chunk_keys(sorted(entry2db))
    taxa = chunk_keys(sorted(mysql.get_taxa(my_uri, method="basic")))
    proteomes = chunk_keys(sorted(mysql.get_proteomes(my_uri)))
    sets = chunk_keys(sorted(mysql.get_sets(my_uri)))
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
        taxa_chunk.append((tax_id, {"protein": {acc}}))

        _entries = {
            m["entry_ac"]
            for m in protein2matches.get(acc, [])
            if m["entry_ac"]
        }

        # Add source databases
        # TODO: remove `if in` check after debug
        _entries = [
            (entry_ac, entry2db[entry_ac])
            for entry_ac in _entries if entry_ac in entry2db
        ]

        for pdbe_id in protein2pdb.get(acc, []):
            structures_chunk.append((pdbe_id, {"protein": {acc}}))

            for entry_ac, entry_db in _entries:
                structures_chunk.append((pdbe_id,
                                         {"entry": {entry_db: entry_ac}}))
                entries_chunk.append((entry_ac,
                                      {"structure": pdbe_id}))

        upid = protein2proteome.get(acc)
        if upid:
            proteomes_chunk.append((upid, {"protein": {acc}}))

            for entry_ac, entry_db in _entries:
                proteomes_chunk.append((upid,
                                        {"entry": {entry_db: entry_ac}}))
                entries_chunk.append((entry_ac, {"proteome": upid}))

        for entry_ac, entry_db in _entries:
            entries_chunk.append((entry_ac, {"protein": {acc}}))
            entries_chunk.append((entry_ac, {"taxonomy": {tax_id}}))
            taxa_chunk.append((tax_id, {"entry": {entry_db: entry_ac}}))

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
            logging.info('{:>12} ({:.0f} proteins/sec)'.format(
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

    logging.info('{:>12} ({:.0f} proteins/sec)'.format(
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

    logging.info("complete")
