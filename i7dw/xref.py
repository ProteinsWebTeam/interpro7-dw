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
    with disk.KVStore(filepath, **kwargs) as store:
        while True:
            chunks = queue.get()
            if chunks is None:
                break

            for chunk in chunks:
                store.add(*chunk)

            store.flush()

        logging.info("{} filled".format(filepath))

        tmp_usage = store.close()

    logging.info("{} ready (temporary files: {} bytes)".format(
        filepath, tmp_usage
    ))


def count_xrefs(my_uri, proteins_f, prot_matches_f, proteomes_f,
                entries_kvf, taxa_kvf, proteomes_kvf, sets_kvf,
                structures_kvf, chunk_size=100000, compress=False,
                tmpdir=None):

    try:
        os.makedirs(tmpdir, exist_ok=True)
    except (FileNotFoundError, TypeError):
        tmpdir = None

    entry2db = {
        acc: e["database"]
        for acc, e
        in mysql.get_entries(my_uri).items()
    }
    entries = sorted(entry2db)
    taxa = sorted(mysql.get_taxa(my_uri, slim=True))
    proteomes = sorted(mysql.get_proteomes(my_uri))
    sets = sorted(mysql.get_sets(my_uri, by_members=False))
    structures = mysql.get_structures(my_uri)

    prot2struct = {}
    for s in structures.values():
        for p in s["proteins"]:
            if p in prot2struct:
                prot2struct[p].add(s["accession"])
            else:
                prot2struct[p] = {s["accession"]}

    structures = sorted(structures)

    l_entries = []
    q_entries = Queue(maxsize=1)
    p_entries = Process(target=feed_store,
                        args=(entries_kvf, q_entries),
                        kwargs={"compress": compress,
                                "ids": entries,
                                "tmpdir": tmpdir})

    l_taxa = []
    q_taxa = Queue(maxsize=1)
    p_taxa = Process(target=feed_store,
                     args=(taxa_kvf, q_taxa),
                     kwargs={"compress": compress,
                             "ids": taxa,
                             "tmpdir": tmpdir})

    l_proteomes = []
    q_proteomes = Queue(maxsize=1)
    p_proteomes = Process(target=feed_store,
                          args=(proteomes_kvf, q_proteomes),
                          kwargs={"compress": compress,
                                  "ids": proteomes,
                                  "tmpdir": tmpdir})

    l_sets = []
    q_sets = Queue(maxsize=1)
    p_sets = Process(target=feed_store,
                     args=(sets_kvf, q_sets),
                     kwargs={"compress": compress,
                             "ids": sets,
                             "tmpdir": tmpdir})

    l_structures = []
    q_structures = Queue(maxsize=1)
    p_structures = Process(target=feed_store,
                           args=(structures_kvf, q_structures),
                           kwargs={"compress": compress,
                                   "ids": structures,
                                   "tmpdir": tmpdir})

    for p in (p_entries, p_taxa, p_proteomes, p_sets, p_structures):
        p.start()

    proteins_s = disk.Store(proteins_f)
    prot_matches_s = disk.Store(prot_matches_f)
    proteomes_s = disk.Store(proteomes_f)

    n_proteins = 0
    ts = time.time()
    for acc, protein in proteins_s.iter():
        tax_id = protein["taxon"]
        matches = prot_matches_s.get(acc, [])
        proteomes = proteomes_s.get(acc, [])
        structures = prot2struct.get(acc, [])

        l_taxa.append((tax_id, "protein", acc))

        _entries = set()
        for m in matches:
            _entries.add(m["method_ac"])
            if m["entry_ac"]:
                _entries.add(m["entry_ac"])

        # Add source databases
        _entries = [
            (entry_ac, entry2db[entry_ac])
            for entry_ac in _entries
        ]

        for pdbe_id in structures:
            l_structures.append((pdbe_id, "protein", acc))

            for entry_ac, entry_db in _entries:
                l_structures.append((pdbe_id, "entry", entry_db, entry_ac))
                l_entries.append((entry_ac, "structure", pdbe_id))

        for upid in proteomes:
            l_proteomes.append((upid, "protein", acc))

            for entry_ac, entry_db in _entries:
                l_proteomes.append((upid, "entry", entry_db, entry_ac))
                l_entries.append((entry_ac, "proteome", upid))

        for entry_ac, entry_db in _entries:
            l_entries.append((entry_ac, "protein", acc))
            l_entries.append((entry_ac, "taxon", tax_id))
            l_taxa.append((tax_id, "entry", entry_db, entry_ac))

        n_proteins += 1

        if not n_proteins % chunk_size:
            q_entries.put(l_entries)
            l_entries = []
            q_taxa.put(l_taxa)
            l_taxa = []
            q_proteomes.put(l_proteomes)
            l_proteomes = []
            q_sets.put(l_sets)
            l_sets = []
            q_structures.put(l_structures)
            l_structures = []

        if not n_proteins % 1000000:
            logging.info('{:>12} ({:.0f} proteins/sec)'.format(
                n_proteins, n_proteins // (time.time() - ts)
            ))

    logging.info('{:>12} ({:.0f} proteins/sec)'.format(
        n_proteins, n_proteins // (time.time() - ts)
    ))

    proteins_s.close()
    prot_matches_s.close()
    proteomes_s.close()

    for q in (q_entries, q_taxa, q_proteomes,  q_sets, q_structures):
        q.put(None)
        
    for p in (p_entries, p_taxa, p_proteomes,  p_sets, p_structures):
        p.join()

    logging.info("complete")
