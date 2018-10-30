#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
import os
import time
from multiprocessing import Process, Queue

from . import mysql
from .. import io

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s: %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)


def feed_store(filepath: str, queue: Queue, **kwargs: dict):
    with io.Store(filepath, **kwargs) as store:
        while True:
            chunk = queue.get()
            if chunk is None:
                break

            for args in chunk:
                store.update_from_seq(*args)

            store.flush()

        logging.info("temporary files ({}): {:,} bytes".format(
            os.path.basename(filepath), store.getsize()
        ))

        store.save()


def chunk_keys(keys: list, chunk_size: int) -> list:
    return [keys[i] for i in range(0, len(keys), chunk_size)]


def cross(my_uri: str, src_proteins: str, src_matches: str,
          src_proteomes: str, dst_entries: str, dst_proteomes: str,
          dst_structures: str, dst_taxa: str, chunk_size: int=10000,
          processes=4, tmpdir: str=None):

    logging.info("starting")

    """
    We do not keep track of set ---> entities 
    because these can be found from the union of the set's members
    """

    if tmpdir is not None:
        os.makedirs(tmpdir, exist_ok=True)

    # Get entries, initiate matches
    entry_matches = {}
    entry_database = {}
    integrated = {}
    for acc, e in mysql.get_entries(my_uri).items():
        entry_matches[acc] = 0
        entry_database[acc] = e["database"]

        if e["integrated"]:
            integrated[acc] = e["integrated"]

    entries_data = []
    entries_queue = Queue(maxsize=1)
    entries_proc = Process(target=feed_store,
                           args=(dst_entries, entries_queue),
                           kwargs={
                               "keys": chunk_keys(sorted(entry_matches), 10),
                               "tmpdir": tmpdir
                           })
    entries_proc.start()

    proteomes_data = []
    proteomes_queue = Queue(maxsize=1)
    proteomes_proc = Process(target=feed_store,
                             args=(dst_proteomes, proteomes_queue),
                             kwargs={
                                 "keys": chunk_keys(sorted(
                                     mysql.get_proteomes(my_uri)
                                 ), 10),
                                 "tmpdir": tmpdir
                             })
    proteomes_proc.start()

    taxa_data = []
    taxa_queue = Queue(maxsize=1)
    taxa_proc = Process(target=feed_store,
                        args=(dst_taxa, taxa_queue),
                        kwargs={
                            "keys": chunk_keys(sorted(
                                mysql.get_taxa(my_uri)
                            ), 10),
                            "tmpdir": tmpdir
                        })
    taxa_proc.start()

    # Set members
    entry_set = {
        entry_ac: set_ac
        for set_ac, s in mysql.get_sets(my_uri).items()
        for entry_ac in s["members"]
    }

    # Protein -> PDBe structure
    structure2others = {}
    protein2pdb = {}
    for pdb_id, s in mysql.get_structures(my_uri).items():
        for protein_ac in s["proteins"]:
            if protein_ac in protein2pdb:
                protein2pdb[protein_ac].add(pdb_id)
            else:
                protein2pdb[protein_ac] = {pdb_id}

        structure2others[pdb_id] = {
            "domains": set(),
            "entries": set(),
            "proteins": set(),
            "proteomes": set(),
            "sets": set(),
            "taxa": set()
        }

    # Get lineages so we can propagate relationships to a taxon's parents
    lineages = {
        # "lineage" stored as a string in MySQL (string include the taxon)
        tax_id: t["lineage"].strip().split()
        for tax_id, t in mysql.get_taxa(my_uri, lineage=True).items()
    }

    proteins = io.Store(src_proteins)
    protein2matches = io.Store(src_matches)
    protein2proteome = io.Store(src_proteomes)

    n_proteins = 0
    ts = time.time()
    for protein_ac, protein in proteins:
        tax_id = protein["taxon"]
        protein_id = protein["identifier"]
        matches = protein2matches.get(protein_ac, [])
        upid = protein2proteome.get(protein_ac)
        protein_structures = protein2pdb.get(protein_ac, [])

        # Create domain architecture, and count protein matches
        protein_entries = {}
        dom_entries = set()
        dom_arch = []
        for m in matches:
            method_ac = m["method_ac"]
            entry_ac = integrated.get(method_ac)

            if method_ac in protein_entries:
                protein_entries[method_ac] += 1
            else:
                protein_entries[method_ac] = 1

            if entry_ac:
                if entry_ac in protein_entries:
                    protein_entries[entry_ac] += 1
                else:
                    protein_entries[entry_ac] = 1

            if entry_database[method_ac] == "pfam":
                dom_entries.add(method_ac)
                if entry_ac:
                    dom_entries.add(entry_ac)
                    dom_arch.append("{}:{}".format(method_ac, entry_ac))
                else:
                    dom_arch.append("{}".format(method_ac))

        if dom_arch:
            dom_arch = '-'.join(dom_arch)

        # update matches and add domain architectures to entries
        protein_sets = set()
        for entry_ac, n_matches in protein_entries.items():
            entry_matches[entry_ac] += n_matches

            if entry_ac in dom_entries:
                # Has a domain architecture
                entries_data.append((entry_ac, "domains", dom_arch))

            if entry_ac in entry_set:
                protein_sets.add(entry_set[entry_ac])

        if upid:
            for entry_ac in protein_entries:
                has_domain = entry_ac in dom_entries
                database = entry_database[entry_ac]

                # Proteome <---> entries
                proteomes_data.append((upid, "entries", database, entry_ac))
                entries_data.append((entry_ac, "proteomes", upid))

                if has_domain:
                    # Proteome ---> IDA
                    proteomes_data.append((upid, "domains", dom_arch))

                # Entry ---> protein
                entries_data.append((entry_ac, "proteins",
                                     (protein_ac, protein_id)))

                # Entry <---> taxon
                entries_data.append((entry_ac, "taxa", tax_id))
                for _tax_id in lineages[tax_id]:
                    taxa_data.append((_tax_id, "entries", database, entry_ac))

                    if has_domain:
                        # Taxon ---> domain
                        taxa_data.append((_tax_id, "domains", dom_arch))

                # Entry <---> structure
                for pdb_id in protein_structures:
                    s = structure2others[pdb_id]
                    entries_data.append((entry_ac, "structures", pdb_id))
                    s["entries"].add((database, entry_ac))

                    if has_domain:
                        # Structure ---> domain
                        s["domains"].add(dom_arch)

            # Proteome <---> protein
            proteomes_data.append((upid, "proteins", protein_ac))

            # Proteome <---> taxon and taxon ---> protein
            proteomes_data.append((upid, "taxa", tax_id))
            for _tax_id in lineages[tax_id]:
                taxa_data.append((_tax_id, "proteomes", upid))
                taxa_data.append((_tax_id, "proteins", protein_ac))

            for set_ac in protein_sets:
                # Proteome ---> set
                proteomes_data.append((upid, "sets", set_ac))

                # Taxon --> set
                for _tax_id in lineages[tax_id]:
                    taxa_data.append((_tax_id, "sets", set_ac))

            for pdb_id in protein_structures:
                s = structure2others[pdb_id]

                # Proteome <---> structure
                proteomes_data.append((upid, "structures", pdb_id))
                s["proteomes"].add(upid)

                # Structure -> protein
                s["proteins"].add(protein_ac)

                # Structure ---> set
                for set_ac in protein_sets:
                    s["sets"].add(set_ac)

                # Structure <---> taxon
                s["taxa"].add(tax_id)
                for _tax_id in lineages[tax_id]:
                    taxa_data.append((_tax_id, "structures", pdb_id))
        else:
            for entry_ac in protein_entries:
                has_domain = entry_ac in dom_entries
                database = entry_database[entry_ac]

                # Entry ---> protein
                entries_data.append((entry_ac, "proteins",
                                     (protein_ac, protein_id)))

                # Entry <---> taxon
                entries_data.append((entry_ac, "taxa", tax_id))
                for _tax_id in lineages[tax_id]:
                    taxa_data.append((_tax_id, "entries", database, entry_ac))

                    if has_domain:
                        # Taxon ---> domain
                        taxa_data.append((_tax_id, "domains", dom_arch))

                # Entry <---> structure
                for pdb_id in protein_structures:
                    s = structure2others[pdb_id]
                    entries_data.append((entry_ac, "structures", pdb_id))
                    s["entries"].add((database, entry_ac))

                    if has_domain:
                        # Structure ---> domain
                        s["domains"].add(dom_arch)

            # taxon ---> protein
            for _tax_id in lineages[tax_id]:
                taxa_data.append((_tax_id, "proteins", protein_ac))

            for set_ac in protein_sets:
                # Taxon --> set
                for _tax_id in lineages[tax_id]:
                    taxa_data.append((_tax_id, "sets", set_ac))

            for pdb_id in protein_structures:
                s = structure2others[pdb_id]

                # Structure -> protein
                s["proteins"].add(protein_ac)

                # Structure ---> set
                for set_ac in protein_sets:
                    s["sets"].add(set_ac)

                # Structure <---> taxon
                s["taxa"].add(tax_id)
                for _tax_id in lineages[tax_id]:
                    taxa_data.append((_tax_id, "structures", pdb_id))

        n_proteins += 1
        if not n_proteins % chunk_size:
            entries_queue.put(entries_data)
            entries_data = []
            proteomes_queue.put(proteomes_data)
            proteomes_data = []
            taxa_queue.put(taxa_data)
            taxa_data = []

            logging.info('{:>12,} ({:.0f} proteins/sec)'.format(
                n_proteins, n_proteins / (time.time() - ts)
            ))

    proteins.close()
    protein2matches.close()
    protein2proteome.close()

    entries_queue.put(entries_data)
    entries_queue.put(None)
    proteomes_queue.put(proteomes_data)
    proteomes_queue.put(None)
    taxa_queue.put(taxa_data)
    taxa_queue.put(None)

    logging.info("merging")
    with io.Store(dst_entries) as store:
        store.reload()
        store.merge(processes=processes)

    with io.Store(dst_proteomes) as store:
        store.reload()
        store.merge(processes=processes)

    with io.Store(dst_taxa) as store:
        store.reload()
        store.merge(processes=processes)

    keys = chunk_keys(sorted(structure2others), 100)
    with io.Store(dst_structures, keys, tmpdir) as store:
        for pdb_id, v in structure2others.items():
            store.update(pdb_id, v)
        size = store.merge(processes=processes)

        logging.info("temporary files ({}): {:,} bytes".format(
            os.path.basename(dst_structures), size
        ))

    logging.info("complete")


def export(my_uri: str, src_proteins: str, src_matches: str,
           src_proteomes: str, dst_entries: str, dst_taxa: str,
           dst_proteomes: str, dst_sets: str, dst_structures: str,
           chunk_size: int=10000, tmpdir: str=None):

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

    entry_keys = chunk_keys(sorted(entries), chunk_size=10)
    taxa = chunk_keys(sorted(mysql.get_taxa(my_uri)), chunk_size=10)
    proteomes = chunk_keys(sorted(mysql.get_proteomes(my_uri)), chunk_size=10)

    sets = mysql.get_sets(my_uri)
    entry2set = {}
    for acc, s in sets.items():
        for entry_ac in s["members"]:
            entry2set[entry_ac] = acc
    sets = chunk_keys(sorted(sets), chunk_size=10)

    structures = mysql.get_structures(my_uri)
    protein2pdb = {}
    for pdb_id, s in structures.items():
        for acc in s["proteins"]:
            if acc in protein2pdb:
                protein2pdb[acc].add(pdb_id)
            else:
                protein2pdb[acc] = {pdb_id}
    structures = chunk_keys(sorted(structures), chunk_size=10)

    entries_chunk = []
    entries_queue = Queue(maxsize=1)
    entries_proc = Process(target=feed_store,
                           args=(dst_entries, entries_queue),
                           kwargs={"keys": entry_keys,
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

    taxa_chunk = []
    taxa_queue = Queue(maxsize=1)
    taxa_proc = Process(target=feed_store,
                        args=(dst_taxa, taxa_queue),
                        kwargs={"keys": taxa,
                                "tmpdir": tmpdir})

    for p in (entries_proc, taxa_proc, proteomes_proc, sets_proc,
              structures_proc):
        p.start()

    # Get lineages AFTER forking processes to reduce the memory usage
    lineages = {
        # Reverse the list and exclude the first element (previously last)
        # because it's the current taxon ID
        # -2: first element to include in the list is the second from the end
        # -1: negative step (reverse)
        tax_id: t["lineage"].strip().split()[-2::-1]
        for tax_id, t in mysql.get_taxa(my_uri, lineage=True).items()
    }

    proteins = io.Store(src_proteins)
    protein2matches = io.Store(src_matches)
    protein2proteome = io.Store(src_proteomes)

    n_proteins = 0
    ts = time.time()
    for acc, protein in proteins:
        tax_id = protein["taxon"]

        # Taxon ---> protein
        taxa_chunk.append((tax_id, "proteins", acc))
        for parent_id in lineages[tax_id]:
            taxa_chunk.append((parent_id, "proteins", acc))

        protein_entries = {}
        dom_entries = set()
        dom_arch = []
        for m in protein2matches.get(acc, []):
            method_ac = m["method_ac"]
            entry_ac = integrated.get(method_ac)

            if method_ac in protein_entries:
                protein_entries[method_ac] += 1
            else:
                protein_entries[method_ac] = 1

            if entry_ac:
                if entry_ac in protein_entries:
                    protein_entries[entry_ac] += 1
                else:
                    protein_entries[entry_ac] = 1

            if entries[method_ac]["database"] == "pfam":
                dom_entries.add(method_ac)
                if entry_ac:
                    dom_entries.add(entry_ac)
                    dom_arch.append("{}:{}".format(method_ac, entry_ac))
                else:
                    dom_arch.append("{}".format(method_ac))

        if dom_arch:
            dom_arch = '-'.join(dom_arch)

        for entry_ac in protein_entries:
            e = entries[entry_ac]
            e["matches"] += protein_entries[entry_ac]
            protein_entries[entry_ac] = e["database"]

            if entry_ac in dom_entries:
                # Entry ---> Domain architecture
                entries_chunk.append((entry_ac, "domains", dom_arch))

        upid = protein2proteome.get(acc)
        pdbe_ids = protein2pdb.get(acc, [])

        if upid:
            # Proteome <---> protein
            proteomes_chunk.append((upid, "proteins", acc))

            # Proteome <---> taxon
            proteomes_chunk.append((upid, "taxa", tax_id))
            taxa_chunk.append((tax_id, "proteomes", upid))
            for parent_id in lineages[tax_id]:
                taxa_chunk.append((parent_id, "proteomes", upid))

            for pdbe_id in pdbe_ids:
                # Structure ---> protein
                structures_chunk.append((pdbe_id, "proteins", acc))

                # Structure <---> taxon
                structures_chunk.append((pdbe_id, "taxa", tax_id))
                taxa_chunk.append((tax_id, "structures", pdbe_id))
                for parent_id in lineages[tax_id]:
                    taxa_chunk.append((parent_id, "structures", pdbe_id))

                # Structure <---> proteome
                structures_chunk.append((pdbe_id, "proteomes", upid))
                proteomes_chunk.append((upid, "structures", pdbe_id))

                _sets = set()
                for entry_ac, entry_db in protein_entries.items():
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
            for entry_ac, entry_db in protein_entries.items():
                # Entry ---> Protein
                entries_chunk.append(
                    (entry_ac, "proteins", (acc, protein["identifier"]))
                )

                # Entry <---> taxon
                entries_chunk.append((entry_ac, "taxa", tax_id))
                taxa_chunk.append((tax_id, "entries", entry_db, entry_ac))
                for parent_id in lineages[tax_id]:
                    taxa_chunk.append((parent_id, "entries", entry_db,
                                       entry_ac))

                # Proteome <---> entry
                proteomes_chunk.append((upid, "entries", entry_db, entry_ac))
                entries_chunk.append((entry_ac, "proteomes", upid))

                if entry_ac in entry2set:
                    set_ac = entry2set[entry_ac]
                    _sets.add(set_ac)

                    # Set <---> entry
                    sets_chunk.append((set_ac, "entries", entry_db, entry_ac))
                    entries_chunk.append((entry_ac, "sets", set_ac))

            for set_ac in _sets:
                # Set ---> protein
                sets_chunk.append((set_ac, "proteins", acc))

                # Proteome <---> set
                proteomes_chunk.append((upid, "sets", set_ac))
                sets_chunk.append((set_ac, "proteomes", upid))

                # Taxon <---> set
                taxa_chunk.append((tax_id, "sets", set_ac))
                sets_chunk.append((set_ac, "taxa", tax_id))
                for parent_id in lineages[tax_id]:
                    taxa_chunk.append((parent_id, "sets", set_ac))

        else:
            for pdbe_id in pdbe_ids:
                # Structure ---> protein
                structures_chunk.append((pdbe_id, "proteins", acc))

                # Structure <---> taxon
                structures_chunk.append((pdbe_id, "taxa", tax_id))
                taxa_chunk.append((tax_id, "structures", pdbe_id))
                for parent_id in lineages[tax_id]:
                    taxa_chunk.append((parent_id, "structures", pdbe_id))

                _sets = set()
                for entry_ac, entry_db in protein_entries.items():
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
            for entry_ac, entry_db in protein_entries.items():
                # Entry ---> Protein
                entries_chunk.append(
                    (entry_ac, "proteins", (acc, protein["identifier"]))
                )

                # Entry <---> taxon
                entries_chunk.append((entry_ac, "taxa", tax_id))
                taxa_chunk.append((tax_id, "entries", entry_db, entry_ac))
                for parent_id in lineages[tax_id]:
                    taxa_chunk.append((parent_id, "entries", entry_db,
                                       entry_ac))

                if entry_ac in entry2set:
                    set_ac = entry2set[entry_ac]
                    _sets.add(set_ac)

                    # Set <---> entry
                    sets_chunk.append((set_ac, "entries", entry_db, entry_ac))
                    entries_chunk.append((entry_ac, "sets", set_ac))

            for set_ac in _sets:
                # Set ---> protein
                sets_chunk.append((set_ac, "proteins", acc))

                # Taxon <---> set
                taxa_chunk.append((tax_id, "sets", set_ac))
                sets_chunk.append((set_ac, "taxa", tax_id))
                for parent_id in lineages[tax_id]:
                    taxa_chunk.append((parent_id, "sets", set_ac))

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

    # Adding match counts
    for entry_ac, e in entries.items():
        entries_chunk.append((entry_ac, "matches", e["matches"]))

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

    logging.info("merging")
    for f in (dst_entries, dst_taxa, dst_proteomes, dst_sets, dst_structures):
        with io.Store(f) as store:
            store.reload()
            store.merge(processes=6)

    logging.info("complete")
