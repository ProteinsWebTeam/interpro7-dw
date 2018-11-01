#!/usr/bin/env python
# -*- coding: utf-8 -*-

import json
import logging
import os
import time
from multiprocessing import Process, Queue
from tempfile import mkstemp

from . import mysql
from .. import dbms, io

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s: %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)


def create_store(filepath: str, queue: Queue, **kwargs: dict):
    with io.Store(filepath, **kwargs) as store:
        while True:
            chunk = queue.get()
            if chunk is None:
                break

            for args in chunk:
                store.update_from_seq(*args)

            store.flush()

        store.save()


def chunk_keys(keys: list, chunk_size: int) -> list:
    return [keys[i] for i in range(0, len(keys), chunk_size)]


def update_disk(my_uri: str, src_proteins: str, src_matches: str,
                src_proteomes: str, dst_entries: str, flush: int=100000,
                processes=5, tmpdir: str=None):
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
    entries_proc = Process(target=create_store,
                           args=(dst_entries, entries_queue),
                           kwargs={
                               "keys": chunk_keys(sorted(entry_matches), 10),
                               "tmpdir": tmpdir
                           })
    entries_proc.start()

    fd, dst_proteomes = mkstemp(dir=tmpdir)
    os.close(fd)
    proteomes_data = []
    proteomes_queue = Queue(maxsize=1)
    proteomes_proc = Process(target=create_store,
                             args=(dst_proteomes, proteomes_queue),
                             kwargs={
                                 "keys": chunk_keys(sorted(
                                     mysql.get_proteomes(my_uri)), 100
                                 ),
                                 "tmpdir": tmpdir
                             })
    proteomes_proc.start()

    fd, dst_structures = mkstemp(dir=tmpdir)
    os.close(fd)
    structures_data = []
    structures_queue = Queue(maxsize=1)
    structures_proc = Process(target=create_store,
                              args=(dst_structures, structures_queue),
                              kwargs={
                                  "keys": chunk_keys(sorted(
                                      mysql.get_structures(my_uri)), 100
                                  ),
                                  "tmpdir": tmpdir
                              })
    structures_proc.start()

    fd, dst_taxa = mkstemp(dir=tmpdir)
    os.close(fd)
    taxa_data = []
    taxa_queue = Queue(maxsize=1)
    taxa_proc = Process(target=create_store,
                        args=(dst_taxa, taxa_queue),
                        kwargs={
                            "keys": chunk_keys(sorted(
                                mysql.get_taxa(my_uri)), 100
                            ),
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
    protein2pdb = {}
    for pdb_id, s in mysql.get_structures(my_uri).items():
        for protein_ac in s["proteins"]:
            if protein_ac in protein2pdb:
                protein2pdb[protein_ac].add(pdb_id)
            else:
                protein2pdb[protein_ac] = {pdb_id}

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
        else:
            dom_arch = None

        # update matches and add domain architectures to entries
        protein_sets = set()
        for entry_ac, n_matches in protein_entries.items():
            entry_matches[entry_ac] += n_matches

            if entry_ac in dom_entries:
                # Has a domain architecture
                entries_data.append((entry_ac, "domains", dom_arch))

            if entry_ac in entry_set:
                protein_sets.add(entry_set[entry_ac])

        # Taxon ---> protein
        taxa_data.append((tax_id, "proteins", protein_ac))

        if upid:
            # Proteome ---> protein
            proteomes_data.append((upid, "proteins", protein_ac))

            # Proteome <---> taxon
            proteomes_data.append((upid, "taxa", tax_id))
            taxa_data.append((tax_id, "proteomes", upid))

            # ---> Domain architecture
            if dom_arch:
                proteomes_data.append((upid, "domains", dom_arch))
                taxa_data.append((tax_id, "domains", dom_arch))

            for entry_ac in protein_entries:
                database = entry_database[entry_ac]

                # Proteome <---> entries
                proteomes_data.append((upid, "entries", database, entry_ac))
                entries_data.append((entry_ac, "proteomes", upid))

                # Entry ---> protein
                entries_data.append((entry_ac, "proteins", (protein_ac,
                                                            protein_id)))

                # Entry <---> taxon
                entries_data.append((entry_ac, "taxa", tax_id))
                taxa_data.append((tax_id, "entries", database, entry_ac))

                # Entry <---> structure
                for pdb_id in protein_structures:
                    entries_data.append((entry_ac, "structures", pdb_id))
                    structures_data.append((pdb_id, "entries", database,
                                            entry_ac))

            for set_ac in protein_sets:
                # Proteome ---> set
                proteomes_data.append((upid, "sets", set_ac))

                # Taxon --> set
                taxa_data.append((tax_id, "sets", set_ac))

            for pdb_id in protein_structures:
                # Structure ---> protein
                structures_data.append((pdb_id, "proteins", protein_ac))

                # Structure <---> taxon
                structures_data.append((pdb_id, "taxa", tax_id))
                taxa_data.append((tax_id, "structures", pdb_id))

                # Structure <---> proteome
                proteomes_data.append((upid, "structures", pdb_id))
                structures_data.append((pdb_id, "proteomes", upid))

                if dom_arch:
                    structures_data.append((pdb_id, "domains", dom_arch))

                # Structure ---> set
                for set_ac in protein_sets:
                    structures_data.append((pdb_id, "sets", set_ac))
        else:
            # ---> Domain architecture
            if dom_arch:
                taxa_data.append((tax_id, "domains", dom_arch))

            for entry_ac in protein_entries:
                database = entry_database[entry_ac]

                # Entry ---> protein
                entries_data.append((entry_ac, "proteins", (protein_ac,
                                                            protein_id)))

                # Entry <---> taxon
                entries_data.append((entry_ac, "taxa", tax_id))
                taxa_data.append((tax_id, "entries", database, entry_ac))

                # Entry <---> structure
                for pdb_id in protein_structures:
                    entries_data.append((entry_ac, "structures", pdb_id))
                    structures_data.append((pdb_id, "entries", database,
                                            entry_ac))

            for set_ac in protein_sets:
                # Taxon --> set
                taxa_data.append((tax_id, "sets", set_ac))

            for pdb_id in protein_structures:
                # Structure ---> protein
                structures_data.append((pdb_id, "proteins", protein_ac))

                # Structure <---> taxon
                structures_data.append((pdb_id, "taxa", tax_id))
                taxa_data.append((tax_id, "structures", pdb_id))

                if dom_arch:
                    structures_data.append((pdb_id, "domains", dom_arch))

                # Structure ---> set
                for set_ac in protein_sets:
                    structures_data.append((pdb_id, "sets", set_ac))

        n_proteins += 1
        if not n_proteins % flush:
            entries_queue.put(entries_data)
            entries_data = []
            proteomes_queue.put(proteomes_data)
            proteomes_data = []
            structures_queue.put(structures_data)
            structures_data = []
            taxa_queue.put(taxa_data)
            taxa_data = []

        if not n_proteins % 1000000:
            logging.info('{:>12,} ({:.0f} proteins/sec)'.format(
                n_proteins, n_proteins / (time.time() - ts)
            ))

    proteins.close()
    protein2matches.close()
    protein2proteome.close()

    entries_queue.put(entries_data)
    entries_queue.put(None)
    entries_data = None
    proteomes_queue.put(proteomes_data)
    proteomes_data = None
    proteomes_queue.put(None)
    structures_queue.put(structures_data)
    structures_data = None
    structures_queue.put(None)
    taxa_queue.put(taxa_data)
    taxa_data = None
    taxa_queue.put(None)

    logging.info('{:>12,} ({:.0f} proteins/sec)'.format(
        n_proteins, n_proteins / (time.time() - ts)
    ))

    entries_proc.join()
    proteomes_proc.join()
    structures_proc.join()
    taxa_proc.join()

    with io.Store(dst_entries) as store:
        store.reload()
        size = store.merge(processes=processes)
        logging.info("temporary files ({}): {:,} bytes".format(
            os.path.basename(dst_entries), size
        ))

    with io.Store(dst_proteomes) as store:
        store.reload()
        size = store.merge(processes=processes)
        logging.info("temporary files ({}): {:,} bytes".format(
            os.path.basename(dst_proteomes), size
        ))

    with io.Store(dst_structures) as store:
        store.reload()
        size = store.merge(processes=processes)
        logging.info("temporary files ({}): {:,} bytes".format(
            os.path.basename(dst_structures), size
        ))

    with io.Store(dst_taxa) as store:
        store.reload()
        size = store.merge(processes=processes)
        logging.info("temporary files ({}): {:,} bytes".format(
            os.path.basename(dst_taxa), size
        ))

    return


def update_mem(my_uri: str, src_proteins: str, src_matches: str,
               src_proteomes: str, dst_entries: str, flush: int=100000,
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
    entries_proc = Process(target=create_store,
                           args=(dst_entries, entries_queue),
                           kwargs={
                               "keys": chunk_keys(sorted(entry_matches), 10),
                               "tmpdir": tmpdir
                           })
    entries_proc.start()

    # Proteomes
    proteome2others = {
        upid: {
            "domains": set(),
            "entries": {},
            "proteins": 0,
            "sets": set(),
            "structures": set(),
            "taxa": set()
        } for upid in mysql.get_proteomes(my_uri)
    }

    # Set members
    entry_set = {
        entry_ac: set_ac
        for set_ac, s in mysql.get_sets(my_uri).items()
        for entry_ac in s["members"]
    }

    # Protein -> PDBe structure
    protein2pdb = {}
    structure2others = {}
    for pdb_id, s in mysql.get_structures(my_uri).items():
        for protein_ac in s["proteins"]:
            if protein_ac in protein2pdb:
                protein2pdb[protein_ac].add(pdb_id)
            else:
                protein2pdb[protein_ac] = {pdb_id}

        structure2others[pdb_id] = {
            "domains": set(),
            "entries": {},
            "proteins": 0,
            "proteomes": set(),
            "sets": set(),
            "taxa": set()
        }

    taxon2others = {
        tax_id: {
            "domains": set(),
            "entries": {},
            "proteins": 0,
            "proteins_total": 0,
            "proteomes": set(),
            "sets": set(),
            "structures": set()
        } for tax_id in mysql.get_taxa(my_uri)
    }

    proteins = io.Store(src_proteins)
    protein2matches = io.Store(src_matches)
    protein2proteome = io.Store(src_proteomes)

    n_proteins = 0
    ts = time.time()
    for protein_ac, protein in proteins:
        tax_id = protein["taxon"]
        try:
            taxon = taxon2others[tax_id]
        except KeyError:  # TODO: remove after debug
            continue

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
            proteome = proteome2others[upid]

            for entry_ac in protein_entries:
                has_domain = entry_ac in dom_entries
                database = entry_database[entry_ac]

                # Proteome <---> entries
                if database in proteome["entries"]:
                    proteome["entries"][database].add(entry_ac)
                else:
                    proteome["entries"][database] = {entry_ac}
                entries_data.append((entry_ac, "proteomes", upid))

                if has_domain:
                    # Proteome ---> IDA
                    proteome["domains"].add(dom_arch)

                # Entry ---> protein
                entries_data.append((entry_ac, "proteins", (protein_ac, protein_id)))

                # Entry <---> taxon
                entries_data.append((entry_ac, "taxa", tax_id))
                if database in taxon["entries"]:
                    taxon["entries"][database].add(entry_ac)
                else:
                    taxon["entries"][database] = {entry_ac}

                if has_domain:
                    # Taxon ---> domain
                    taxon["domains"].add(dom_arch)

                # Entry <---> structure
                for pdb_id in protein_structures:
                    s = structure2others[pdb_id]
                    if database in s["entries"]:
                        s["entries"][database].add(entry_ac)
                    else:
                        s["entries"][database] = {entry_ac}
                    entries_data.append((entry_ac, "structures", pdb_id))

                    if has_domain:
                        # Structure ---> domain
                        s["domains"].add(dom_arch)

            # Proteome <---> protein
            proteome["proteins"] += 1

            # Proteome <---> taxon and taxon ---> protein
            proteome["taxa"].add(tax_id)
            taxon["proteomes"].add(upid)
            taxon["proteins"] += 1

            for set_ac in protein_sets:
                # Proteome ---> set
                proteome["sets"].add(set_ac)

                # Taxon --> set
                taxon["sets"].add(set_ac)

            for pdb_id in protein_structures:
                s = structure2others[pdb_id]

                # Proteome <---> structure
                proteome["structures"].add(pdb_id)
                s["proteomes"].add(upid)

                # Structure -> protein
                s["proteins"] += 1

                # Structure ---> set
                for set_ac in protein_sets:
                    s["sets"].add(set_ac)

                # Structure <---> taxon
                s["taxa"].add(tax_id)
                taxon["structures"].add(pdb_id)
        else:
            for entry_ac in protein_entries:
                has_domain = entry_ac in dom_entries
                database = entry_database[entry_ac]

                # Entry ---> protein
                entries_data.append((entry_ac, "proteins", (protein_ac, protein_id)))

                # Entry <---> taxon
                entries_data.append((entry_ac, "taxa", tax_id))
                if database in taxon["entries"]:
                    taxon["entries"][database].add(entry_ac)
                else:
                    taxon["entries"][database] = {entry_ac}

                if has_domain:
                    # Taxon ---> domain
                    taxon["domains"].add(dom_arch)

                # Entry <---> structure
                for pdb_id in protein_structures:
                    s = structure2others[pdb_id]
                    if database in s["entries"]:
                        s["entries"][database].add(entry_ac)
                    else:
                        s["entries"][database] = {entry_ac}
                    entries_data.append((entry_ac, "structures", pdb_id))

                    if has_domain:
                        # Structure ---> domain
                        s["domains"].add(dom_arch)

            # taxon ---> protein
            taxon["proteins"] += 1

            for set_ac in protein_sets:
                # Taxon --> set
                taxon["sets"].add(set_ac)

            for pdb_id in protein_structures:
                s = structure2others[pdb_id]

                # Structure -> protein
                s["proteins"] += 1

                # Structure ---> set
                for set_ac in protein_sets:
                    s["sets"].add(set_ac)

                # Structure <---> taxon
                s["taxa"].add(tax_id)
                taxon["structures"].add(pdb_id)

        n_proteins += 1
        if not n_proteins % flush:
            entries_queue.put(entries_data)
            entries_data = []

        if not n_proteins % 1000000:
            logging.info('{:>12,} ({:.0f} proteins/sec)'.format(
                n_proteins, n_proteins / (time.time() - ts)
            ))

    proteins.close()
    protein2matches.close()
    protein2proteome.close()

    entries_queue.put(entries_data)
    entries_queue.put(None)

    logging.info('{:>12,} ({:.0f} proteins/sec)'.format(
        n_proteins, n_proteins / (time.time() - ts)
    ))

    logging.info("updating tables")
    con, cur = dbms.connect(my_uri)
    for upid, xrefs in proteome2others.items():
        counts = reduce(xrefs)
        counts["entries"]["total"] = sum(counts["entries"].values())

        cur.execute(
            """
            UPDATE webfront_proteome
            SET counts = %s
            WHERE accession = %s
            """,
            (json.dumps(counts), upid)
        )
    proteome2others = None

    for pdb_id, xrefs in structure2others.items():
        counts = reduce(xrefs)
        counts["entries"]["total"] = sum(counts["entries"].values())

        cur.execute(
            """
            UPDATE webfront_structure
            SET counts = %s
            WHERE accession = %s
            """,
            (json.dumps(counts), pdb_id)
        )
    structure2others = None

    for tax_id, t in mysql.get_taxa(my_uri, lineage=True).items():
        # "lineage" stored as a string in MySQL (string include the taxon)
        lineage = t["lineage"].strip().split()[-2::-1]

        # propagate relationships to parents
        t = taxon2others[tax_id]
        for parent_id in lineage:
            p = taxon2others[parent_id]
            p["domains"] |= t["domains"]
            p["proteomes"] |= t["proteomes"]
            p["sets"] |= t["sets"]
            p["structures"] |= t["structures"]

            for db, db_entries in t["entries"].items():
                if db in p["entries"]:
                    p["entries"][db] |= db_entries
                else:
                    p["entries"][db] = db_entries

            p["proteins_total"] += t["proteins"]

        del t["proteins"]

    for tax_id, xrefs in taxon2others.items():
        counts = reduce(xrefs)
        counts["proteins"] = counts.pop("proteins_total")
        counts["entries"]["total"] = sum(counts["entries"].values())

        cur.execute(
            """
            UPDATE webfront_taxonomy
            SET counts = %s
            WHERE accession = %s
            """,
            (json.dumps(counts), tax_id)
        )
    taxon2others = None

    entries_proc.join()

    with io.Store(dst_entries) as store:
        store.merge(processes=processes)

        for entry_ac, xrefs in store.items(processes):
            counts = reduce(xrefs)
            counts["matches"] = entry_matches.pop(entry_ac)

            cur.execute(
                """
                UPDATE webfront_entry
                SET counts = %s
                WHERE accession = %s
                """,
                (json.dumps(counts), entry_ac)
            )

    for entry_ac in entry_matches:
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
                "sets": 0,
                "structures": 0,
                "taxa": 0
            }), entry_ac)
        )

    con.commit()
    cur.close()
    con.close()
    logging.info("complete")


def reduce(src: dict):
    dst = {}
    for k, v in src.items():
        if isinstance(v, dict):
            dst[k] = reduce(v)
        elif isinstance(v, set):
            dst[k] = len(v)
        else:
            dst[k] = v

    return dst
