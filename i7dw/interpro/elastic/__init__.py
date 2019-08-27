# -*- coding: utf-8 -*-

import os
import shutil
from multiprocessing import Queue
from tempfile import mkdtemp
from typing import List

from ... import io, logger
from .. import mysql
from . import utils


def init(path: str):
    try:
        shutil.rmtree(path)
    except FileNotFoundError:
        pass
    finally:
        os.makedirs(path)
        open(os.path.join(path, utils.LOADING_FILE), "w").close()


def write_documents(ora_ipr: str, my_ipr: str, src_proteins: str,
                    src_names: str, src_comments: str, src_proteomes: str,
                    src_matches: str, outdir: str, **kwargs):
    processes = kwargs.get("processes", 1)
    chunk_size = kwargs.get("chunk_size", 10000)
    limit = kwargs.get("limit", 0)

    processes = max(1, processes-1)  # minus one for parent process
    task_queue = Queue(processes)
    done_queue = Queue()
    workers = []
    for _ in range(processes):
        p = utils.DocumentProducer(ora_ipr, my_ipr, task_queue, done_queue,
                                   mkdtemp(dir=outdir))
        p.start()
        workers.append(p)

    # MySQL data
    logger.info("loading data from MySQL")
    taxa = mysql.taxonomy.get_taxa(my_ipr, lineage=True)
    integrated = {}
    entry_accessions = set()
    for entry_ac, e in mysql.entry.get_entries(my_ipr).items():
        entry_accessions.add(entry_ac)
        if e["integrated"]:
            integrated[entry_ac] = e["integrated"]

    # Open stores
    proteins = io.Store(src_proteins)
    protein2names = io.Store(src_names)
    protein2comments = io.Store(src_comments)
    protein2proteome = io.Store(src_proteomes)
    protein2matches = io.Store(src_matches)

    tax_ids = set(taxa.keys())

    logger.info("starting")
    n_proteins = 0
    chunk = []
    entries_with_matches = set()
    for acc, protein in proteins:
        tax_id = protein["taxon"]
        taxon = taxa[tax_id]

        name, other_names = protein2names.get(acc, (None, None))
        matches = protein2matches.get(acc, [])

        # Enqueue protein
        chunk.append((
            acc,
            protein["identifier"],
            name,
            "reviewed" if protein["is_reviewed"] else "unreviewed",
            protein["is_fragment"],
            protein["length"],
            protein2comments.get(acc, []),
            matches,
            protein2proteome.get(acc),
            taxon
        ))

        if len(chunk) == chunk_size:
            task_queue.put(("protein", chunk))
            chunk = []

        # Keep track of taxa associated to at least one protein
        try:
            tax_ids.remove(tax_id)
        except KeyError:
            pass

        # Keep track of entries with protein matches
        for m in matches:
            method_ac = m["method_ac"]
            entries_with_matches.add(method_ac)

            if method_ac in integrated:
                entries_with_matches.add(integrated[method_ac])

        n_proteins += 1
        if n_proteins == limit:
            break
        elif not n_proteins % 10000000:
            logger.info("{:>12,}".format(n_proteins))

    if chunk:
        task_queue.put(("protein", chunk))

    logger.info("{:>12,}".format(n_proteins))

    # Add entries without matches
    chunk = [
        (entry_ac,)
        for entry_ac in entry_accessions - entries_with_matches
    ]
    for i in range(0, len(chunk), chunk_size):
        task_queue.put(("entry", chunk[i:i+chunk_size]))

    # Add taxa without proteins
    chunk = [(taxa[tax_id],) for tax_id in tax_ids]
    for i in range(0, len(chunk), chunk_size):
        task_queue.put(("taxonomy", chunk[i:i+chunk_size]))

    # Poison pill
    for _ in workers:
        task_queue.put(None)

    # Closing stores
    proteins.close()
    protein2names.close()
    protein2comments.close()
    protein2proteome.close()
    protein2matches.close()

    n_docs = sum([done_queue.get() for _ in workers])

    # Wait for workers to finish
    for p in workers:
        p.join()

    # Delete loading file so Loaders know that all files are generated
    utils.set_ready(outdir)

    logger.info("complete: {:,} documents".format(n_docs))


def index_documents(uri: str, hosts: List[str], src: str, **kwargs):
    # Load databases (base of indices) - we exclude MobiDBlite
    indices = [utils.NODB_INDEX]
    for database in mysql.database.get_databases(uri).keys():
        if database != utils.MOBIDBLITE:
            indices.append(database)

    if kwargs.get("body_path"):
        # Create indices
        utils.create_indices(hosts, indices, kwargs.pop("body_path"), **kwargs)

    alias = kwargs.pop("alias", None)
    if utils.index_documents(hosts, src, **kwargs) and alias:
        utils.update_alias(hosts, indices, alias, **kwargs)


def update_alias(uri: str, hosts: List[str], **kwargs):
    # Load databases (base of indices) - we exclude MobiDBlite
    indices = [utils.NODB_INDEX]
    for database in mysql.database.get_databases(uri).keys():
        if database != utils.MOBIDBLITE:
            indices.append(database)

    utils.update_alias(hosts, indices, "current", **kwargs)
