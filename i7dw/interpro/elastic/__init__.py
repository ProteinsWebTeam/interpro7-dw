# -*- coding: utf-8 -*-

import os
import shutil
from multiprocessing import Queue
from tempfile import mkdtemp
from typing import List

from i7dw import io, logger
from i7dw.interpro import mysql
from . import utils


def init(path: str):
    try:
        shutil.rmtree(path)
    except FileNotFoundError:
        pass
    finally:
        os.makedirs(path)
        open(os.path.join(path, utils.LOADING_FILE), "w").close()


def write_documents(url: str, src_comments: str, src_matches: str,
                    src_names: str, src_proteins: str, src_proteomes: str,
                    outdir: str, chunk_size: int=10000, processes: int=1):
    task_queue = Queue(processes)
    done_queue = Queue()
    workers = []
    for _ in range(max(1, processes-1)):
        _outdir = mkdtemp(dir=outdir)
        p = utils.DocumentProducer(url, task_queue, done_queue, _outdir)
        p.start()
        workers.append(p)

    # MySQL data
    logger.info("loading data")
    taxa = {t["id"]: t for t in mysql.taxonomy.iter_taxa(url, lineage=True)}

    # Open stores
    comments = io.Store(src_comments)
    matches = io.Store(src_matches)
    names = io.Store(src_names)
    proteins = io.Store(src_proteins)
    proteomes = io.Store(src_proteomes)

    logger.info("starting")
    n_proteins = 0
    chunk = []
    entries_with_matches = set()
    taxa_with_proteins = set()
    for protein_acc, protein_info in proteins:
        tax_id = protein_info["taxon"]

        name, other_names = names.get(protein_acc, (None, None))
        protein_matches = matches.get(protein_acc, {})

        # Enqueue protein
        chunk.append((
            protein_acc,
            protein_info["identifier"],
            name,
            "reviewed" if protein_info["is_reviewed"] else "unreviewed",
            protein_info["is_fragment"],
            protein_info["length"],
            comments.get(protein_acc, []),
            protein_matches,
            proteomes.get(protein_acc),
            taxa[tax_id]
        ))

        if len(chunk) == chunk_size:
            task_queue.put(("protein", chunk))
            chunk = []

        # Keep track of taxa having at least one protein
        taxa_with_proteins.add(tax_id)

        # Keep track of entries with at least one protein match
        entries_with_matches |= set(protein_matches.keys())

        n_proteins += 1
        if not n_proteins % 10000000:
            logger.info(f"{n_proteins:>12,}")

    if chunk:
        task_queue.put(("protein", chunk))

    logger.info(f"{n_proteins:>12,}")

    # Add entries without matches
    entries = set(mysql.entries.get_entries(url).keys())
    chunk = []
    for entry_acc in entries - entries_with_matches:
        chunk.append((entry_acc,))

    for i in range(0, len(chunk), chunk_size):
        task_queue.put(("entry", chunk[i:i+chunk_size]))

    # Add taxa without proteins
    chunk = []
    for taxon in taxa.values():
        if taxon["id"] not in taxa_with_proteins:
            chunk.append((taxon,))

    for i in range(0, len(chunk), chunk_size):
        task_queue.put(("taxon", chunk[i:i+chunk_size]))

    # Poison pill
    for _ in workers:
        task_queue.put(None)

    # Closing stores
    proteins.close()
    names.close()
    comments.close()
    proteomes.close()
    matches.close()

    n_docs = sum([done_queue.get() for _ in workers])

    # Wait for workers to finish
    for p in workers:
        p.join()

    # Delete loading file so Loaders know that all files are generated
    utils.set_ready(outdir)

    logger.info(f"complete: {n_docs:,} documents")


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
