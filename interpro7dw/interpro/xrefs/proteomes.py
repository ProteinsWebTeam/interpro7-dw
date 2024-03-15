import math
import multiprocessing as mp
import pickle
import shelve

from interpro7dw.interpro.utils import copy_dict
from interpro7dw.utils import logger
from interpro7dw.utils.store import BasicStore, Directory, KVStore
from .utils import dump_to_tmp, unpack_taxon2pdb


_BASE_XREFS = {
    "proteins": {
        "all": 0,
        "databases": {}
    },
    "structures": {
        "all": set(),
        "databases": {}
    }
}


def _process(proteins_file: str, matches_file: str,
             proteomes_file: str, structures_file: str, uniprot2pdb_file: str,
             pdbmatches_file: str, start: str, stop: str | None,
             workdir: Directory, queue: mp.Queue):
    with open(uniprot2pdb_file, "rb") as fh:
        uniprot2pdb = pickle.load(fh)

    proteins_store = KVStore(proteins_file)
    matches_store = KVStore(matches_file)
    proteomes_store = KVStore(proteomes_file)

    i = 0
    tmp_stores = {}
    xrefs = {}
    proteome2taxon = {}
    for protein_acc, proteome_id in proteomes_store.range(start, stop):
        if proteome_id not in proteome2taxon:
            protein = proteins_store[protein_acc]
            proteome2taxon[proteome_id] = protein["taxid"]

        if proteome_id in xrefs:
            proteome_xrefs = xrefs[proteome_id]
        else:
            proteome_xrefs = xrefs[proteome_id] = {}
            copy_dict(_BASE_XREFS, proteome_xrefs)

        proteome_xrefs["proteins"]["all"] += 1
        signatures, entries = matches_store.get(protein_acc, ({}, {}))
        databases = set()
        for obj in [signatures, entries]:
            for entry_acc, entry in obj.items():
                database = entry["database"]

                try:
                    db = proteome_xrefs["proteins"]["databases"][database]
                except KeyError:
                    db = proteome_xrefs["proteins"]["databases"][database] = {
                        "count": 0,
                        "entries": {}
                    }

                if database not in databases:
                    # Counts the protein once per database
                    databases.add(database)
                    db["count"] += 1

                try:
                    db["entries"][entry_acc] += 1
                except KeyError:
                    db["entries"][entry_acc] = 1

        # Add structures, regardless of entry matches
        for pdb_chain in uniprot2pdb.get(protein_acc, {}):
            pdb_id, chain = pdb_chain.split("_")
            proteome_xrefs["structures"]["all"].add(pdb_id)

        i += 1
        if i == 1e5:
            dump_to_tmp(xrefs, tmp_stores, workdir)
            queue.put((False, i))
            i = 0

    del uniprot2pdb
    proteins_store.close()
    matches_store.close()
    proteomes_store.close()
    dump_to_tmp(xrefs, tmp_stores, workdir)
    queue.put((False, i))

    # Add mapping between proteomes and PDB structures matched by entries
    taxon2pdb = unpack_taxon2pdb(structures_file)
    with shelve.open(pdbmatches_file, writeback=False) as chains:
        for proteome_id, taxon_id in proteome2taxon.items():

            proteome_xrefs = xrefs[proteome_id] = {}
            copy_dict(_BASE_XREFS, proteome_xrefs)
            proteome_structures = proteome_xrefs["structures"]["databases"]

            for pdb_chain in taxon2pdb.get(taxon_id, []):
                pdb_entry = chains.get(pdb_chain)
                if pdb_entry:
                    pdb_id, chain_id = pdb_chain.split("_")

                    for entry_acc, entry in pdb_entry["matches"].items():
                        database = entry["database"]

                        try:
                            db = proteome_structures[database]
                        except KeyError:
                            db = proteome_structures[database] = {}

                        try:
                            db[entry_acc].add(pdb_id)
                        except KeyError:
                            db[entry_acc] = {pdb_id}

    dump_to_tmp(xrefs, tmp_stores, workdir)
    queue.put((True, tmp_stores))


def export_xrefs(proteins_file: str, matches_file: str,
                 proteomes_file: str, structures_file: str, uniprot2pdb_file: str,
                 pdbmatches_file: str, proteomeinfo_file: str, output: str,
                 processes: int = 8, tempdir: str | None = None):
    """Export proteome cross-references, that is:
        - proteins
        - taxa
        - PDBe structures
        - InterPro entries and member database signatures
        - domain organisations
        - clans

    :param proteins_file: KVStore file of protein info.
    :param matches_file: KVStore file of protein matches.
    :param proteomes_file: KVStore file of protein-proteome mapping.
    :param structures_file: File of PDBe structures.
    :param uniprot2pdb_file: File of protein-structures mappings.
    :param pdbmatches_file: File of PDB matches.
    :param proteomeinfo_file: File of reference proteomes.
    :param output: Output BasicStore file
    :param processes: Number of workers
    :param tempdir: Temporary directory
    """
    logger.info("iterating proteins")
    processes = max(1, processes - 1)
    with KVStore(proteins_file) as store:
        keys = store.get_keys()

    chunksize = math.ceil(len(keys) / processes)
    queue = mp.Queue()
    workers = []
    for i in range(processes):
        start = keys[i * chunksize]
        try:
            stop = keys[(i + 1) * chunksize]
        except IndexError:
            stop = None

        workdir = Directory(tempdir=tempdir)
        p = mp.Process(target=_process,
                       args=(proteins_file, matches_file,
                             proteomes_file, structures_file, uniprot2pdb_file,
                             pdbmatches_file, start, stop, workdir, queue))
        p.start()
        workers.append((p, workdir))

    proteome2stores = {}
    progress = 0
    milestone = step = 1e7
    work_done = 0
    while work_done < len(workers):
        is_done, obj = queue.get()
        if is_done:
            work_done += 1
            for proteome_id, store in obj.items():
                if proteome_id in proteome2stores:
                    proteome2stores[proteome_id].append(store)
                else:
                    proteome2stores[proteome_id] = [store]
        else:
            progress += obj
            if progress >= milestone:
                logger.info(f"{progress:>15,.0f}")
                milestone += step

    logger.info(f"{progress:>15,.0f}")
    for p, workdir in workers:
        p.join()

    logger.info("writing final file")
    with open(proteomeinfo_file, "rb") as fh:
        proteomes = set(pickle.load(fh).keys())

    with BasicStore(output, mode="w") as store:
        while proteome2stores:
            proteome_id, proteome_stores = proteome2stores.popitem()
            proteomes.remove(proteome_id)

            # Merge cross-references
            proteome_xrefs = {}
            for proteome_store in proteome_stores:
                for xrefs in proteome_store:
                    copy_dict(xrefs, proteome_xrefs, concat_or_incr=True)

            store.write((proteome_id, proteome_xrefs))

        logger.info(f"{len(proteomes)} proteomes without cross-references")
        for proteome_id in proteomes:
            store.write((proteome_id, _BASE_XREFS))

    size = 0
    for p, workdir, in workers:
        size += workdir.get_size()
        workdir.remove()

    logger.info(f"temporary files: {size / 1024 ** 2:.0f} MB")
    logger.info("done")
