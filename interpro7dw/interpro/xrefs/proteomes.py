import math
import multiprocessing as mp
import pickle
from typing import Optional

from interpro7dw.interpro.utils import copy_dict
from interpro7dw.utils import logger
from interpro7dw.utils.store import BasicStore, Directory, KVStore
from .utils import dump_to_tmp


def _process(member2clan: dict, proteins_file: str, matches_file: str,
             proteomes_file: str, domorgs_file: str, structures_file: str,
             start: str, stop: Optional[str], workdir: Directory,
             queue: mp.Queue):
    with open(structures_file, "rb") as fh:
        protein2structures = pickle.load(fh)

    proteins_store = KVStore(proteins_file)
    matches_store = KVStore(matches_file)
    proteomes_store = KVStore(proteomes_file)
    domorgs_store = KVStore(domorgs_file)

    i = 0
    it = enumerate(proteomes_store.range(start, stop))
    tmp_stores = {}
    xrefs = {}
    for i, (protein_acc, proteome_id) in it:
        protein = proteins_store[protein_acc]
        taxon_id = protein["taxid"]

        if proteome_id in xrefs:
            proteome_xrefs = xrefs[proteome_id]
        else:
            proteome_xrefs = xrefs[proteome_id] = {
                "dom_orgs": set(),
                "entries": {},
                "proteins": 0,
                "sets": set(),
                "structures": set(),
                "taxa": set()
            }

        proteome_xrefs["proteins"] += 1

        try:
            domain = domorgs_store[protein_acc]
        except KeyError:
            pass
        else:
            proteome_xrefs["dom_orgs"].add(domain["id"])

        signatures, entries = matches_store.get(protein_acc, ({}, {}))
        for signature_acc, signature in signatures.items():
            database = signature["database"]

            if database in proteome_xrefs["entries"]:
                proteome_xrefs["entries"][database].add(signature_acc)
            else:
                proteome_xrefs["entries"][database] = {signature_acc}

            if signature_acc in member2clan:
                proteome_xrefs["sets"].add(member2clan[signature_acc])

        structures = protein2structures.get(protein_acc, {})
        proteome_xrefs["structures"] |= set(structures.keys())

        proteome_xrefs["taxa"].add(taxon_id)

        if (i + 1) % 1e5 == 0:
            dump_to_tmp(xrefs, tmp_stores, workdir)
            queue.put((False, 1e5))

    dump_to_tmp(xrefs, tmp_stores, workdir)
    proteins_store.close()
    matches_store.close()
    proteomes_store.close()
    domorgs_store.close()

    queue.put((False, i % 1e5))
    queue.put((True, tmp_stores))


def export_xrefs(clans_file: str, proteins_file: str, matches_file: str,
                 proteomes_file: str, domorgs_file: str, structures_file: str,
                 proteomeinfo_file: str, output: str, processes: int = 8,
                 tempdir: Optional[str] = None):
    """Export proteome cross-references, that is:
        - proteins
        - taxa
        - PDBe structures
        - InterPro entries and member database signatures
        - domain organisations
        - clans

    :param clans_file: File of clan information.
    :param proteins_file: KVStore file of protein info.
    :param matches_file: KVStore file of protein matches.
    :param proteomes_file: KVStore file of protein-proteome mapping.
    :param domorgs_file: KVStore file of domain organisations.
    :param structures_file: File of protein-structures mappings.
    :param proteomeinfo_file: File of reference proteomes.
    :param output: Output BasicStore file
    :param processes: Number of workers
    :param tempdir: Temporary directory
    """
    logger.info("loading clan members")
    member2clan = {}
    with open(clans_file, "rb") as fh:
        for clan_acc, clan in pickle.load(fh).items():
            for entry_acc, _, _ in clan["members"]:
                member2clan[entry_acc] = clan_acc

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
                       args=(member2clan, proteins_file, matches_file,
                             proteomes_file, domorgs_file, structures_file,
                             start, stop, workdir, queue))
        p.start()
        workers.append((p, workdir))

    proteome2stores = {}
    num_proteins = 0
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
            num_proteins += obj
            if num_proteins % 1e7 == 0:
                logger.info(f"{num_proteins:>15,.0f}")

    logger.info(f"{num_proteins:>15,.0f}")
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
            store.write((proteome_id, {
                "dom_orgs": set(),
                "entries": {},
                "proteins": 0,
                "sets": set(),
                "structures": set(),
                "taxa": set()
            }))

    size = 0
    for p, workdir, in workers:
        size += workdir.get_size()
        workdir.remove()

    logger.info(f"temporary files: {size / 1024 ** 2:.0f} MB")
    logger.info("done")
