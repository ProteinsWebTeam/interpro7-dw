import math
import multiprocessing as mp
import pickle
from typing import Optional

from interpro7dw.interpro.utils import copy_dict
from interpro7dw.utils import logger
from interpro7dw.utils.store import BasicStore, Directory, KVStore
from .utils import dump_to_tmp, unpack_pdb_matches


def _process(member2clan: dict, proteins_file: str, matches_file: str,
             proteomes_file: str, domorgs_file: str,
             pdb2matches_file: str, start: str, stop: Optional[str],
             workdir: Directory, queue: mp.Queue):
    entry2structures = unpack_pdb_matches(pdb2matches_file)
    proteins_store = KVStore(proteins_file)
    matches_store = KVStore(matches_file)
    proteomes_store = KVStore(proteomes_file)
    domorgs_store = KVStore(domorgs_file)

    i = 0
    tmp_stores = {}
    xrefs = {}
    for protein_acc, (signatures, entries) in matches_store.range(start, stop):
        protein = proteins_store[protein_acc]
        taxon_id = protein["taxid"]
        proteome_id = proteomes_store.get(protein_acc)

        try:
            domain = domorgs_store[protein_acc]
        except KeyError:
            domain_id = None
            domain_members = []
        else:
            domain_id = domain["id"]
            domain_members = domain["members"]

        for signature_acc, signature in signatures.items():
            if signature_acc not in member2clan:
                continue

            clan_acc, database = member2clan[signature_acc]
            if clan_acc in xrefs:
                clan_xrefs = xrefs[clan_acc]
            else:
                clan_xrefs = xrefs[clan_acc] = {
                    "dom_orgs": set(),
                    "entries": {
                        "all": set()
                    },
                    "proteins": [],
                    "proteomes": set(),
                    "structures": set(),
                    "taxa": set()
                }

            if signature_acc in domain_members:
                clan_xrefs["dom_orgs"].add(domain_id)

            clan_xrefs["entries"]["all"].add(signature_acc)
            if database in clan_xrefs["entries"]:
                clan_xrefs["entries"][database].add(signature_acc)
            else:
                clan_xrefs["entries"][database] = {signature_acc}

            clan_xrefs["proteins"].append(protein_acc)

            if proteome_id:
                clan_xrefs["proteomes"].add(proteome_id)

            """
            Use `pop()` instead of `get()` so we add structures once
            per signature 
            """
            for pdb_id in entry2structures.pop(signature_acc, []):
                clan_xrefs["structures"].add(pdb_id)

            clan_xrefs["taxa"].add(taxon_id)

        i += 1
        if i == 1e5:
            dump_to_tmp(xrefs, tmp_stores, workdir)
            queue.put((False, i))
            i = 0

    dump_to_tmp(xrefs, tmp_stores, workdir)
    proteins_store.close()
    matches_store.close()
    proteomes_store.close()
    domorgs_store.close()

    queue.put((False, i))
    queue.put((True, tmp_stores))


def export_xrefs(clans_file: str, proteins_file: str, matches_file: str,
                 proteomes_file: str, domorgs_file: str,
                 pdb2matches_file: str, output: str,
                 processes: int = 8, tempdir: Optional[str] = None):
    logger.info("loading clan members")
    clans = {}
    member2clan = {}
    with open(clans_file, "rb") as fh:
        for clan_acc, clan in pickle.load(fh).items():
            database = clan["database"]
            members = set()

            for entry_acc, _, _, _, _ in clan["members"]:
                member2clan[entry_acc] = (clan_acc, database)
                members.add(entry_acc)

            clans[clan_acc] = (database, members)

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
                             proteomes_file, domorgs_file, pdb2matches_file,
                             start, stop, workdir, queue))
        p.start()
        workers.append((p, workdir))

    clan2stores = {}
    progress = 0
    milestone = step = 1e7
    work_done = 0
    while work_done < len(workers):
        is_done, obj = queue.get()
        if is_done:
            work_done += 1
            for clan_acc, clan_store in obj.items():
                if clan_acc in clan2stores:
                    clan2stores[clan_acc].append(clan_store)
                else:
                    clan2stores[clan_acc] = [clan_store]
        else:
            progress += obj
            if progress >= milestone:
                logger.info(f"{progress:>15,.0f}")
                milestone += step

    logger.info(f"{progress:>15,.0f}")
    for p, workdir in workers:
        p.join()

    logger.info("writing final file")
    with BasicStore(output, mode="w") as store:
        progress = 0
        total = len(clan2stores)
        milestone = step = math.ceil(0.1 * total)
        while clan2stores:
            clan_acc, clan_stores = clan2stores.popitem()
            clans.pop(clan_acc)

            # Merge cross-references
            clan_xrefs = {}
            for clan_store in clan_stores:
                for xrefs in clan_store:
                    copy_dict(xrefs, clan_xrefs, concat_or_incr=True)

            store.write((clan_acc, clan_xrefs))

            progress += 1
            if progress == milestone:
                logger.info(f"{progress:>15,.0f} / {total:,}")
                milestone += step

        logger.info(f"{progress:>15,.0f} / {total:,}")

        logger.info(f"{len(clans)} clans without cross-references")
        for clan_acc, (database, members) in clans.items():
            store.write((clan_acc, {
                "dom_orgs": set(),
                "entries": {
                    "all": members,
                    database: members
                },
                "proteins": [],
                "proteomes": set(),
                "structures": set(),
                "taxa": set()
            }))

    size = 0
    for p, workdir, in workers:
        size += workdir.get_size()
        workdir.remove()

    logger.info(f"temporary files: {size / 1024 ** 2:.0f} MB")
    logger.info("done")
