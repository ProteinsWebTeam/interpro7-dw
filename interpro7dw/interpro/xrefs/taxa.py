import math
import multiprocessing as mp
import pickle
from typing import Optional

from interpro7dw.interpro.utils import copy_dict, overlaps_pdb_chain
from interpro7dw.utils import logger
from interpro7dw.utils.store import BasicStore, Directory, KVStore
from .utils import dump_to_tmp, load_protein2structures


_BASE_XREFS = {
    "entries": {},
    "proteins": {
        "all": 0,
        "by_database": {},  # grouped by entry database
        "by_entries": {}  # grouped by entry
    },
    "proteomes": set(),
    "structures": {
        "all": set(),
        "by_entries": {}  # overlapping with these entries
    }
}


def _process_proteins(proteins_file: str, matches_file: str,
                      proteomes_file: str, structures_file: str,
                      start: str, stop: Optional[str], workdir: Directory,
                      queue: mp.Queue):
    with open(structures_file, "rb") as fh:
        protein2structures = pickle.load(fh)

    proteins_store = KVStore(proteins_file)
    matches_store = KVStore(matches_file)
    proteomes_store = KVStore(proteomes_file)

    i = 0
    it = enumerate(proteins_store.range(start, stop))
    tmp_stores = {}
    xrefs = {}
    for i, (protein_acc, protein) in it:
        taxon_id = protein["taxid"]
        proteome_id = proteomes_store.get(protein_acc)
        structures = protein2structures.get(protein_acc, {})
        signatures, entries = matches_store.get(protein_acc, ({}, {}))

        if taxon_id in xrefs:
            taxon_xrefs = xrefs[taxon_id]
        else:
            taxon_xrefs = xrefs[taxon_id] = {}
            copy_dict(_BASE_XREFS, taxon_xrefs)

        taxon_xrefs["proteins"]["all"] += 1
        if proteome_id:
            taxon_xrefs["proteomes"].add(proteome_id)

        # Add structures, regardless of entry matches
        taxon_xrefs["structures"]["all"] |= set(structures.keys())

        databases = set()
        for obj in [signatures, entries]:
            for entry_acc, entry in obj.items():
                database = entry["database"]

                if database in taxon_xrefs["entries"]:
                    taxon_xrefs["entries"][database].add(entry_acc)
                else:
                    taxon_xrefs["entries"][database] = {entry_acc}

                if database not in databases:
                    # Counts the protein once per database
                    databases.add(database)
                    if database in taxon_xrefs["proteins"]["by_database"]:
                        taxon_xrefs["proteins"]["by_database"][database] += 1
                    else:
                        taxon_xrefs["proteins"]["by_database"][database] = 1

                if entry_acc in taxon_xrefs["proteins"]["by_entries"]:
                    taxon_xrefs["proteins"]["by_entries"][entry_acc] += 1
                else:
                    taxon_xrefs["proteins"]["by_entries"][entry_acc] = 1

                by_entries = taxon_xrefs["structures"]["by_entries"]
                for pdbe_id, chains in structures.items():
                    for chain_id, segments in chains.items():
                        if overlaps_pdb_chain(entry["locations"], segments):
                            if entry_acc in by_entries:
                                by_entries[entry_acc].add(pdbe_id)
                            else:
                                by_entries[entry_acc] = {pdbe_id}

                            break  # Ignore other chains

        if (i + 1) % 1e6 == 0:
            dump_to_tmp(xrefs, tmp_stores, workdir)
            queue.put((False, 1e6))

    dump_to_tmp(xrefs, tmp_stores, workdir)
    proteins_store.close()
    matches_store.close()
    proteomes_store.close()

    queue.put((False, i % 1e6))
    queue.put((True, tmp_stores))


def export_xrefs2(proteins_file: str, matches_file: str, proteomes_file: str,
                  structures_file: str, taxa_file: str, output: str,
                  processes: int = 8, tempdir: Optional[str] = None):
    """Export taxonomic cross-references, that is:
        - proteins
        - proteomes
        - InterPro entries and member database signatures
        - PDBe structures

    :param proteins_file: KVStore file of protein info
    :param matches_file: KVStore file of protein matches
    :param proteomes_file: KVStore file of protein-proteome mapping
    :param structures_file: File of protein-structures mapping
    :param taxa_file: File of taxonomic information
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
        p = mp.Process(target=_process_proteins,
                       args=(proteins_file, matches_file, proteomes_file,
                             structures_file, start, stop, workdir, queue))
        p.start()
        workers.append((p, workdir))

    taxon2stores = {}
    num_proteins = 0
    work_done = 0
    while work_done < len(workers):
        is_done, obj = queue.get()
        if is_done:
            work_done += 1
            for taxon_id, taxon_store in obj.items():
                if taxon_id in taxon2stores:
                    taxon2stores[taxon_id].append(taxon_store)
                else:
                    taxon2stores[taxon_id] = [taxon_store]
        else:
            num_proteins += obj
            if num_proteins % 1e7 == 0:
                logger.info(f"{num_proteins:>15,.0f}")

    logger.info(f"{num_proteins:>15,.0f}")
    for p, workdir in workers:
        p.join()

    logger.info("loading taxonomy")
    lineages = {}
    final_stores = {}
    tempdir = Directory(tempdir=tempdir)
    with open(taxa_file, "rb") as fh:
        for taxon_id, taxon in pickle.load(fh).items():
            lineages[taxon_id] = taxon["lineage"]

            # Init (empty) final store
            final_stores[taxon_id] = BasicStore(tempdir.mktemp(), mode="a")

    logger.info("propagating to ancestors")
    while taxon2stores:
        taxon_id, taxon_stores = taxon2stores.popitem()

        # Merge cross-references from temporary store
        taxon_xrefs = {}
        for taxon_store in taxon_stores:
            for chunk in taxon_store:
                copy_dict(chunk, taxon_xrefs, concat_or_incr=True)

        # Propagate to lineage (including current taxon!) in final stores
        for node_id in lineages[taxon_id]:
            final_stores[node_id].append(taxon_xrefs)

    # Delete workers' temp directories
    size = 0
    for p, workdir in workers:
        size += workdir.get_size()
        workdir.remove()

    logger.info("writing final file")
    with BasicStore(output, mode="w") as store:
        while final_stores:
            taxon_id, taxon_store = final_stores.popitem()

            # Default cross-references
            taxon_xrefs = {}
            copy_dict(_BASE_XREFS, taxon_xrefs)

            # Merge cross-references from descendants
            for chunk in taxon_store:
                copy_dict(chunk, taxon_xrefs, concat_or_incr=True)

            store.write((taxon_id, taxon_xrefs))

    size += tempdir.get_size()
    logger.info(f"temporary files: {size / 1024 ** 2:.0f} MB")

    tempdir.remove()
    logger.info("done")


def export_xrefs(proteins_file: str, matches_file: str, proteomes_file: str,
                 structures_file: str, taxa_file: str, output: str,
                 tempdir: Optional[str] = None):
    """Export taxonomic cross-references, that is:
        - proteins
        - proteomes
        - InterPro entries and member database signatures
        - PDBe structures

    :param proteins_file: KVStore file of protein info.
    :param matches_file: KVStore file of protein matches.
    :param proteomes_file: KVStore file of protein-proteome mapping.
    :param structures_file: File of PDBe structures.
    :param taxa_file: File of taxonomic information.
    :param output: Output BasicStore file
    :param tempdir: Temporary directory
    """
    logger.info("loading PDBe structures")
    protein2structures = load_protein2structures(structures_file)

    logger.info("iterating proteins")
    proteins_store = KVStore(proteins_file)
    matches_store = KVStore(matches_file)
    proteomes_store = KVStore(proteomes_file)

    i = 0
    tempdir = Directory(tempdir=tempdir)
    tmp_stores = {}
    xrefs = {}
    base_xrefs = {
        "entries": {},
        "proteins": {
            "all": 0,
            "by_database": {},   # grouped by entry database
            "by_entries": {}     # grouped by entry
        },
        "proteomes": set(),
        "structures": {
            "all": set(),
            "by_entries": {}     # overlapping with these entries
        }
    }
    for i, (protein_acc, protein) in enumerate(proteins_store.items()):
        taxon_id = protein["taxid"]
        proteome_id = proteomes_store.get(protein_acc)
        structures = protein2structures.get(protein_acc, {})
        signatures, entries = matches_store.get(protein_acc, ({}, {}))

        if taxon_id in xrefs:
            taxon_xrefs = xrefs[taxon_id]
        else:
            taxon_xrefs = xrefs[taxon_id] = {}
            copy_dict(base_xrefs, taxon_xrefs)

        taxon_xrefs["proteins"]["all"] += 1
        if proteome_id:
            taxon_xrefs["proteomes"].add(proteome_id)

        # Add structures, regardless of entry matches
        taxon_xrefs["structures"]["all"] |= set(structures.keys())

        databases = set()
        for obj in [signatures, entries]:
            for entry_acc, entry in obj.items():
                database = entry["database"]

                if database in taxon_xrefs["entries"]:
                    taxon_xrefs["entries"][database].add(entry_acc)
                else:
                    taxon_xrefs["entries"][database] = {entry_acc}

                if database not in databases:
                    # Counts the protein once per database
                    databases.add(database)
                    if database in taxon_xrefs["proteins"]["by_database"]:
                        taxon_xrefs["proteins"]["by_database"][database] += 1
                    else:
                        taxon_xrefs["proteins"]["by_database"][database] = 1

                if entry_acc in taxon_xrefs["proteins"]["by_entries"]:
                    taxon_xrefs["proteins"]["by_entries"][entry_acc] += 1
                else:
                    taxon_xrefs["proteins"]["by_entries"][entry_acc] = 1

                by_entries = taxon_xrefs["structures"]["by_entries"]
                for pdbe_id, chains in structures.items():
                    for chain_id, segments in chains.items():
                        if overlaps_pdb_chain(entry["locations"], segments):
                            if entry_acc in by_entries:
                                by_entries[entry_acc].add(pdbe_id)
                            else:
                                by_entries[entry_acc] = {pdbe_id}

                            break  # Ignore other chains

        if (i + 1) % 1e6 == 0:
            dump_to_tmp(xrefs, tmp_stores, tempdir)

            if (i + 1) % 1e7 == 0:
                logger.info(f"{i + 1:>15,}")

    dump_to_tmp(xrefs, tmp_stores, tempdir)
    logger.info(f"{i + 1:>15,}")

    proteins_store.close()
    matches_store.close()
    proteomes_store.close()

    # Free memory
    protein2structures.clear()

    logger.info("loading taxonomy")
    lineages = {}
    final_stores = {}
    with open(taxa_file, "rb") as fh:
        for taxon_id, taxon in pickle.load(fh).items():
            lineages[taxon_id] = taxon["lineage"]

            # Init (empty) final store
            final_stores[taxon_id] = BasicStore(tempdir.mktemp(), mode="a")

    logger.info("propagating to ancestors")
    while tmp_stores:
        taxon_id, taxon_store = tmp_stores.popitem()

        # Merge cross-references from temporary store
        taxon_xrefs = {}
        for chunk in taxon_store:
            copy_dict(chunk, taxon_xrefs, concat_or_incr=True)

        # Propagate to lineage (including current taxon!) in final stores
        for node_id in lineages[taxon_id]:
            final_stores[node_id].append(taxon_xrefs)

    logger.info("writing final file")
    with BasicStore(output, mode="w") as store:
        while final_stores:
            taxon_id, taxon_store = final_stores.popitem()

            # Default cross-references
            taxon_xrefs = {}
            copy_dict(base_xrefs, taxon_xrefs)

            # Merge cross-references from descendants
            for chunk in taxon_store:
                copy_dict(chunk, taxon_xrefs, concat_or_incr=True)

            store.write((taxon_id, taxon_xrefs))

    logger.info(f"temporary files: {tempdir.get_size() / 1024 ** 2:.0f} MB")
    tempdir.remove()

    logger.info("done")
