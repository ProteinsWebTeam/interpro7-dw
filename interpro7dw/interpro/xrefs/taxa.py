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
    "proteomes": set(),
    "structures": {
        "all": set(),
        "databases": {}
    }
}


def _process(proteins_file: str, matches_file: str, proteomes_file: str,
             structures_file: str, uniprot2pdb_file: str, pdbmatches_file: str,
             start: str, stop: str | None, workdir: Directory,
             queue: mp.Queue):
    with open(uniprot2pdb_file, "rb") as fh:
        uniprot2pdb = pickle.load(fh)

    proteins_store = KVStore(proteins_file)
    matches_store = KVStore(matches_file)
    proteomes_store = KVStore(proteomes_file)

    i = 0
    all_taxa = set()
    tmp_stores = {}
    xrefs = {}
    for protein_acc, protein in proteins_store.range(start, stop):
        taxon_id = protein["taxid"]
        proteome_id = proteomes_store.get(protein_acc)

        if taxon_id in xrefs:
            taxon_xrefs = xrefs[taxon_id]
        else:
            taxon_xrefs = xrefs[taxon_id] = {}
            copy_dict(_BASE_XREFS, taxon_xrefs)
            all_taxa.add(taxon_id)

        taxon_xrefs["proteins"]["all"] += 1
        databases = set()
        for match in matches_store.get(protein_acc, []):
            match_acc = match["accession"]
            match_db = match["database"]

            try:
                db = taxon_xrefs["proteins"]["databases"][match_db]
            except KeyError:
                db = taxon_xrefs["proteins"]["databases"][match_db] = {
                    "count": 0,
                    "entries": {}
                }

            if match_db not in databases:
                # Counts the protein once per database
                databases.add(match_db)
                db["count"] += 1

            try:
                db["entries"][match_acc] += 1
            except KeyError:
                db["entries"][match_acc] = 1

        if proteome_id:
            taxon_xrefs["proteomes"].add(proteome_id)

        # Add structures, regardless of entry matches
        for pdb_chain in uniprot2pdb.get(protein_acc, {}):
            pdb_id, chain = pdb_chain.split("_")
            taxon_xrefs["structures"]["all"].add(pdb_id)

        i += 1
        if i == 1e6:
            dump_to_tmp(xrefs, tmp_stores, workdir)
            queue.put((False, i))
            i = 0

    del uniprot2pdb
    proteins_store.close()
    matches_store.close()
    proteomes_store.close()
    dump_to_tmp(xrefs, tmp_stores, workdir)
    queue.put((False, i))

    # Add mapping between taxa and PDB structures matched by entries
    taxon2pdb = unpack_taxon2pdb(structures_file)
    with shelve.open(pdbmatches_file, writeback=False) as chains:
        for taxon_id in all_taxa:

            taxon_xrefs = xrefs[taxon_id] = {}
            copy_dict(_BASE_XREFS, taxon_xrefs)
            taxon_structures = taxon_xrefs["structures"]["databases"]

            for pdb_chain in taxon2pdb.get(taxon_id, []):
                pdb_entry = chains.get(pdb_chain)
                if pdb_entry:
                    pdb_id, chain_id = pdb_chain.split("_")

                    for entry_acc, entry in pdb_entry["matches"].items():
                        database = entry["database"]

                        try:
                            db = taxon_structures[database]
                        except KeyError:
                            db = taxon_structures[database] = {}

                        try:
                            db[entry_acc].add(pdb_id)
                        except KeyError:
                            db[entry_acc] = {pdb_id}

    dump_to_tmp(xrefs, tmp_stores, workdir)
    queue.put((True, tmp_stores))


def export_xrefs(proteins_file: str, matches_file: str, proteomes_file: str,
                 structures_file: str, uniprot2pdb_file: str,
                 pdbmatches_file: str, taxa_file: str, output: str,
                 processes: int = 8, tempdir: str | None = None):
    """Export taxonomic cross-references, that is:
        - proteins
        - proteomes
        - InterPro entries and member database signatures
        - PDBe structures

    :param proteins_file: KVStore file of protein info
    :param matches_file: KVStore file of protein matches
    :param proteomes_file: KVStore file of protein-proteome mapping
    :param structures_file: File of PDBe structures.
    :param uniprot2pdb_file: File of UniProt-PDB mappings
    :param pdbmatches_file: File of PDB matches
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
        p = mp.Process(target=_process,
                       args=(proteins_file, matches_file, proteomes_file,
                             structures_file, uniprot2pdb_file,
                             pdbmatches_file, start, stop, workdir, queue))
        p.start()
        workers.append((p, workdir))

    taxon2stores = {}
    progress = 0
    milestone = step = 1e7
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
            progress += obj
            if progress >= milestone:
                logger.info(f"{progress:>15,.0f}")
                milestone += step

    logger.info(f"{progress:>15,.0f}")
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
            final_stores[taxon_id] = BasicStore(tempdir.mktemp(), mode="a",
                                                compresslevel=9)

    logger.info("propagating to ancestors")
    progress = 0
    total = len(taxon2stores)
    milestone = step = math.ceil(0.1 * total)
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

        progress += 1
        if progress == milestone:
            logger.info(f"{progress:>15,.0f} / {total:,}")
            milestone += step

    logger.info(f"{progress:>15,.0f} / {total:,}")

    # Delete workers' temp directories
    size = 0
    for p, workdir in workers:
        size += workdir.get_size()
        workdir.remove()

    logger.info("writing final file")
    with BasicStore(output, mode="w") as store:
        progress = 0
        total = len(final_stores)
        milestone = step = math.ceil(0.1 * total)
        while final_stores:
            taxon_id, taxon_store = final_stores.popitem()

            # Default cross-references
            taxon_xrefs = {}
            copy_dict(_BASE_XREFS, taxon_xrefs)

            # Merge cross-references from descendants
            for chunk in taxon_store:
                copy_dict(chunk, taxon_xrefs, concat_or_incr=True)

            store.write((taxon_id, taxon_xrefs))

            progress += 1
            if progress == milestone:
                logger.info(f"{progress:>15,.0f} / {total:,}")
                milestone += step

    logger.info(f"{progress:>15,.0f} / {total:,}")

    size += tempdir.get_size()
    logger.info(f"temporary files: {size / 1024 ** 2:.0f} MB")

    tempdir.remove()
    logger.info("done")
