import pickle
from typing import Optional

from interpro7dw.interpro.utils import copy_dict
from interpro7dw.utils import logger
from interpro7dw.utils.store import BasicStore, Directory, KVStore
from .utils import dump_to_tmp, load_protein2structures


def export_xrefs(clans_file: str, proteins_file: str, matches_file: str,
                 proteomes_file: str, domorgs_file: str, structures_file: str,
                 proteomeinfo_file: str, output: str,
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
    :param structures_file: File of PDBe structures.
    :param proteomeinfo_file: File of reference proteomes.
    :param output: Output BasicStore file
    :param tempdir: Temporary directory
    """
    logger.info("loading clan members")
    entry2clan = {}
    with open(clans_file, "rb") as fh:
        for clan_acc, clan in pickle.load(fh).items():
            for entry_acc, _, _ in clan["members"]:
                entry2clan[entry_acc] = clan_acc

    logger.info("loading PDBe structures")
    protein2structures = load_protein2structures(structures_file)

    logger.info("iterating proteins")
    proteins_store = KVStore(proteins_file)
    matches_store = KVStore(matches_file)
    proteomes_store = KVStore(proteomes_file)
    domorgs_store = KVStore(domorgs_file)

    i = 0
    tempdir = Directory(tempdir=tempdir)
    tmp_stores = {}
    xrefs = {}
    for i, (protein_acc, proteome_id) in enumerate(proteomes_store.items()):
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

        for obj in [signatures, entries]:
            for entry_acc, entry in obj.items():
                database = entry["database"]

                if database in proteome_xrefs["entries"]:
                    proteome_xrefs["entries"][database].add(entry_acc)
                else:
                    proteome_xrefs["entries"][database] = {entry_acc}

                if entry_acc in entry2clan:
                    proteome_xrefs["sets"].add(entry2clan[entry_acc])

        structures = protein2structures.get(protein_acc, {})
        proteome_xrefs["structures"] |= set(structures.keys())

        proteome_xrefs["taxa"].add(taxon_id)

        if (i + 1) % 1e4 == 0:
            dump_to_tmp(xrefs, tmp_stores, tempdir)

            if (i + 1) % 1e7 == 0:
                logger.info(f"{i + 1:>15,}")

    dump_to_tmp(xrefs, tmp_stores, tempdir)
    logger.info(f"{i + 1:>15,}")

    proteins_store.close()
    matches_store.close()
    proteomes_store.close()
    domorgs_store.close()

    # Free memory
    entry2clan.clear()
    protein2structures.clear()

    logger.info("writing final file")
    with open(proteomeinfo_file, "rb") as fh:
        proteomes = set(pickle.load(fh).keys())

    with BasicStore(output, mode="w") as store:
        for proteome_id, proteome_store in tmp_stores.items():
            proteomes.remove(proteome_id)

            # Merge cross-references
            proteome_xrefs = {}
            for xrefs in proteome_store:
                copy_dict(xrefs, proteome_xrefs, concat_or_incr=True)

            store.write((proteome_id, proteome_xrefs))

        logger.info(f"{len(proteomes)} proteomes without cross-references")
        for proteome_id in proteomes:
            store.add((proteome_id, {
                "dom_orgs": set(),
                "entries": {},
                "proteins": 0,
                "sets": set(),
                "structures": set(),
                "taxa": set()
            }))

    logger.info(f"temporary files: {tempdir.get_size() / 1024 ** 2:.0f} MB")
    tempdir.remove()

    logger.info("done")
