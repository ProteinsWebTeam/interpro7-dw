import pickle
from typing import Optional

from interpro7dw.interpro.utils import copy_dict, overlaps_pdb_chain
from interpro7dw.utils import logger
from interpro7dw.utils.store import BasicStore, Directory, KVStore
from .utils import dump_to_tmp, load_protein2structures


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
