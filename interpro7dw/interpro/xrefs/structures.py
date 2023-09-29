import pickle
import shelve

from interpro7dw.utils import logger
from interpro7dw.utils.store import BasicStore, KVStore


def export_xrefs(clans_file: str, proteins_file: str, proteomes_file: str,
                 domorgs_file: str, structures_file: str, pdbmatches_file: str,
                 uniprot2pdb_file: str, output: str):
    """Export PDBe structures cross-references, that is:
        - proteins
        - taxa
        - proteomes
        - InterPro entries and member database signatures
        - domain organisations
        - clans

    :param clans_file: File of clan information.
    :param proteins_file: KVStore file of protein info.
    :param proteomes_file: KVStore file of protein-proteome mapping.
    :param domorgs_file: KVStore file of domain organisations.
    :param structures_file: File of PDBe structures.
    :param pdbmatches_file: File of PDB matches.
    :param uniprot2pdb_file: File of protein-structures mapping.
    :param output: Output BasicStore file
    """
    logger.info("loading clan members")
    entry2clan = {}
    with open(clans_file, "rb") as fh:
        for clan_acc, clan in pickle.load(fh).items():
            for entry_acc, _, _, _, _ in clan["members"]:
                entry2clan[entry_acc] = clan_acc

    logger.info("loading protein-structures mappings")
    with open(uniprot2pdb_file, "rb") as fh:
        uniprot2pdb = pickle.load(fh)

    logger.info("initializing cross-references")
    xrefs = {}
    with open(structures_file, "rb") as fh:
        for pdb_id in pickle.load(fh):
            xrefs[pdb_id] = {
                "dom_orgs": set(),
                "entries": {},
                "proteomes": set(),
                "proteins": 0,
                "sets": set(),
                "taxa": set()
            }

    with shelve.open(pdbmatches_file, writeback=False) as d:
        for pdb_chain, pdb_entry in d.items():
            pdb_id, chain_id = pdb_chain.split("_")

            try:
                obj = xrefs[pdb_id]
            except KeyError:
                """
                Match against a structure that is not NMR, X-ray, or EM.
                Only NMR, X-ray, and EM structures are in `structures_file`
                """
                continue

            databases = obj["entries"]

            for entry_acc, match in pdb_entry["matches"].items():
                database = match["database"]

                try:
                    databases[database].add(entry_acc)
                except KeyError:
                    databases[database] = {entry_acc}

                try:
                    clan_acc = entry2clan[entry_acc]
                except KeyError:
                    pass
                else:
                    obj["sets"].add(clan_acc)

    logger.info("iterating proteins")
    proteins_store = KVStore(proteins_file)
    proteomes_store = KVStore(proteomes_file)
    domorgs_store = KVStore(domorgs_file)

    i = 0
    for i, protein_acc in enumerate(sorted(uniprot2pdb)):
        structures = uniprot2pdb[protein_acc]

        try:
            protein = proteins_store[protein_acc]
        except KeyError:
            continue

        taxon_id = protein["taxid"]
        proteome_id = proteomes_store.get(protein_acc)

        try:
            domain = domorgs_store[protein_acc]
        except KeyError:
            domain_id = None
        else:
            domain_id = domain["id"]

        for pdb_chain in structures:
            pdb_id, chain = pdb_chain.split("_")
            struct_xrefs = xrefs[pdb_id]

            if domain_id:
                struct_xrefs["dom_orgs"].add(domain_id)

            if proteome_id:
                struct_xrefs["proteomes"].add(proteome_id)

            struct_xrefs["proteins"] += 1
            struct_xrefs["taxa"].add(taxon_id)

        if (i + 1) % 1e4 == 0:
            logger.info(f"{i + 1:>15,}")

    logger.info(f"{i + 1:>15,}")

    proteins_store.close()
    proteomes_store.close()
    domorgs_store.close()

    # Free memory
    entry2clan.clear()

    with BasicStore(output, mode="w") as store:
        for obj in xrefs.items():
            store.write(obj)

    logger.info("done")
