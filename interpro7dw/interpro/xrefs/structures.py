import pickle

from interpro7dw.interpro.utils import overlaps_pdb_chain
from interpro7dw.utils import logger
from interpro7dw.utils.store import BasicStore, KVStore


def export_xrefs(clans_file: str, proteins_file: str, matches_file: str,
                 proteomes_file: str, domorgs_file: str,
                 structureinfo_file: str, structures_file: str, output: str):
    """Export PDBe structures cross-references, that is:
        - proteins
        - taxa
        - proteomes
        - InterPro entries and member database signatures
        - domain organisations
        - clans

    :param clans_file: File of clan information.
    :param proteins_file: KVStore file of protein info.
    :param matches_file: KVStore file of protein matches.
    :param proteomes_file: KVStore file of protein-proteome mapping.
    :param domorgs_file: KVStore file of domain organisations.
    :param structureinfo_file: File of PDBe structures.
    :param structures_file: File of protein-structures mapping.
    :param output: Output BasicStore file
    """
    logger.info("loading clan members")
    entry2clan = {}
    with open(clans_file, "rb") as fh:
        for clan_acc, clan in pickle.load(fh).items():
            for entry_acc, _, _ in clan["members"]:
                entry2clan[entry_acc] = clan_acc

    logger.info("loading PDBe structures")
    with open(structures_file, "rb") as fh:
        protein2structures = pickle.load(fh)

    logger.info("initializing cross-references")
    xrefs = {}
    with open(structureinfo_file, "rb") as fh:
        for pdbe_id in pickle.load(fh)["entries"]:
            xrefs[pdbe_id] = {
                "dom_orgs": set(),
                "entries": {},
                "proteomes": set(),
                "proteins": 0,
                "sets": set(),
                "taxa": set()
            }

    logger.info("iterating proteins")
    proteins_store = KVStore(proteins_file)
    matches_store = KVStore(matches_file)
    proteomes_store = KVStore(proteomes_file)
    domorgs_store = KVStore(domorgs_file)

    i = 0
    for i, protein_acc in enumerate(sorted(protein2structures)):
        structures = protein2structures[protein_acc]

        protein = proteins_store[protein_acc]
        taxon_id = protein["taxid"]
        proteome_id = proteomes_store.get(protein_acc)
        signatures, entries = matches_store.get(protein_acc, ({}, {}))

        try:
            domain = domorgs_store[protein_acc]
        except KeyError:
            domain_id = None
            domain_members = []
        else:
            domain_id = domain["id"]
            domain_members = domain["members"]

        for pdbe_id, chains in structures.items():
            struct_xrefs = xrefs[pdbe_id]

            if domain_id:
                struct_xrefs["dom_orgs"].add(domain_id)

            if proteome_id:
                struct_xrefs["proteomes"].add(proteome_id)

            struct_xrefs["proteins"] += 1
            struct_xrefs["taxa"].add(taxon_id)

            for obj in [signatures, entries]:
                for entry_acc, entry in obj.items():
                    database = entry["database"]
                    clan_acc = entry2clan.get(entry_acc)

                    for chain_id, segments in chains.items():
                        if overlaps_pdb_chain(entry["locations"], segments):
                            if database in struct_xrefs["entries"]:
                                struct_xrefs["entries"][database].add(entry_acc)
                            else:
                                struct_xrefs["entries"][database] = {entry_acc}

                            if clan_acc:
                                struct_xrefs["sets"].add(clan_acc)

                            break  # ignore other chains

        if (i + 1) % 1e4 == 0:
            logger.info(f"{i + 1:>15,}")

    logger.info(f"{i + 1:>15,}")

    proteins_store.close()
    matches_store.close()
    proteomes_store.close()
    domorgs_store.close()

    # Free memory
    entry2clan.clear()
    protein2structures.clear()

    with BasicStore(output, mode="w") as store:
        for obj in xrefs.items():
            store.write(obj)

    logger.info("done")
