import shelve
from interpro7dw.utils.store import BasicStore, Directory


def dump_to_tmp(xrefs: dict, stores: dict, outdir: Directory,
                compresslevel: int = 6):
    while xrefs:
        item_acc, item_xrefs = xrefs.popitem()

        try:
            store = stores[item_acc]
        except KeyError:
            file = outdir.mktemp()
            store = stores[item_acc] = BasicStore(file,
                                                  mode="a",
                                                  compresslevel=compresslevel)

        store.append(item_xrefs)


def unpack_pdb_matches(file: str) -> dict[str, set[str]]:
    entry2structures = {}
    with shelve.open(file) as matches:
        for pdb_chain, pdb_protein in matches.items():
            pdb_id, chain_id = pdb_chain.split("_")

            for entry_acc in pdb_protein["matches"]:
                try:
                    entry2structures[entry_acc].add(pdb_id)
                except KeyError:
                    entry2structures[entry_acc] = {pdb_id}

    return entry2structures
