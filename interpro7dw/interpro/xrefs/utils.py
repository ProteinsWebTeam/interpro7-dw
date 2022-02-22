import pickle

from interpro7dw.utils.store import BasicStore, Directory


def dump_to_tmp(xrefs: dict, stores: dict, outdir: Directory):
    while xrefs:
        item_acc, item_xrefs = xrefs.popitem()

        try:
            store = stores[item_acc]
        except KeyError:
            file = outdir.mktemp()
            store = stores[item_acc] = BasicStore(file, mode="a")

        store.append(item_xrefs)


def load_protein2structures(file: str) -> dict:
    data = {}
    with open(file, "rb") as fh:
        # See pdbe.export_structures for structure of pickled object
        for e in pickle.load(fh)["entries"].values():
            pdbe_id = e["id"]
            for protein_acc, chains in e["proteins"].items():
                try:
                    data[protein_acc][pdbe_id] = chains
                except KeyError:
                    data[protein_acc] = {pdbe_id: chains}

    return data
