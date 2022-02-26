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
