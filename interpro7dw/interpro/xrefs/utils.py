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


def unpack_pdb_matches(file: str) -> dict[str, dict]:
    entry2structures = {}
    with shelve.open(file) as matches:
        for pdb_chain, pdb_protein in matches.items():
            pdb_id, chain_id = pdb_chain.split("_")
            length = pdb_protein["length"]

            for entry_acc, match in pdb_protein["matches"].items():
                coverage = [0] * length

                for location in match["locations"]:
                    for fragment in location["fragments"]:
                        for i in range(fragment["start"] - 1, fragment["end"]):
                            coverage[i] = 1

                try:
                    entry_structures = entry2structures[entry_acc]
                except KeyError:
                    entry_structures = entry2structures[entry_acc] = {}

                try:
                    structure = entry_structures[pdb_id]
                except KeyError:
                    structure = entry_structures[pdb_id] = {
                        "length": 0,
                        "coverage": 0
                    }

                structure["length"] += length
                structure["coverage"] += sum(coverage)

    return entry2structures
