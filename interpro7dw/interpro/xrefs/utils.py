import pickle
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


def unpack_entry2structures(file: str) -> dict[str, dict]:
    entry2structures = {}
    with shelve.open(file, writeback=False) as matches:
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


def unpack_struct2entries(file: str) -> dict[str, dict[str, set[str]]]:
    struct2entries = {}
    with shelve.open(file, writeback=False) as d:
        for pdb_chain, pdb_entry in d.items():
            pdb_id, chain_id = pdb_chain.split("_")

            for entry_acc, entry in pdb_entry["matches"].items():
                database = entry["database"]

                if pdb_id in struct2entries:
                    dbs = struct2entries[pdb_id]
                    try:
                        dbs[database].add(entry_acc)
                    except KeyError:
                        dbs[database] = {entry_acc}
                else:
                    struct2entries[pdb_id] = {database: {entry_acc}}

    return struct2entries


def unpack_taxon2pdb(file: str) -> dict[str, set[str]]:
    taxon2pdb = {}
    with open(file, "rb") as fh:
        for s in pickle.load(fh).values():
            pdb_id = s["id"]
            for chain_id, taxon_id in s["taxonomy"].items():
                pdb_chain = f"{pdb_id}_{chain_id}"
                try:
                    taxon2pdb[taxon_id].add(pdb_chain)
                except KeyError:
                    taxon2pdb[taxon_id] = {pdb_chain}

    return taxon2pdb
