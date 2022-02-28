import hashlib
import math
import multiprocessing as mp
from typing import Optional

from interpro7dw.utils import logger
from interpro7dw.utils.store import BasicStore, Directory, KVStoreBuilder
from interpro7dw.utils.store import KVStore


def _process(proteins_file: str, matches_file: str, start: str,
             stop: Optional[str], output: str, queue: mp.Queue):
    proteins_store = KVStore(proteins_file)
    matches_store = KVStore(matches_file)

    all_domains = {}
    i = 0
    with BasicStore(output, "w") as store:
        it = matches_store.range(start, stop)
        for protein_acc, (signatures, entries) in it:
            i += 1
            if i == 1e6:
                queue.put(i)
                i = 0

            protein = proteins_store[protein_acc]

            locations = []
            for signature_acc, signature in signatures.items():
                if signature["database"].lower() != "pfam":
                    continue  # only consider Pfam matches

                for loc in signature["locations"]:
                    locations.append({
                        "pfam": signature_acc,
                        "interpro": signature["entry"],
                        # We do not consider fragmented locations
                        "start": loc["fragments"][0]["start"],
                        "end": max(f["end"] for f in loc["fragments"])
                    })

            if not locations:
                continue  # No Pfam matches: no domain organisation

            # Sort by position
            locations.sort(key=lambda l: (l["start"], l["end"]))

            domains = []
            members = set()
            for loc in locations:
                pfam_acc = loc["pfam"]
                interpro_acc = loc["interpro"]

                if interpro_acc:
                    domains.append(f"{pfam_acc}:{interpro_acc}")
                    members.add(interpro_acc)
                else:
                    domains.append(pfam_acc)

                members.add(pfam_acc)

            dom_str = '-'.join(domains)
            dom_id = hashlib.sha1(dom_str.encode("utf-8")).hexdigest()
            store.write((protein_acc, dom_id, members))

            # string (YYYY-MM-DD) which is enough to compare dates
            date = protein["date"]
            length = protein["length"]

            # Selects the oldest protein to represent
            # this domain organisation.
            try:
                domain = all_domains[dom_id]
            except KeyError:
                all_domains[dom_id] = {
                    "id": dom_id,
                    "key": dom_str,
                    "protein": protein_acc,
                    "date": date,
                    "length": length,
                    "locations": locations,
                    "count": 1
                }
            else:
                domain["count"] += 1
                if date < domain["date"]:
                    domain.update({
                        "protein": protein_acc,
                        "date": date,
                        "length": length,
                        "locations": locations,
                    })

    queue.put(i)
    queue.put(list(all_domains.values()))


def export(proteins_file: str, matches_file: str, output: str,
           processes: int = 8, tempdir: Optional[str] = None):
    """Calculates and exports the domain architectures/organisations of
    UniProt entries based on the Pfam matches.

    :param proteins_file: The file containing protein information.
    :param matches_file: The file containing protein matches.
    :param output: The output domain organisation file.
    :param processes: Number of workers.
    :param tempdir: Temporary directory .
    """
    logger.info("iterating proteins")
    tempdir = Directory(tempdir=tempdir)

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

        tempfile = tempdir.mktemp()
        p = mp.Process(target=_process,
                       args=(proteins_file, matches_file, start, stop,
                             tempfile, queue))
        p.start()
        workers.append((p, tempfile))

    progress = 0
    milestone = step = 1e7
    done = 0
    domains = {}
    while done < len(workers):
        obj = queue.get()
        if isinstance(obj, list):
            done += 1

            for domain in obj:
                domain_id = domain["id"]
                if domain_id not in domains:
                    domains[domain_id] = domain
                    continue

                # Add number of proteins
                domains[domain_id]["count"] += domain["count"]

                if domain["date"] < domains[domain_id]["date"]:
                    # Update representative protein
                    domains[domain_id].update({
                        "protein": domain["protein"],
                        "date": domain["date"],
                        "length": domain["length"],
                        "locations": domain["locations"],
                    })
        else:
            progress += obj
            if progress >= milestone:
                logger.info(f"{progress:>15,.0f}")
                milestone += step

    logger.info(f"{progress:>15,.0f}")

    logger.info("exporting domain organisations")
    with KVStoreBuilder(output, keys=[], cachesize=10000) as dst:
        i = 0
        for p, tempfile in workers:
            p.join()

            with BasicStore(tempfile, "r") as src:
                for protein_acc, domain_id, members in src:
                    domain = domains[domain_id].copy()
                    domain["members"] = members
                    dst.append(protein_acc, domain)

                    i += 1
                    if i % 1e7 == 0:
                        logger.info(f"{i:>15,}")

        dst.close()
        logger.info(f"{i:>15,}")

    logger.info(f"temporary files: {tempdir.get_size() / 1024 ** 2:.0f} MB")
    tempdir.remove()

    logger.info("done")
