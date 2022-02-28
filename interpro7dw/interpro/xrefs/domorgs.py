import hashlib
import os

from interpro7dw.utils import logger
from interpro7dw.utils.store import BasicStore, KVStoreBuilder, KVStore


def export(proteins_file: str, matches_file: str, output: str):
    """Calculates and exports the domain architectures/organisations of
    UniProt entries based on the Pfam matches.

    :param proteins_file: The file containing protein information.
    :param matches_file: The file containing protein matches.
    :param output: The output domain organisation file.
    """
    logger.info("starting")

    tmpfile = f"{output}.tmp"
    ps = KVStore(proteins_file)
    ms = KVStore(matches_file)
    ds = BasicStore(tmpfile, "w")
    all_domains = {}
    i = 0
    for i, (prot_acc, (signatures, entries)) in enumerate(ms.items()):
        protein = ps[prot_acc]

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

        if (i + 1) % 1e7 == 0:
            logger.info(f"{i + 1:>15,}")

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
        ds.write((prot_acc, dom_id, members))

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
                "protein": prot_acc,
                "date": date,
                "length": length,
                "locations": locations,
                "count": 1
            }
        else:
            domain["count"] += 1
            if date < domain["date"]:
                domain.update({
                    "protein": prot_acc,
                    "date": date,
                    "length": length,
                    "locations": locations,
                })

    logger.info(f"{i + 1:>15,}")

    for s in [ps, ms, ds]:
        s.close()

    logger.info("exporting domain organisations")
    src = BasicStore(tmpfile, "r")
    dst = KVStoreBuilder(output, keys=[], cachesize=10000)

    i = 0
    for i, (prot_acc, dom_id, members) in enumerate(src):
        domain = all_domains[dom_id].copy()
        domain["members"] = members
        dst.append(prot_acc, domain)

        if (i + 1) % 1e7 == 0:
            logger.info(f"{i + 1:>15,}")

    logger.info(f"{i + 1:>15,}")

    src.close()
    dst.close()
    os.unlink(tmpfile)

    logger.info("done")
