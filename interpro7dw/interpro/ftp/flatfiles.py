import gzip
import os
from datetime import datetime

from interpro7dw.utils import logger
from interpro7dw.utils.store import Store, loadobj


_LIST = "entry.list"
_NAMES = "names.dat"
_SHORT_NAMES = "short_names.dat"
_INTERPRO2GO = "interpro2go"
_HIERARCHY = "ParentChildTreeFile.txt"
_UNIPROT2INTERPRO = "protein2ipr.dat.gz"


def _write_node(node, fh, level):
    fh.write(f"{'-'*2*level}{node['accession']}::{node['name']}::\n")

    for child in node["children"]:
        _write_node(child, fh, level+1)


def export(entries_file: str, matches_file: str, outdir: str):
    os.makedirs(outdir, exist_ok=True)

    logger.info("loading entries")
    entries = []
    integrated = {}
    for e in loadobj(entries_file).values():
        if e.database == "interpro" and e.is_public:
            entries.append(e)

            for signatures in e.integrates.values():
                for signature_acc in signatures:
                    integrated[signature_acc] = (e.accession, e.name)

    logger.info(f"writing {_LIST}")
    with open(os.path.join(outdir, _LIST), "wt") as fh:
        fh.write("ENTRY_AC\tENTRY_TYPE\tENTRY_NAME\n")

        for e in sorted(entries, key=lambda x: (x.type, x.accession)):
            fh.write(f"{e.accession}\t{e.type}\t{e.name}\n")

    logger.info(f"writing {_NAMES}")
    with open(os.path.join(outdir, _NAMES), "wt") as fh:
        for e in sorted(entries, key=lambda x: x.accession):
            fh.write(f"{e.accession}\t{e.name}\n")

    logger.info(f"writing {_SHORT_NAMES}")
    with open(os.path.join(outdir, _SHORT_NAMES), "wt") as fh:
        for e in sorted(entries, key=lambda x: x.accession):
            fh.write(f"{e.accession}\t{e.short_name}\n")

    logger.info(f"writing {_INTERPRO2GO}")
    with open(os.path.join(outdir, _INTERPRO2GO), "wt") as fh:
        fh.write(f"!date: {datetime.now():%Y/%m/%d %H:%M:%S}\n")
        fh.write("!Mapping of InterPro entries to GO\n")
        fh.write("!\n")

        for e in sorted(entries, key=lambda x: x.accession):
            for term in e.go_terms:
                fh.write(f"InterPro:{e.accession} {e.name} > "
                         f"GO:{term['name']} ; {term['identifier']}\n")

    logger.info(f"writing {_HIERARCHY}")
    with open(os.path.join(outdir, _HIERARCHY), "wt") as fh:
        for e in sorted(entries, key=lambda x: x.accession):
            root = e.hierarchy["accession"]
            if root == e.accession and e.hierarchy["children"]:
                _write_node(e.hierarchy, fh, level=0)

    logger.info("writing protein2ipr.dat.gz")
    filepath = os.path.join(outdir, "protein2ipr.dat.gz")
    with gzip.open(filepath, "wt") as fh, Store(matches_file, "r") as store:
        i = 0
        for uniprot_acc, protein_entries in store.items():
            matches = []
            for signature_acc in sorted(protein_entries):
                try:
                    entry_acc, entry_name = integrated[signature_acc]
                except KeyError:
                    # Not integrated signature or InterPro entry
                    continue

                locations = protein_entries[signature_acc]

                for loc in locations:
                    matches.append((
                        uniprot_acc,
                        entry_acc,
                        entry_name,
                        signature_acc,
                        # We do not consider fragmented locations
                        loc["fragments"][0]["start"],
                        max(f["end"] for f in loc["fragments"])
                    ))

            for m in sorted(matches):
                fh.write('\t'.join(map(str, m)) + '\n')

            i += 1
            if not i % 1e7:
                logger.debug(f"{i:>12,}")

        logger.debug(f"{i:>12,}")

    logger.info("complete")
