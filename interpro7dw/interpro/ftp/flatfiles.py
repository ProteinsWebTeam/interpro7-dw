import gzip
import os
import pickle
from datetime import datetime

from interpro7dw.interpro.oracle.entries import Entry
from interpro7dw.utils import logger
from interpro7dw.utils.store import KVStore


_CITATION = "Paysan-Lafosse et al. (2023) Nucl. Acids Res. 51:D418–D427"
_LIST = "entry.list"
_NAMES = "names.dat"
_SHORT_NAMES = "short_names.dat"
_INTERPRO2GO = "interpro2go"
_HIERARCHY = "ParentChildTreeFile.txt"
_UNIPROT2INTERPRO = "protein2ipr.dat.gz"


def gen_hierarchy_tree(accession: str, entries: dict[str, Entry],
                       hierarchy: dict[str, list[str]], level: int = 0):
    entry = entries[accession]
    yield f"{'-'*2*level}{entry.accession}::{entry.name}::"

    for child_acc in sorted(hierarchy.get(entry.accession, [])):
        yield from gen_hierarchy_tree(child_acc, entries, hierarchy, level+1)


def export(entries_file: str, matches_file: str, outdir: str):
    os.makedirs(outdir, exist_ok=True)

    logger.info("loading entries")
    entries = {}
    with open(entries_file, "rb") as fh:
        for e in pickle.load(fh).values():
            if e.database.lower() == "interpro" and e.deletion_date is None:
                entries[e.accession] = e

    logger.info(f"writing {_LIST}")
    with open(os.path.join(outdir, _LIST), "wt") as fh:
        fh.write("ENTRY_AC\tENTRY_TYPE\tENTRY_NAME\n")

        for e in sorted(entries.values(), key=lambda x: (x.type, x.accession)):
            fh.write(f"{e.accession}\t{e.type}\t{e.name}\n")

    logger.info(f"writing {_NAMES}")
    with open(os.path.join(outdir, _NAMES), "wt") as fh:
        for e in sorted(entries.values(), key=lambda x: x.accession):
            fh.write(f"{e.accession}\t{e.name}\n")

    logger.info(f"writing {_SHORT_NAMES}")
    with open(os.path.join(outdir, _SHORT_NAMES), "wt") as fh:
        for e in sorted(entries.values(), key=lambda x: x.accession):
            fh.write(f"{e.accession}\t{e.short_name}\n")

    logger.info(f"writing {_INTERPRO2GO}")
    with open(os.path.join(outdir, _INTERPRO2GO), "wt") as fh:
        fh.write(f"!date: {datetime.now():%Y/%m/%d %H:%M:%S}\n")
        fh.write("!Mapping of InterPro entries to GO\n")
        fh.write("!external resource: http://www.ebi.ac.uk/interpro\n")
        fh.write(f"!citation: {_CITATION}\n")
        fh.write("!contact:interhelp@ebi.ac.uk")
        fh.write("!\n")

        for e in sorted(entries.values(), key=lambda x: x.accession):
            for term in e.go_terms:
                fh.write(f"InterPro:{e.accession} {e.name} > "
                         f"GO:{term['name']} ; {term['identifier']}\n")

    logger.info(f"writing {_HIERARCHY}")
    hierarchy = {}
    for e in entries.values():
        if e.parent:
            if e.parent in hierarchy:
                hierarchy[e.parent].append(e.accession)
            else:
                hierarchy[e.parent] = [e.accession]

    with open(os.path.join(outdir, _HIERARCHY), "wt") as fh:
        for e in sorted(entries.values(), key=lambda x: x.accession):
            e: Entry
            if e.accession in hierarchy:
                for line in gen_hierarchy_tree(accession=e.accession,
                                               entries=entries,
                                               hierarchy=hierarchy):
                    fh.write(line + '\n')

    logger.info("writing protein2ipr.dat.gz")
    filepath = os.path.join(outdir, "protein2ipr.dat.gz")
    with gzip.open(filepath, "wt") as fh, KVStore(matches_file) as store:
        i = 0
        for protein_acc, matches in store.items():
            ipr_matches = []
            for match in matches:
                if match["database"].lower() != "interpro":
                    # We only want member database signature
                    interpro_acc = match["entry"]

                    if interpro_acc is None:
                        # Not integrated
                        continue

                    for loc in match["locations"]:
                        ipr_matches.append((
                            protein_acc,
                            interpro_acc,
                            entries[interpro_acc]["name"],
                            match["accession"],
                            # We do not consider fragmented locations
                            loc["fragments"][0]["start"],
                            max(f["end"] for f in loc["fragments"])
                        ))

            for m in sorted(ipr_matches):
                fh.write('\t'.join(map(str, m)) + '\n')

            i += 1
            if i % 1e7 == 0:
                logger.debug(f"{i:>15,}")

        logger.debug(f"{i:>15,}")

    logger.info("done")
