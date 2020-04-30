# -*- coding: utf-8 -*-

import gzip
import os
from datetime import datetime

from interpro7dw import logger
from interpro7dw.utils import Store, loadobj


def _write_node(node, fh, level):
    fh.write(f"{'-'*2*level}{node['accession']}::{node['name']}\n")

    for child in node["children"]:
        _write_node(child, fh, level+1)


def export(p_entries: str, p_uniprot2matches: str, outdir: str):
    logger.info("loading entries")
    entries = []
    integrated = {}
    for e in loadobj(p_entries).values():
        if e.database == "interpro" and not e.is_deleted:
            entries.append(e)

            for signatures in e.integrates.values():
                for signature_acc in signatures:
                    integrated[signature_acc] = (e.accession, e.name)

    logger.info("writing entry.list")
    with open(os.path.join(outdir, "entry.list"), "wt") as fh:
        fh.write("ENTRY_AC\tENTRY_TYPE\tENTRY_NAME\n")

        for e in sorted(entries, key=lambda e: (e.type, e.accession)):
            fh.write(f"{e.accession}\t{e.type}\t{e.name}\n")

    logger.info("writing names.dat")
    with open(os.path.join(outdir, "names.dat"), "wt") as fh:
        for e in sorted(entries, key=lambda e: e.accession):
            fh.write(f"{e.accession}\t{e.name}\n")

    logger.info("writing short_names.dat")
    with open(os.path.join(outdir, "short_names.dat"), "wt") as fh:
        for e in sorted(entries, key=lambda e: e.accession):
            fh.write(f"{e.accession}\t{e.short_name}\n")

    logger.info("writing interpro2go")
    with open(os.path.join(outdir, "interpro2go"), "wt") as fh:
        fh.write(f"!date: {datetime.now():%Y/%m/%d %H:%M:%S}\n")
        fh.write("!Mapping of InterPro entries to GO\n")
        fh.write("!\n")

        for e in sorted(entries, key=lambda e: e.accession):
            for term in e.go_terms:
                fh.write(f"InterPro:{e.accession} {e.name} > "
                         f"GO:{term['name']} ; {term['identifier']}\n")

    logger.info("writing ParentChildTreeFile.txt")
    with open(os.path.join(outdir, "ParentChildTreeFile.txt"), "wt") as fh:
        for e in sorted(entries, key=lambda e: e.accession):
            root = e.hierarchy["accession"]
            if root == e.accession and e.hierarchy["children"]:
                _write_node(e.hierarchy, fh, level=0)

    logger.info("writing protein2ipr.dat.gz")
    filepath = os.path.join(outdir, "protein2ipr.dat.gz")
    with gzip.open(filepath, "wt") as fh, Store(p_uniprot2matches) as sh:
        i = 0
        for uniprot_acc, protein_entries in sh.items():
            matches = []
            for signature_acc in sorted(protein_entries):
                try:
                    interpro_acc, name = integrated[signature_acc]
                except KeyError:
                    # Not integrated signature or InterPro entry
                    continue

                locations = protein_entries[signature_acc]

                for loc in locations:
                    matches.append((
                        uniprot_acc,
                        interpro_acc,
                        name,
                        signature_acc,
                        # We do not consider fragmented locations
                        loc["fragments"][0]["start"],
                        max(f["end"] for f in loc["fragments"])
                    ))

            for m in sorted(matches):
                fh.write('\t'.join(map(str, m)) + '\n')

            i += 1
            if not i % 10000000:
                logger.info(f"{i:>12,}")

        logger.info(f"{i:>12,}")

    logger.info("complete")
