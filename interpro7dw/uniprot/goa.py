import os
import pickle
import shelve
from datetime import datetime

import cx_Oracle

from interpro7dw.utils.store import BasicStore, KVStore, copy_files


_PDB2INTERPRO2GO2 = "pdb2interpro2go.tsv"
_INTERPRO2GO2UNIPROT = "interpro2go2uniprot.tsv"
_TREEGRAFTER2GO2UNIPROT = "treegrafter2go2uniprot.tsv"


def get_terms(uri: str) -> dict[str, tuple]:
    con = cx_Oracle.connect(uri)
    cur = con.cursor()
    cur.execute(
        """
        SELECT GT.GO_ID, GT.NAME, GT.CATEGORY, GC.TERM_NAME, GC.SORT_ORDER
        FROM GO.TERMS GT
        INNER JOIN GO.CV_CATEGORIES GC
          ON GT.CATEGORY = GC.CODE
        """
    )
    terms = {row[0]: row[1:] for row in cur}
    cur.close()
    con.close()
    return terms


def export(databases_file: str, entries_file: str, matches_file:str,
           structures_file: str, pdb2matches_file: str, uniprot2pdb_file: str,
           entry2xrefs_file: str, outdir: str):
    os.makedirs(outdir, exist_ok=True)

    with open(entries_file, "rb") as fh:
        entries = pickle.load(fh)

    _export_ipr2go2uni(entries,
                       entry2xrefs_file,
                       os.path.join(outdir, _INTERPRO2GO2UNIPROT))

    _export_pthr2go2uni(entries,
                        matches_file,
                        os.path.join(outdir, _TREEGRAFTER2GO2UNIPROT))

    _export_pdb2ipr2go(entries,
                       structures_file,
                       pdb2matches_file,
                       uniprot2pdb_file,
                       os.path.join(outdir, _PDB2INTERPRO2GO2))

    release_version = release_date = None
    with open(databases_file, "rb") as fh:
        for db in pickle.load(fh).values():
            if db["name"].lower() == "interpro":
                release_version = db["release"]["version"]
                release_date = db["release"]["date"]
                break

    if release_version is None:
        raise RuntimeError("missing release version/date for InterPro")

    file = os.path.join(outdir, "release.txt")
    with open(file, "wt") as fh:
        fh.write(f"InterPro version:    {release_version}\n")
        fh.write(f"Release date:        {release_date:%A, %d %B %Y}\n")
        fh.write(f"Generated on:        {datetime.now():%Y-%m-%d %H:%M}\n")


def _export_pdb2ipr2go(entries: dict, structures_file: str,
                       pdb2matches_file: str, uniprot2pdb_file: str,
                       output: str):
    with open(structures_file, "rb") as fh:
        structures = pickle.load(fh)

    pdb2uniprot = {}
    with open(uniprot2pdb_file, "rb") as fh:
        for protein_acc, pdb_entries in pickle.load(fh).items():
            for pdb_chain in pdb_entries:
                try:
                    pdb2uniprot[pdb_chain].add(protein_acc)
                except KeyError:
                    pdb2uniprot[pdb_chain] = {protein_acc}

    with (shelve.open(pdb2matches_file, writeback=False) as d,
          open(output, "wt") as fh):
        fh.write("#PDBe ID\tchain\tTaxon ID\t"
                 "InterPro accession\tGO ID\tUniProt accession\n")

        for pdb_chain, pdb_entry in d.items():
            pdb_id, chain = pdb_chain.split("_")

            try:
                structure = structures[pdb_id]
            except KeyError:
                continue

            taxon_id = structure["taxonomy"].get(chain)
            if not taxon_id:
                continue

            # If not proteins: use empty field
            proteins = pdb2uniprot.get(pdb_chain, [""])

            for entry_acc in pdb_entry["matches"]:
                entry = entries[entry_acc]

                if not entry.public:
                    continue

                for term in entry.go_terms:
                    go_id = term["identifier"]

                    for protein_acc in proteins:
                        fh.write(f"{pdb_id}\t{chain}\t"
                                 f"{taxon_id}\t{entry_acc}\t"
                                 f"{go_id}\t{protein_acc}\n")


def _export_ipr2go2uni(entries: dict, xrefs_file: str,
                       output: str = _INTERPRO2GO2UNIPROT):
    with BasicStore(xrefs_file, mode="r") as sh, open(output, "wt") as fh:
        fh.write("#InterPro accession\tGO ID\tUniProt accession\n")

        for accession, entry_xrefs in sh:
            entry = entries[accession]

            if entry.database.lower() == "interpro" and entry.public:
                for term in entry.go_terms:
                    go_id = term["identifier"]

                    for uniprot_acc, _, _ in entry_xrefs["proteins"]:
                        fh.write(f"{accession}\t{go_id}\t{uniprot_acc}\n")


def _export_pthr2go2uni(entries: dict, matches_file: str,
                        output: str = _TREEGRAFTER2GO2UNIPROT):
    with KVStore(matches_file) as kvs, open(output, "wt") as fh:
        fh.write("#PANTHER accession\tInterPro accession\t"
                 "GO ID\tUniProt accession\n")

        for protein_acc, (signatures, _) in kvs.items():
            subfamilies = set()
            for signature_acc, match in signatures.items():
                if match["database"].lower() == "panther":
                    entry_acc = match["entry"] or "-"

                    for loc in match["locations"]:
                        try:
                            subfam = loc["subfamily"]
                        except KeyError:
                            # PANTHER match, but no subfamily
                            # (no SF associated to the grated node)
                            continue
                        else:
                            subfamilies.add((subfam["accession"], entry_acc))

            for subfam_acc, interpro_acc in subfamilies:
                subfam = entries[subfam_acc]

                for term in subfam.go_terms:
                    go_id = term["identifier"]

                    fh.write(f"{subfam_acc}\t{interpro_acc}\t{go_id}\t"
                             f"{protein_acc}\n")


def publish(src: str, dst: str):
    copy_files(src, dst)
