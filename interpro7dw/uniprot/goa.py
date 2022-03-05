import os
import pickle
import shutil
from datetime import datetime

import cx_Oracle

from interpro7dw.utils.store import BasicStore


_PDB2INTERPRO2GO2 = "pdb2interpro2go.tsv"
_INTERPRO2GO2UNIPROT = "interpro2go2uniprot.tsv"


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


def export(databases_file: str, entries_file: str, structures_file: str,
           pdb_matches_file: str, entry2xrefs_file: str, outdir: str):
    os.makedirs(outdir, exist_ok=True)

    with open(entries_file, "rb") as fh:
        entries = pickle.load(fh)

    file = os.path.join(outdir, _INTERPRO2GO2UNIPROT)
    _export_ipr2go2uni(entries, entry2xrefs_file, file)

    file = os.path.join(outdir, _PDB2INTERPRO2GO2)
    _export_pdb2ipr2go(entries, structures_file, pdb_matches_file, file)

    release_version = release_date = None
    with open(databases_file, "rb") as fh:
        for db in pickle.load(fh).values():
            if db["name"] == "interpro":
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

    # TODO: review if necessary
    # os.chmod(output, 0o774)


def _export_pdb2ipr2go(entries: dict, structures_file: str, matches_file: str,
                       output: str):
    with open(structures_file, "rb") as fh:
        pdb2taxonomy = pickle.load(fh)["taxonomy"]

    with BasicStore(matches_file, mode="r") as sh, open(output, "wt") as fh:
        fh.write("#PDBe ID\tchain\tTaxon ID\t"
                 "InterPro accession\tGO ID\tUniProt accession\n")

        for pdb_key, matches, proteins in sh:
            # pdb_key: PDB ID + '_' + chain ID
            try:
                structure = pdb2taxonomy[pdb_key]
            except KeyError:
                continue

            pdb_id = structure["id"]
            chain_id = structure["chain"]

            for taxon_id in structure["taxa"]:
                for signature_acc, entry_acc, pos_start, pos_end in matches:
                    for term in entries[entry_acc].go_terms:
                        go_id = term["identifier"]

                        # If no UniProt proteins: use empty field
                        for protein_acc in proteins or [""]:
                            fh.write(f"{pdb_id}\t{chain_id}\t"
                                     f"{taxon_id}\t{entry_acc}\t"
                                     f"{go_id}\t{protein_acc}\n")

    # TODO: review if necessary
    # os.chmod(output, 0o774)


def _export_ipr2go2uni(entries: dict, xrefs_file: str, output: str):
    with BasicStore(xrefs_file, mode="r") as sh, open(output, "wt") as fh:
        fh.write("#InterPro accession\tGO ID\tUniProt accession\n")

        for accession, entry_xrefs in sh:
            entry = entries[accession]
            if entry.database.lower() == "interpro":
                for term in entry.go_terms:
                    go_id = term["identifier"]

                    for uniprot_acc, uniprot_id in entry_xrefs["proteins"]:
                        fh.write(f"{accession}\t{go_id}\t{uniprot_acc}\n")

    # TODO: review if necessary
    # os.chmod(output, 0o774)


def publish(src: str, dst: str):
    os.umask(0o002)
    os.makedirs(dst, mode=0o775, exist_ok=True)

    for name in os.listdir(src):
        path = os.path.join(dst, name)
        try:
            os.unlink(path)
        except FileNotFoundError:
            pass
        finally:
            shutil.copy(os.path.join(src, name), path)
