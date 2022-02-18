import os
import shutil
from datetime import datetime
from typing import Dict

import cx_Oracle
import MySQLdb

from interpro7dw import pdbe
from interpro7dw.utils.mysql import url2dict
from interpro7dw.utils.store import SimpleStore, loadobj


_PDB2INTERPRO2GO2 = "pdb2interpro2go.tsv"
_INTERPRO2GO2UNIPROT = "interpro2go2uniprot.tsv"


def get_terms(uri: str) -> Dict[str, tuple]:
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


def export(ipr_pro_uri: str, ipr_stg_uri: str, pdbe_uri: str,
           entries_file: str, entry2xrefs_file: str, outdir: str):
    os.makedirs(outdir, exist_ok=True)

    file = os.path.join(outdir, _PDB2INTERPRO2GO2)
    _export_pdb2ipr2go(ipr_pro_uri, pdbe_uri, entries_file, file)

    file = os.path.join(outdir, _PDB2INTERPRO2GO2)
    _export_ipr2go2uniprot(entries_file, entry2xrefs_file, file)

    con = MySQLdb.connect(**url2dict(ipr_stg_uri), charset="utf8mb4")
    cur = con.cursor()
    cur.execute(
        """
        SELECT version, release_date 
        FROM webfront_database 
        WHERE name='interpro'
        """
    )
    version, date = cur.fetchone()
    cur.close()
    con.close()

    file = os.path.join(outdir, "release.txt")
    with open(file, "wt") as fh:
        fh.write(f"InterPro version:    {version}\n")
        fh.write(f"Release date:        {date:%A, %d %B %Y}\n")
        fh.write(f"Generated on:        {datetime.now():%Y-%m-%d %H:%M}\n")

    os.chmod(file, 0o775)


def _export_pdb2ipr2go(ipr_pro_uri: str, pdbe_uri: str, entries_file: str,
                       output: str):
    entries = loadobj(entries_file)

    con = cx_Oracle.connect(ipr_pro_uri)
    cur = con.cursor()

    # PDBe sequences from UniParc
    cur.execute(
        """
        SELECT UPI, AC
        FROM UNIPARC.XREF
        WHERE DBID = 21
        AND DELETED = 'N'
        """
    )
    sequences = {}
    for upi, pdb_acc in cur:
        if upi in sequences:
            sequences[upi]["structures"].add(pdb_acc)
        else:
            sequences[upi] = {
                "structures": {pdb_acc},
                "entries": set()
            }

    # PDBe matches
    cur.execute(
        """
        SELECT DISTINCT UPI, METHOD_AC
        FROM IPRSCAN.MV_IPRSCAN
        WHERE UPI IN (
            SELECT UPI
            FROM UNIPARC.XREF
            WHERE DBID = 21
            AND DELETED = 'N'
        )
        """
    )
    for upi, signature_acc in cur:
        try:
            entry = entries[signature_acc]
        except KeyError:
            continue

        if not entry.integrated_in:
            # Non-integrated signature
            continue

        interpro_entry = entries[entry.integrated_in]
        if interpro_entry.go_terms:
            sequences[upi]["entries"].add(interpro_entry.accession)

    # UniProt accessions
    cur.execute(
        """
        SELECT DISTINCT A.AC, B.AC
        FROM UNIPARC.XREF A
        LEFT OUTER JOIN UNIPARC.XREF B ON A.UPI = B.UPI
        WHERE A.DBID = 21
        AND A.DELETED = 'N'
        AND B.DBID IN (2, 3)
        AND B.DELETED = 'N'
        """
    )
    pdb2uniprot = {}
    for pdb_acc, protein_acc in cur:
        if not protein_acc:
            continue
        elif pdb_acc in pdb2uniprot:
            pdb2uniprot[pdb_acc].add(protein_acc)
        else:
            pdb2uniprot[pdb_acc] = {protein_acc}

    cur.close()
    con.close()

    # PDBe taxonomy
    structures = pdbe.get_chain_taxonomy(pdbe_uri)

    with open(output, "wt") as fh:
        fh.write("#PDBe ID\tchain\tTaxon ID\t"
                 "InterPro accession\tGO ID\tUniProt accession\n")

        for seq in sequences.values():
            for pdb_acc in seq["structures"]:
                try:
                    s = structures[pdb_acc]
                except KeyError:
                    # Structure does not exist in PDBe database
                    continue

                pdb_id = s["id"]
                chain = s["chain"]
                proteins = pdb2uniprot.get(pdb_acc, {''})

                for tax_id in s["taxa"]:
                    for entry_acc in seq["entries"]:
                        for term in entries[entry_acc].go_terms:
                            go_id = term["identifier"]

                            for protein_acc in proteins:
                                fh.write(f"{pdb_id}\t{chain}\t{tax_id}\t"
                                         f"{entry_acc}\t{go_id}\t"
                                         f"{protein_acc}\n")

    os.chmod(output, 0o775)


def _export_ipr2go2uniprot(entries_file: str, entry2xrefs_file: str,
                           output: str):
    entries = loadobj(entries_file)
    with SimpleStore(entry2xrefs_file) as store, open(output, "wt") as fh:
        fh.write("#InterPro accession\tGO ID\tUniProt accession\n")

        for accession, entry_xrefs in store:
            entry = entries[accession]
            if entry.database == "interpro":
                for term in entry.go_terms:
                    go_id = term["identifier"]

                    for uniprot_acc, uniprot_id in entry_xrefs["proteins"]:
                        fh.write(f"{accession}\t{go_id}\t{uniprot_acc}\n")

    os.chmod(output, 0o775)


def publish(src: str, dst: str):
    os.makedirs(dst, exist_ok=True)

    for name in os.listdir(src):
        path = os.path.join(dst, name)
        try:
            os.unlink(path)
        except FileNotFoundError:
            pass
        finally:
            shutil.copy(os.path.join(src, name), path)
