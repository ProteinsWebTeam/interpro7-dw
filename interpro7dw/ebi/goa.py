# -*- coding: utf-8 -*-

import os
from datetime import datetime

import cx_Oracle
import MySQLdb

from interpro7dw import logger
from interpro7dw.ebi import pdbe
from interpro7dw.utils import url2dict


def _export_pdb2interpro2go2uniprot(cur: cx_Oracle.Cursor, output: str):
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

    # Integrated signatures whose entry is checked and has GO terms
    cur.execute(
        """
        SELECT METHOD_AC, ENTRY_AC
        FROM INTERPRO.ENTRY2METHOD
        WHERE ENTRY_AC IN (
          SELECT ENTRY_AC
          FROM INTERPRO.ENTRY
          WHERE CHECKED = 'Y'
          INTERSECT
          SELECT DISTINCT ENTRY_AC
          FROM INTERPRO.INTERPRO2GO
        )
        """
    )
    signatures = dict(cur.fetchall())

    # GO terms in InterPro
    cur.execute(
        """
        SELECT ENTRY_AC, GO_ID
        FROM INTERPRO.INTERPRO2GO
        WHERE ENTRY_AC IN (
          SELECT ENTRY_AC
          FROM INTERPRO.ENTRY
          WHERE CHECKED = 'Y'
        )
        """
    )
    entries = {}
    for entry_acc, go_id in cur:
        if entry_acc in entries:
            entries[entry_acc].add(go_id)
        else:
            entries[entry_acc] = {go_id}

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
            entry_acc = signatures[signature_acc]
        except KeyError:
            pass
        else:
            sequences[upi]["entries"].add(entry_acc)

    # PDBe taxonomy
    structures = pdbe.get_chain_taxonomy(cur)

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

    tmp_path = f"{output}.tmp"
    with open(tmp_path, "wt") as fh:
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
                        for go_id in entries[entry_acc]:
                            for protein_acc in proteins:
                                fh.write(f"{pdb_id}\t{chain}\t{tax_id}\t"
                                         f"{entry_acc}\t{go_id}\t"
                                         f"{protein_acc}\n")

    try:
        os.remove(output)
    except FileNotFoundError:
        pass
    finally:
        os.rename(tmp_path, output)
        os.chmod(output, 0o775)


def _export_interpro2go2uniprot(cur: cx_Oracle.Cursor, output: str):
    cur.execute(
        """
        SELECT DISTINCT IG.ENTRY_AC, IG.GO_ID, M.PROTEIN_AC
        FROM INTERPRO.MATCH M
        INNER JOIN INTERPRO.ENTRY2METHOD EM
          ON M.METHOD_AC = EM.METHOD_AC
        INNER JOIN INTERPRO.ENTRY E
          ON EM.ENTRY_AC = E.ENTRY_AC
        INNER JOIN INTERPRO.INTERPRO2GO IG
          ON E.ENTRY_AC = IG.ENTRY_AC
        WHERE E.CHECKED = 'Y'
        """
    )

    tmp_path = f"{output}.tmp"
    with open(tmp_path, "wt") as fh:
        fh.write("#InterPro accession\tGO ID\tUniProt accession\n")

        for row in cur:
            fh.write('\t'.join(row) + '\n')

    try:
        os.remove(output)
    except FileNotFoundError:
        pass
    finally:
        os.rename(tmp_path, output)
        os.chmod(output, 0o775)


def export(pro_url: str, stg_url: str, outdir: str):
    os.makedirs(outdir, exist_ok=True)

    con = cx_Oracle.connect(pro_url)
    cur = con.cursor()

    logger.info("exporting PDB-InterPro-GO-UniProt mapping")
    filepath = os.path.join(outdir, "pdb2interpro2go.tsv")
    _export_pdb2interpro2go2uniprot(cur, filepath)

    logger.info("exporting InterPro-GO-UniProt mapping")
    filepath = os.path.join(outdir, "interpro2go2uniprot.tsv")
    _export_interpro2go2uniprot(cur, filepath)
    cur.close()
    con.close()

    logger.info("exporting release info")
    con = MySQLdb.connect(**url2dict(stg_url))
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

    filepath = os.path.join(outdir, "release.txt")
    with open(filepath, "wt") as fh:
        fh.write(f"InterPro version:    {version}\n")
        fh.write(f"Release date:        {date:%A, %d %B %Y}\n")
        fh.write(f"Generated on:        {datetime.now():%Y-%m-%d %H:%M}\n")

    os.chmod(filepath, 0o775)
    logger.info("complete")
