# -*- coding: utf-8 -*-

import os
from datetime import datetime
from typing import Generator

import cx_Oracle

from i7dw import logger, pdbe
from i7dw.interpro.mysql.databases import get_databases


def get_terms(cursor: cx_Oracle.Cursor) -> Generator[tuple, None, None]:
    cursor.execute(
        """
        SELECT
          I2G.ENTRY_AC, GT.GO_ID, GT.NAME, GT.CATEGORY, GC.TERM_NAME
        FROM INTERPRO.INTERPRO2GO I2G
        INNER JOIN GO.TERMS@GOAPRO GT
          ON I2G.GO_ID = GT.GO_ID
        INNER JOIN GO.CV_CATEGORIES@GOAPRO GC
          ON GT.CATEGORY = GC.CODE
        """
    )

    for row in cursor:
        yield row


def export_mappings(my_url: str, ora_url: str, outdir: str):
    databases = get_databases(my_url)
    interpro = databases["interpro"]
    version = interpro["version"]
    release_date = interpro["release_date"]

    con = cx_Oracle.connect(ora_url)
    cur = con.cursor()

    logger.info("exporting PDB-InterPro-GO-UniProt mapping")
    logger.debug("\tloading PDBe sequences from UniParc")
    cur.execute(
        """
        SELECT UPI, AC
        FROM UNIPARC.XREF
        WHERE DBID = 21
        AND DELETED = 'N'
        """
    )
    sequences = {}
    for upi, pdbe_acc in cur:
        if upi in sequences:
            sequences[upi]["structures"].add(pdbe_acc)
        else:
            sequences[upi] = {
                "structures": {pdbe_acc},
                "entries": set()
            }

    logger.debug("\tloading integrated signatures")
    # Only integrated signatures whose entry is checked and has GO terms
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

    logger.debug("\tloading GO terms in InterPro")
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

    logger.debug("\tloading PDBe matches")
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

    logger.debug("\tloading PDBe taxonomy")
    structures = pdbe.get_chain_taxonomy(cur)

    logger.debug("\tloading UniProt accessions")
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
    for pdbe_acc, protein_acc in cur:
        if not protein_acc:
            continue
        elif pdbe_acc in pdb2uniprot:
            pdb2uniprot[pdbe_acc].add(protein_acc)
        else:
            pdb2uniprot[pdbe_acc] = {protein_acc}

    os.makedirs(outdir, exist_ok=True)

    logger.debug("\twriting 'pdb2interpro2go.tsv'")
    dst = os.path.join(outdir, "pdb2interpro2go.tsv")
    with open(dst + ".tmp", "wt") as fh:
        fh.write("#PDBe ID\tchain\tTaxon ID\t"
                 "InterPro accession\tGO ID\tUniProt accession\n")

        for seq in sequences.values():
            for pdbe_acc in seq["structures"]:
                try:
                    s = structures[pdbe_acc]
                except KeyError:
                    # Structure does not exist in PDBe database
                    continue

                pdbe_id = s["id"]
                chain = s["chain"]
                proteins = pdb2uniprot.get(pdbe_acc, {''})

                for tax_id in s["taxa"]:
                    for entry_acc in seq["entries"]:
                        for go_id in entries[entry_acc]:
                            for protein_acc in proteins:
                                fh.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(
                                    pdbe_id, chain, tax_id, entry_acc,
                                    go_id, protein_acc
                                ))

    try:
        os.remove(dst)
    except FileNotFoundError:
        pass
    finally:
        os.rename(dst + ".tmp", dst)
        os.chmod(dst, 0o777)

    logger.info("exporting InterPro-GO-UniProt mapping")
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

    dst = os.path.join(outdir, "interpro2go2uniprot.tsv")
    with open(dst + ".tmp", "wt") as fh:
        fh.write("#InterPro accession\tGO ID\tUniProt accession\n")

        for row in cur:
            fh.write('\t'.join(row) + '\n')

    try:
        os.remove(dst)
    except FileNotFoundError:
        pass
    finally:
        os.rename(dst + ".tmp", dst)
        os.chmod(dst, 0o777)

    cur.close()
    con.close()

    dst = os.path.join(outdir, "release.txt")
    with open(dst, "wt") as fh:
        fh.write("InterPro version:        "
                 "{}\n".format(version))

        fh.write("Release date:            "
                 "{:%A, %d %B %Y}\n".format(release_date))

        fh.write("Generated on:            "
                 "{:%Y-%m-%d %H:%M}\n".format(datetime.now()))
    os.chmod(dst, 0o777)

    logger.info("complete")
