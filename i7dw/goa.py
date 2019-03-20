import logging
import os
from datetime import datetime
from typing import Generator

from cx_Oracle import Cursor

from . import dbms
from .interpro import mysql

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s: %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)


def get_terms(cursor: Cursor) -> Generator[tuple, None, None]:
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


def export_mapping_files(my_url: str, ora_url: str, outdir: str):
    databases = mysql.get_entry_databases(my_url)
    interpro = databases["interpro"]
    version = interpro["version"]
    release_date = interpro["release_date"]

    logging.info("exporting PDB-InterPro-GO[-UniProt] mapping")
    con, cur = dbms.connect(ora_url)
    cur.execute(
        """
        SELECT DISTINCT 
          UPPER(M.AC), 
          ASYM.ENTRY_ID, 
          ASYM.AUTH_ASYM_ID, 
          SRC.TAX_ID,
          E.ENTRY_AC,
          E.GO_ID,
          X.AC
        FROM (
            /* Select signature matches against PDBe sequences */
            SELECT DISTINCT M.METHOD_AC, M.UPI, X.AC
            FROM IPRSCAN.MV_IPRSCAN M
            INNER JOIN UNIPARC.XREF X
              ON M.UPI = X.UPI AND X.DBID = 21 AND X.DELETED = 'N'        
        ) M
        INNER JOIN (
          /* Select GO terms associated to integrated signatures */
          SELECT EM.METHOD_AC, E.ENTRY_AC, EG.GO_ID
          FROM INTERPRO.ENTRY E
          INNER JOIN INTERPRO.ENTRY2METHOD EM 
            ON E.ENTRY_AC = EM.ENTRY_AC
          INNER JOIN INTERPRO.INTERPRO2GO EG 
            ON E.ENTRY_AC = EG.ENTRY_AC
          WHERE E.CHECKED = 'Y'
        ) E ON M.METHOD_AC = E.METHOD_AC
        /* Select PDB ID and chain */
        INNER JOIN PDBE.STRUCT_ASYM@PDBE_LIVE ASYM
          ON M.AC = ASYM.ENTRY_ID || '_' || ASYM.AUTH_ASYM_ID
        /* Taxonomy identifier of species corresponding to PDB entry  */
        INNER JOIN PDBE.ENTITY_SRC@PDBE_LIVE SRC
          ON ASYM.ENTRY_ID = SRC.ENTRY_ID AND ASYM.ENTITY_ID = SRC.ENTITY_ID
        /* Add UniProt accession if 100% sequence similarity with PDB chain matches */
        LEFT OUTER JOIN UNIPARC.XREF X
          ON M.UPI = X.UPI AND X.DBID IN (2, 3) AND X.DELETED = 'N'


                  
        """
    )

    dst = os.path.join(outdir, "pdb2interpro2go.tsv")
    with open(dst + ".tmp", "wt") as fh:
        fh.write("#PDBe ID\tPDBe accession\tPDBe chain\tTaxon ID\t"
                 "InterPro accession\tGO ID\tUniProt accession\n")

        for row in cur:
            fh.write('\t'.join(map(str, row)) + '\n')

    try:
        os.remove(dst)
    except FileNotFoundError:
        pass
    finally:
        os.rename(dst + ".tmp", dst)
        os.chmod(dst, 0o777)

    logging.info("exporting InterPro-GO-UniProt mapping")
    cur.execute(
        """
        SELECT DISTINCT IG.ENTRY_AC, IG.GO_ID, M.PROTEIN_AC 
        FROM INTERPRO.INTERPRO2GO IG
        INNER JOIN INTERPRO.ENTRY2METHOD EM 
          ON IG.ENTRY_AC = EM.ENTRY_AC
        INNER JOIN INTERPRO.MATCH M 
          ON EM.METHOD_AC = M.METHOD_AC
        WHERE IG.ENTRY_AC IN (
          SELECT ENTRY_AC FROM INTERPRO.ENTRY WHERE CHECKED = 'Y'
        )
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

    with open(os.path.join(outdir, "release.txt"), "wt") as fh:
        fh.write("InterPro version:        "
                 "{}\n".format(version))

        fh.write("Release date:            "
                 "{:%Y-%m-%d:%H:%M}\n".format(release_date))

        fh.write("Generated on:            "
                 "{:%Y-%m-%d:%H:%M}\n".format(datetime.now()))

    logging.info("complete")
