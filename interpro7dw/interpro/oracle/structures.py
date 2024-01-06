import glob
import gzip
import os
import shelve

import oracledb

from interpro7dw.pdbe import get_sequences
from interpro7dw.utils import logger
from interpro7dw.utils.oracle import lob_as_str
from .databases import get_databases_codes
from .entries import load_entries, load_signatures
from .matches import get_fragments, merge_uniprot_matches


def update_pdbe_matches(uri: str):
    logger.info("starting")
    con = oracledb.connect(uri)
    cur = con.cursor()
    cur.execute("TRUNCATE TABLE IPRSCAN.MV_PDB_MATCH REUSE STORAGE")
    cur.execute(
        """
        INSERT INTO IPRSCAN.MV_PDB_MATCH
        SELECT X.AC, E.ENTRY_AC, EM.METHOD_AC, M.SEQ_START, M.SEQ_END
        FROM UNIPARC.XREF X
        INNER JOIN IPRSCAN.MV_IPRSCAN M
          ON X.UPI = M.UPI
        INNER JOIN INTERPRO.ENTRY2METHOD EM 
          ON M.METHOD_AC = EM.METHOD_AC
        INNER JOIN INTERPRO.ENTRY E 
          ON (EM.ENTRY_AC = E.ENTRY_AC AND E.CHECKED = 'Y')
        WHERE X.DBID = 21
          AND X.DELETED = 'N'
        """
    )
    cnt = cur.rowcount
    con.commit()
    cur.close()
    con.close()
    logger.info(f"{cnt:,} rows inserted")


def export_matches(ipr_uri: str, pdbe_uri: str, output: str):
    logger.info("starting")
    for file in glob.glob(f"{output}*"):
        os.unlink(file)

    con = oracledb.connect(ipr_uri)
    cur = con.cursor()
    cur.outputtypehandler = lob_as_str

    with shelve.open(output, writeback=True) as db:
        logger.info("exporting sequences from UniParc")
        cur.execute(
            """
            SELECT X.AC, P.SEQ_SHORT, P.SEQ_LONG
            FROM UNIPARC.XREF X
            INNER JOIN UNIPARC.PROTEIN P ON X.UPI = P.UPI
            WHERE X.DBID = 21 
              AND X.DELETED = 'N' 
            """
        )

        i = 0
        chains = set()
        for pdb_chain, seq_short, seq_long in cur:
            sequence = seq_short or seq_long
            chains.add(pdb_chain)
            db[pdb_chain] = {
                "length": len(sequence),
                "sequence": gzip.compress(sequence.encode("utf-8")),
                "matches": []
            }

            i += 1
            if i % 1000 == 0:
                db.sync()

        db.sync()

        logger.info("exporting sequences from PDBe")
        for pdb_chain, sequence in get_sequences(pdbe_uri):
            if pdb_chain not in chains:
                db[pdb_chain] = {
                    "length": len(sequence),
                    "sequence": gzip.compress(sequence.encode("utf-8")),
                    "matches": []
                }

                i += 1
                if i % 1000 == 0:
                    db.sync()

        db.sync()

        logger.info("exporting matches")
        dbcodes, _ = get_databases_codes(cur)
        params = [f":{i+1}" for i in range(len(dbcodes))]
        cur.execute(
            f"""
            WITH ANALYSES AS (
                SELECT IPRSCAN_SIG_LIB_REL_ID AS ID
                FROM INTERPRO.IPRSCAN2DBCODE
                WHERE DBCODE IN ({','.join(params)})
            )
            SELECT X.AC, M.METHOD_AC, M.SEQ_START, M.SEQ_END, M.FRAGMENTS
            FROM UNIPARC.XREF X
            INNER JOIN IPRSCAN.MV_IPRSCAN M
                ON X.UPI = M.UPI
            WHERE X.DBID = 21 
              AND X.DELETED = 'N'
              AND M.ANALYSIS_ID IN (SELECT ID FROM ANALYSES)      
            """,
            dbcodes
        )

        i = 0
        for pdb_chain, signature_acc, pos_start, pos_end, fragments in cur:
            db[pdb_chain]["matches"].append((
                signature_acc,
                None,  # model accession
                None,  # feature
                None,  # score
                get_fragments(pos_start, pos_end, fragments)
            ))

            i += 1
            if i % 1e3 == 0:
                db.sync()

            if i % 1e6 == 0:
                logger.info(f"{i:>15,}")

        db.sync()
        logger.info(f"{i:>15,}")

    logger.info("loading entries")
    entries = load_entries(cur)
    signatures = load_signatures(cur)

    cur.close()
    con.close()

    logger.info("merging matches")
    with shelve.open(output, writeback=False) as db:
        for pdb_chain, obj in db.items():
            s, e = merge_uniprot_matches(obj["matches"], signatures, entries)
            obj["matches"] = {**s, **e}
            db[pdb_chain] = obj

    logger.info("done")
