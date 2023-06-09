from interpro7dw.utils import logger
from interpro7dw.utils.oracle import lob_as_str
from interpro7dw.utils.store import BasicStore

import oracledb


def export_rosettafold(uri: str, output: str, raise_on_empty: bool = True):
    logger.info("exporting structural models")
    con = oracledb.connect(uri)
    cur = con.cursor()
    cur.execute(
        """
        SELECT EM.METHOD_AC, E.ENTRY_AC
        FROM INTERPRO.ENTRY2METHOD EM
        INNER JOIN INTERPRO.ENTRY E
            ON EM.ENTRY_AC = E.ENTRY_AC
        WHERE E.CHECKED = 'Y'
        """
    )
    integrated = dict(cur.fetchall())

    cur.outputtypehandler = lob_as_str
    cur.execute(
        """
        SELECT METHOD_AC, CONTACTS, PLDDT, STRUCTURE
        FROM INTERPRO.ROSETTAFOLD
        """
    )

    n = 0
    with BasicStore(output, mode="w") as store:
        for sig_acc, cmap_gz, plddt_gz, pdb_gz in cur:
            store.write((sig_acc, integrated.get(sig_acc), cmap_gz, plddt_gz,
                         pdb_gz))
            n += 1

    cur.close()
    con.close()

    if n > 0 or not raise_on_empty:
        logger.info(f"done: {n} models exported")
    else:
        raise RuntimeError("No model exported")


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


def export_pdbe_matches(uri: str, output: str):
    update_pdbe_matches(uri)

    con = oracledb.connect(uri)
    cur = con.cursor()

    logger.info("loading integrated signatures")
    cur.execute(
        """
        SELECT METHOD_AC, ENTRY_AC
        FROM INTERPRO.ENTRY2METHOD
        WHERE ENTRY_AC IN (
            SELECT ENTRY_AC
            FROM INTERPRO.ENTRY
            WHERE CHECKED = 'Y'
        )
        """
    )
    integrated = dict(cur.fetchall())

    logger.info("loading UniProt-PDBe mapping")
    cur.execute(
        """
        SELECT DISTINCT A.AC, B.AC
        FROM UNIPARC.XREF A
        INNER JOIN UNIPARC.XREF B ON A.UPI = B.UPI
        WHERE A.DBID = 21
          AND A.DELETED = 'N'
          AND B.DBID IN (2, 3)
          AND B.DELETED = 'N'
        """
    )

    pdbe2uniprot = {}
    for pdbe_acc, uniprot_acc in cur:
        try:
            pdbe2uniprot[pdbe_acc].append(uniprot_acc)
        except KeyError:
            pdbe2uniprot[pdbe_acc] = [uniprot_acc]

    with BasicStore(output, mode="w") as store:
        logger.info("exporting PDBe matches")
        cur.execute(
            """
            SELECT DISTINCT X.AC, M.METHOD_AC, M.SEQ_START, M.SEQ_END
            FROM UNIPARC.XREF X
            INNER JOIN IPRSCAN.MV_IPRSCAN M 
                ON X.UPI = M.UPI
            WHERE X.DBID = 21 AND X.DELETED = 'N'
            ORDER BY X.AC
            """
        )

        pdbe_id = None
        matches = []
        for _pdbe_id, signature_acc, pos_start, pos_end in cur:
            try:
                entry_acc = integrated[signature_acc]
            except KeyError:
                continue

            if _pdbe_id != pdbe_id:
                if pdbe_id:
                    store.write((pdbe_id,
                                 matches,
                                 pdbe2uniprot.get(pdbe_id, [])))

                pdbe_id = _pdbe_id
                matches = []

            matches.append((signature_acc, entry_acc, pos_start, pos_end))

        if pdbe_id:
            store.write((pdbe_id,
                         matches,
                         pdbe2uniprot.get(pdbe_id, [])))

    cur.close()
    con.close()

    logger.info("done")
