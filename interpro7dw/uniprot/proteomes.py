import pickle

import oracledb


def export_proteomes(uri: str, file: str):
    con = oracledb.connect(uri)
    cur = con.cursor()
    cur.execute(
        """
        SELECT P.UPID, P.PROTEOME_NAME, P.IS_REFERENCE, P.GC_SET_ACC, 
          TO_CHAR(P.PROTEOME_TAXID), SN.NAME
        FROM SPTR.PROTEOME P
        LEFT OUTER JOIN TAXONOMY.SPTR_STRAIN S
          ON P.PROTEOME_TAXID = S.TAX_ID
        LEFT OUTER JOIN TAXONOMY.SPTR_STRAIN_NAME SN
          ON S.STRAIN_ID = SN.STRAIN_ID
        WHERE P.IS_REFERENCE = 1
        """
    )

    proteomes = {}
    for upid, name, is_reference, assembly, taxon_id, strain in cur:
        if upid in proteomes:
            continue

        proteomes[upid] = {
            "id": upid,
            "name": name,
            "is_reference": is_reference != 0,
            "assembly": assembly,
            "taxon_id": taxon_id,
            "strain": strain
        }

    cur.close()
    con.close()

    with open(file, "wb") as fh:
        pickle.dump(proteomes, fh)
