import cx_Oracle
from interpro7dw.utils.oracle import lob_as_str


def export(uri: str, output: str):
    con = cx_Oracle.connect(uri)
    cur = con.cursor()
    cur.execute(
        """
        SELECT E.ENTRY_AC, E.NAME, ET.ABBREV
        FROM INTERPRO.ENTRY E
        INNER JOIN INTERPRO.CV_ENTRY_TYPE ET
          ON E.ENTRY_TYPE = ET.CODE
        WHERE E.CHECKED = 'Y'
        """
    )
    entries = {rec[0]: rec[1:] for rec in cur}

    cur.execute(
        """
        SELECT M.METHOD_AC, M.DESCRIPTION, D.DBSHORT, ET.ABBREV, I2D.EVIDENCE, 
               EM.ENTRY_AC
        FROM INTERPRO.METHOD M
        INNER JOIN INTERPRO.CV_DATABASE D
          ON M.DBCODE = D.DBCODE
        INNER JOIN INTERPRO.CV_ENTRY_TYPE ET
          ON M.SIG_TYPE = ET.CODE
        INNER JOIN INTERPRO.IPRSCAN2DBCODE I2D 
          ON M.DBCODE = I2D.DBCODE          
        LEFT OUTER JOIN INTERPRO.ENTRY2METHOD EM
          ON M.METHOD_AC = EM.METHOD_AC
        """
    )
    signatures = {rec[0]: rec[1:] for rec in cur}

    cur.outputtypehandler = lob_as_str
    cur.execute(
        """
        SELECT V.PROTEIN_AC, V.VARIANT, V.LENGTH, V.CRC64, 
               P.SEQ_SHORT, P.SEQ_LONG
        FROM INTERPRO.VARSPLIC_MASTER V
        INNER JOIN UNIPARC.PROTEIN P ON V.CRC64 = P.CRC64
        """
    )

    isoforms = {}
    for rec in cur:
        variant_acc = rec[0] + '-' + str(rec[1])
        isoforms[variant_acc] = {
            "accession": variant_acc,
            "protein": rec[0],
            "length": rec[2],
            "crc64": rec[3],
            "sequence": rec[4] or rec[5],
            "matches": []
        }

    # PROTEIN_AC is actually PROTEIN-VARIANT (e.g. Q13733-1)
    cur.execute(
        """
        SELECT V.PROTEIN_AC, V.METHOD_AC, E.ENTRY_AC, V.MODEL_AC, V.SCORE, 
               V.POS_FROM, V.POS_TO, V.FRAGMENTS
        FROM INTERPRO.VARSPLIC_MATCH V
        INNER JOIN INTERPRO.METHOD M ON V.METHOD_AC = M.METHOD_AC
        LEFT OUTER JOIN (
            SELECT EM.METHOD_AC, EM.ENTRY_AC 
            FROM INTERPRO.ENTRY2METHOD EM
            INNER JOIN INTERPRO.ENTRY E 
                ON EM.ENTRY_AC = E.ENTRY_AC
            WHERE E.CHECKED = 'Y'        
        ) E ON M.METHOD_AC = E.METHOD_AC
        """
    )

    for rec in cur:
        variant_acc = rec[0]
        signature_acc = rec[1]
        entry_acc = rec[2]
        model_acc = rec[3]
        score = rec[4]
        pos_start = rec[5]
        pos_end = rec[6]
        fragments = rec[7]

        try:
            isoform = isoforms[variant_acc]
        except KeyError:
            continue

        if fragments:
            _fragments = []

            for frag in fragments.split(','):
                # Format: START-END-STATUS
                s, e, t = frag.split('-')
                _fragments.append({
                    "start": int(s),
                    "end": int(e),
                    "dc-status": DC_STATUSES[t]
                })

            fragments = _fragments
        else:
            fragments = [{
                "start": pos_start,
                "end": pos_end,
                "dc-status": DC_STATUSES['S']  # Continuous
            }]

        isoform["matches"].append((signature_acc, model_acc, score,
                                   fragments, entry_acc))

    cur.close()
    con.close()

    with SimpleStore(dst) as store:
        for isoform in isoforms.values():
            isoform["matches"] = _merge_matches(isoform["matches"])
            store.add(isoform)
