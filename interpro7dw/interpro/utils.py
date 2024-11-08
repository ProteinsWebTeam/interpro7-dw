import copy


def copy_dict(src: dict, dst: dict, concat_or_incr: bool = False):
    for key, value in src.items():
        if key in dst:
            if isinstance(value, dict):
                copy_dict(value, dst[key], concat_or_incr)
            elif isinstance(value, (list, tuple)):
                dst[key] += value
            elif isinstance(value, set):
                dst[key] |= value
            elif isinstance(value, (int, float, str)) and concat_or_incr:
                dst[key] += value
            else:
                dst[key] = value
        else:
            dst[key] = copy.deepcopy(value)


def overlaps_pdb_chain(locations: list[dict], segments: list[dict]) -> int:
    """Evaluate in protein matches and chain segments overlap.

    :param locations:
    :param segments:
    :return:
    """
    max_overlap = 0

    for loc in locations:
        # We do not consider fragmented matches
        loc_start = loc["fragments"][0]["start"]
        loc_end = max([f["end"] for f in loc["fragments"]])

        for segment in segments:
            seg_start = segment["protein_start"]
            seg_end = segment["protein_end"]
            overlap = min(loc_end, seg_end) - max(loc_start, seg_start) + 1
            max_overlap = max(max_overlap, overlap)

    return max_overlap


match_complete_sql_query = f""" 

    WITH 
    proteins AS (
    SELECT PROTEIN_AC, NAME, DBCODE, CRC64, LEN,
        TO_CHAR(TIMESTAMP, 'YYYY-MM-DD') AS TIMESTAMP,
        FRAGMENT, TO_CHAR(TAX_ID) AS TAX_ID
    FROM INTERPRO.PROTEIN
    WHERE ROWNUM <= 100
    ORDER BY PROTEIN_AC
    ),

    matches AS (
        -- Limit the number of rows from the MATCH and FEATURE_MATCH tables
        SELECT PROTEIN_AC, METHOD_AC, MODEL_AC, POS_FROM, POS_TO, FRAGMENTS, SCORE, DBCODE, EVIDENCE, STATUS
        FROM INTERPRO.MATCH
    )

    SELECT P.PROTEIN_AC, P.NAME, M.DBCODE, P.CRC64, P.LEN, P.TIMESTAMP, P.FRAGMENT, P.TAX_ID,
        M.METHOD_AC, M.MODEL_AC, M.POS_FROM, M.POS_TO, M.FRAGMENTS, M.SCORE, MN.DESCRIPTION, M.STATUS,
        DB.DBSHORT, CE.ABBREV, CT.ABBREV

    FROM proteins P
    LEFT OUTER JOIN matches M
    ON P.PROTEIN_AC = M.PROTEIN_AC

    LEFT OUTER JOIN INTERPRO.METHOD MN
    ON M.METHOD_AC = MN.METHOD_AC

    LEFT OUTER JOIN INTERPRO.CV_DATABASE DB
    ON M.DBCODE = DB.DBCODE

    LEFT OUTER JOIN CV_EVIDENCE CE
    ON M.EVIDENCE = CE.CODE

    LEFT OUTER JOIN CV_ENTRY_TYPE CT
    ON MN.SIG_TYPE = CT.CODE

    ORDER BY P.PROTEIN_AC

    """
