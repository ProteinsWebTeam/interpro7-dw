import os
from typing import Optional

import cx_Oracle

from interpro7dw.utils import logger
from interpro7dw.utils.oracle import lob_as_str
from interpro7dw.utils.store import BasicStore, KVStoreBuilder, KVStore


DC_STATUSES = {
    # Continuous single chain domain
    "S": "CONTINUOUS",
    # N terminus discontinuous
    "N": "N_TERMINAL_DISC",
    # C terminus discontinuous
    "C": "C_TERMINAL_DISC",
    # N and C terminus discontinuous
    "NC": "NC_TERMINAL_DISC"
}


def _load_entries(cur: cx_Oracle.Cursor) -> dict:
    entries = {}
    cur.execute(
        """
        SELECT E.ENTRY_AC, E.NAME, ET.ABBREV, EE.PARENT_AC
        FROM INTERPRO.ENTRY E
        INNER JOIN INTERPRO.CV_ENTRY_TYPE ET
          ON E.ENTRY_TYPE = ET.CODE
        LEFT OUTER JOIN INTERPRO.ENTRY2ENTRY EE
          ON E.ENTRY_AC = EE.ENTRY_AC AND EE.RELATION = 'TY'
        WHERE E.CHECKED = 'Y'
        """
    )

    for rec in cur:
        entries[rec[0]] = {
            "name": rec[1],
            "type": rec[2],
            "parent": rec[3]
        }

    return entries


def _load_signatures(cur: cx_Oracle.Cursor) -> dict:
    signatures = {}
    cur.execute(
        """
        SELECT M.METHOD_AC, M.NAME, M.DESCRIPTION, D.DBSHORT, ET.ABBREV, 
               EVI.ABBREV, EM.ENTRY_AC
        FROM INTERPRO.METHOD M
        INNER JOIN INTERPRO.CV_DATABASE D
          ON M.DBCODE = D.DBCODE
        INNER JOIN INTERPRO.CV_ENTRY_TYPE ET
          ON M.SIG_TYPE = ET.CODE
        INNER JOIN INTERPRO.IPRSCAN2DBCODE I2D 
          ON M.DBCODE = I2D.DBCODE
        INNER JOIN INTERPRO.CV_EVIDENCE EVI
          ON I2D.EVIDENCE = EVI.CODE
        LEFT OUTER JOIN (
            SELECT E.ENTRY_AC, EM.METHOD_AC
            FROM INTERPRO.ENTRY E
            INNER JOIN INTERPRO.ENTRY2METHOD EM
              ON E.ENTRY_AC = EM.ENTRY_AC
            WHERE E.CHECKED = 'Y'
        ) EM ON M.METHOD_AC = EM.METHOD_AC
        """
    )

    for rec in cur:
        signatures[rec[0]] = {
            "short_name": rec[1],
            "name": rec[2],
            "database": rec[3],
            "type": rec[4],
            "evidence": rec[5],
            "entry": rec[6]
        }

    return signatures


def _get_fragments(pos_start: int, pos_end: int, fragments: str) -> list[dict]:
    if fragments:
        result = []
        for frag in fragments.split(','):
            # Format: START-END-STATUS
            s, e, t = frag.split('-')
            result.append({
                "start": int(s),
                "end": int(e),
                "dc-status": DC_STATUSES[t]
            })

        result.sort(key=lambda x: (x["start"], x["end"]))
    else:
        result = [{
            "start": pos_start,
            "end": pos_end,
            "dc-status": DC_STATUSES['S']  # Continuous
        }]

    return result


def condense_locations(locations: list[list[dict]],
                       min_overlap: int = 0.1) -> list[tuple[int, int]]:
    start = end = None
    condensed = []

    """
    Sort locations using their leftmost fragment
    (assume that fragments are sorted in individual locations)
    """
    for fragments in sorted(locations,
                            key=lambda l: (l[0]["start"], l[0]["end"])):
        """
        1) We do not consider fragmented matches
        2) Fragments are sorted by (start, end):
            * `start` of the first frag is guaranteed to be the leftmost one
            * `end` of the last frag is NOT guaranteed to be the rightmost one
                (e.g. [(5, 100), (6, 80)])
        """
        s = fragments[0]["start"]
        e = max([f["end"] for f in fragments])

        if start is None:
            # First location
            start, end = s, e
            continue
        elif e <= end:
            # Current location within "merged" one: nothing to do
            continue
        elif s <= end:
            # Locations are overlapping (at least one residue)
            overlap = min(end, e) - max(start, s) + 1
            shortest = min(end - start, e - s) + 1

            if overlap >= shortest * min_overlap:
                # Merge
                end = e
                continue

        condensed.append((start, end))
        start, end = s, e

    # Adding last location
    condensed.append((start, end))
    return condensed


def export_uniprot_matches(uri: str, proteins_file: str, output: str,
                           processes: int = 1, tempdir: Optional[str] = None):
    logger.info("starting")

    with KVStore(proteins_file) as store:
        keys = store.get_keys()

    with KVStoreBuilder(output, keys=keys, tempdir=tempdir) as store:
        con = cx_Oracle.connect(uri)
        cur = con.cursor()
        cur.execute(
            """
            SELECT PROTEIN_AC, METHOD_AC, MODEL_AC, POS_FROM, 
                   POS_TO, FRAGMENTS, SCORE
            FROM INTERPRO.MATCH
            """
        )
        i = 0
        for prot_acc, sig_acc, mod_acc, start, end, frags, score in cur:
            store.add(prot_acc, (
                sig_acc,
                mod_acc,
                score,
                _get_fragments(start, end, frags)
            ))

            i += 1
            if i % 1e8 == 0:
                logger.info(f"{i:>15,}")

        logger.info(f"{i:>15,}")

        entries = _load_entries(cur)
        signatures = _load_signatures(cur)

        cur.close()
        con.close()

        size = store.get_size()
        store.build(apply=_merge_uniprot_matches,
                    processes=processes,
                    extraargs=[signatures, entries])

        size = max(size, store.get_size())
        logger.info(f"temporary files: {size / 1024 ** 2:.0f} MB")

    logger.info("done")


def _merge_uniprot_matches(matches: list[tuple], signatures: dict,
                           entries: dict) -> tuple[dict, dict]:
    entry_matches = {}
    signature_matches = {}

    for signature_acc, model_acc, score, fragments in matches:
        if signature_acc in signature_matches:
            match = signature_matches[signature_acc]
        else:
            signature = signatures[signature_acc]
            match = signature_matches[signature_acc] = {
                "name": signature["name"],
                "database": signature["database"],
                "type": signature["type"],
                "evidence": signature["evidence"],
                "entry": signature["entry"],
                "locations": []
            }

            if match["entry"] and match["entry"] not in entry_matches:
                entry = entries[match["entry"]]
                entry_matches[match["entry"]] = {
                    "name": entry["name"],
                    "database": "INTERPRO",
                    "type": entry["type"],
                    "parent": entry["parent"],
                    "locations": []
                }

        match["locations"].append({
            "fragments": fragments,
            "model": model_acc or signature_acc,
            "score": score
        })

        if match["entry"]:
            entry_matches[match["entry"]]["locations"].append(fragments)

    """
    Sort signature locations using their leftmost fragment
    (expects individual locations to be sorted by fragment)
    """
    for match in signature_matches.values():
        match["locations"].sort(key=lambda l: (l["fragments"][0]["start"],
                                               l["fragments"][0]["end"]))

    # Merge overlapping matches
    for match in entry_matches.values():
        condensed = []
        for start, end in condense_locations(match["locations"]):
            condensed.append({
                "fragments": [{
                    "start": start,
                    "end": end,
                    "dc-status": DC_STATUSES['S']
                }],
                "model": None,
                "score": None
            })

        match["locations"] = condensed

    return signatures, entries


def export_residues(uri: str, output: str):
    logger.info("starting")

    with BasicStore(output, mode="w") as store:
        for i, item in enumerate(_iter_residues(uri)):
            store.write(item)

            if (i + 1) % 1e7 == 0:
                logger.info(f"{i + 1:>15,}")

        logger.info(f"{i + 1:>15,}")

    logger.info("done")


def _iter_residues(uri: str):
    con = cx_Oracle.connect(uri)
    cur = con.cursor()

    cur.execute("SELECT DBCODE, DBSHORT FROM INTERPRO.CV_DATABASE")
    databases = dict(cur.fetchall())

    cur.execute("SELECT METHOD_AC, NAME FROM INTERPRO.METHOD")
    signatures = dict(cur.fetchall())

    cur.execute(
        """
        SELECT PROTEIN_AC, METHOD_AC, DBCODE, DESCRIPTION, RESIDUE, 
               RESIDUE_START, RESIDUE_END
        FROM INTERPRO.SITE_MATCH
        ORDER BY PROTEIN_AC
        """
    )

    protein_acc = None
    matches = {}
    for prot_acc, sig_acc, dbcode, descr, res, pos_start, pos_end in cur:
        if prot_acc != protein_acc:
            if protein_acc:
                yield protein_acc, _sort_residues(matches)

            protein_acc = prot_acc
            matches = {}

        try:
            signature = matches[sig_acc]
        except KeyError:
            signature = matches[sig_acc] = {
                "name": signatures.get(sig_acc),
                "database": databases[dbcode],
                "descriptions": {}
            }

        try:
            signature["descriptions"][descr].append((res, pos_start, pos_end))
        except KeyError:
            signature["descriptions"][descr] = [(res, pos_start, pos_end)]

    cur.close()
    con.close()

    if protein_acc:
        yield protein_acc, _sort_residues(matches)


def _sort_residues(matches: dict) -> dict:
    for signature in matches.values():
        for locations in signature["descriptions"].values():
            locations.sort(key=lambda x: (x[1], x[2]))

    return matches


def export_features(uri: str, output: str):
    logger.info("starting")

    with BasicStore(output, mode="w") as store:
        for i, item in enumerate(_iter_features(uri)):
            store.write(item)

            if (i + 1) % 1e7 == 0:
                logger.info(f"{i + 1:>15,}")

        store.close()
        logger.info(f"{i + 1:>15,}")

    logger.info("done")


def _iter_features(uri: str):
    con = cx_Oracle.connect(uri)
    cur = con.cursor()
    cur.execute(
        """
        SELECT M.METHOD_AC, M.NAME, D.DBSHORT, EVI.ABBREV
        FROM INTERPRO.FEATURE_METHOD M
        INNER JOIN INTERPRO.CV_DATABASE D
          ON M.DBCODE = D.DBCODE
        INNER JOIN INTERPRO.IPRSCAN2DBCODE I2D
          ON D.DBCODE = I2D.DBCODE
        INNER JOIN INTERPRO.CV_EVIDENCE EVI
          ON I2D.EVIDENCE = EVI.CODE
        """
    )
    features_info = {row[0]: row[1:] for row in cur}

    cur.execute(
        """
        SELECT PROTEIN_AC, METHOD_AC, POS_FROM, POS_TO, SEQ_FEATURE
        FROM INTERPRO.FEATURE_MATCH
        ORDER BY PROTEIN_AC
        """
    )

    protein_acc = None
    features = {}
    for prot_acc, feat_acc, pos_start, pos_end, seq_feature in cur:
        if prot_acc != protein_acc:
            if protein_acc:
                yield protein_acc, _sort_features(features)

            protein_acc = prot_acc
            features = {}

        try:
            feature = features[feat_acc]
        except KeyError:
            name, database, evidence = features_info[feat_acc]
            feature = features[feat_acc] = {
                "name": name,
                "database": database,
                "evidence": evidence,
                "locations": [(pos_start, pos_end, seq_feature)]
            }

        if seq_feature is None and feature["database"].lower() == "mobidblt":
            seq_feature = "Consensus Disorder Prediction"

        feature["locations"].append((pos_start, pos_end, seq_feature))

    cur.close()
    con.close()

    if protein_acc:
        yield protein_acc, _sort_features(features)


def _sort_features(features: dict) -> dict:
    for feature in features.values():
        feature["locations"].sort(key=lambda x: (x[0], x[1]))

    return features


def export_isoforms(uri: str, output: str):
    con = cx_Oracle.connect(uri)
    cur = con.cursor()

    entries = _load_entries(cur)
    signatures = _load_signatures(cur)

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
        SELECT PROTEIN_AC, METHOD_AC, MODEL_AC, SCORE, POS_FROM, POS_TO, 
               FRAGMENTS
        FROM INTERPRO.VARSPLIC_MATCH V
        """
    )

    for var_acc, sig_acc, model_acc, score, pos_start, pos_end, frags in cur:
        try:
            isoform = isoforms[var_acc]
        except KeyError:
            continue

        fragments = _get_fragments(pos_start, pos_end, frags)
        isoform["matches"].append((sig_acc, model_acc, score, fragments))

    cur.close()
    con.close()

    with BasicStore(output, "w") as store:
        for isoform in isoforms.values():
            matches = isoform.pop("matches")
            isoform["matches"] = _merge_uniprot_matches(matches, signatures,
                                                        entries)
            store.write(isoform)


def export_uniparc_matches(uri: str, proteins_file: str, output: str,
                           processes: int = 8, tempdir: Optional[str] = None):
    logger.info("starting")

    with KVStore(proteins_file) as store:
        keys = store.get_keys()

    with KVStoreBuilder(output, keys=keys, tempdir=tempdir) as store:
        con = cx_Oracle.connect(uri)
        cur = con.cursor()

        cur.execute("SELECT METHOD_AC, NAME FROM INTERPRO.METHOD")
        signatures = {acc for acc, in cur}

        # SEQ_FEATURE -> contains the alignment for ProSite and HAMAP
        cur.execute(
            """
            SELECT UPI, METHOD_AC, MODEL_AC, SEQ_START, SEQ_END, SCORE, 
                   SEQ_FEATURE, FRAGMENTS
            FROM IPRSCAN.MV_IPRSCAN
            """
        )

        i = 0
        for rec in cur:
            if rec[1] in signatures:
                store.add(rec[0], rec[1:])

            i += 1
            if i % 1e8 == 0:
                logger.info(f"{i:>15,}")

        logger.info(f"{i:>15,}")

        entries = _load_entries(cur)
        signatures = _load_signatures(cur)

        cur.close()
        con.close()

        size = store.get_size()
        store.build(apply=_merge_uniparc_matches,
                    processes=processes,
                    extraargs=[signatures, entries])

        size = max(size, store.get_size())
        logger.info(f"temporary files: {size / 1024 ** 2:.0f} MB")

    logger.info("done")


def _merge_uniparc_matches(matches: list[tuple], signatures: dict,
                           entries: dict) -> dict:
    signature_matches = {}
    for sig_acc, mod_acc, start, end, score, aln, fragments in matches:
        if sig_acc in signature_matches:
            match = signature_matches[sig_acc]
        else:
            signature = signatures[sig_acc]
            entry_acc = signature["entry"]
            if entry_acc:
                entry = {
                    "accession": entry_acc,
                    "name": entries[entry_acc]["name"],
                    "type": entries[entry_acc]["type"],
                    "parent": entries[entry_acc]["parent"],
                }
            else:
                entry = None

            match = signature_matches[sig_acc] = {
                "name": signature["short_name"],
                "database": signature["database"],
                "evidence": signature["evidence"],
                "entry": entry,
                "model": mod_acc,
                "locations": []
            }

        match["locations"].append((start, end, score, aln, fragments))

    # Sort locations
    for match in signature_matches.values():
        match["locations"].sort(key=lambda x: (x[0], x[1]))

    return signature_matches
