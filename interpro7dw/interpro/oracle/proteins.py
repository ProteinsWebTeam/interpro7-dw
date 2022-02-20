import gzip
import os
from typing import List, Sequence, Tuple

import cx_Oracle

from interpro7dw.utils import logger
from interpro7dw.utils.oracle import lob_as_str
from interpro7dw.utils.store import copy_dict, loadobj, SimpleStore, Store, StoreDirectLoader


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
                yield protein_acc, _sort_features2(features)

            protein_acc = prot_acc
            features = {}

        name, database, evidence = features_info[feat_acc]

        if database.lower() == "mobidblt" and seq_feature is None:
            seq_feature = "Consensus Disorder Prediction"

        try:
            features[feat_acc]["locations"].append((pos_start, pos_end,
                                                    seq_feature))
        except KeyError:
            features[feat_acc] = {
                "name": name,
                "database": database,
                "evidence": evidence,
                "locations": [(pos_start, pos_end, seq_feature)]
            }

    cur.close()
    con.close()

    if protein_acc:
        yield protein_acc, _sort_features2(features)


def _sort_features2(features: dict) -> dict:
    for feature in features.values():
        feature["locations"].sort()

    return features


def export_features2(uri: str, output: str):
    logger.info("starting")
    with StoreDirectLoader(output) as store:
        for i, (protein_acc, features) in enumerate(_iter_features(uri)):
            store.add(protein_acc, features)

            if (i + 1) % 1e7 == 0:
                logger.info(f"{i + 1:>15,}")

        logger.info(f"{i + 1:>15,}")

        store.close()

    logger.info("done")


def export_features(url: str, src: str, dst: str, **kwargs):
    tempdir = kwargs.get("tempdir")

    logger.info("starting")
    with Store(src, "r") as store:
        keys = store.file_keys

    with Store(dst, "w", keys=keys, tempdir=tempdir) as store:
        con = cx_Oracle.connect(url)
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
        features = {}
        for row in cur:
            features[row[0]] = row[1:]

        cur.execute(
            """
            SELECT PROTEIN_AC, METHOD_AC, POS_FROM, POS_TO, SEQ_FEATURE
            FROM INTERPRO.FEATURE_MATCH
            """
        )

        for i, rec in enumerate(cur):
            protein_acc = rec[0]
            feature_acc = rec[1]
            pos_start = rec[2]
            pos_end = rec[3]
            seq_feature = rec[4]

            name, database, evidence = features[feature_acc]

            if database.lower() == "mobidblt" and seq_feature is None:
                seq_feature = "Consensus Disorder Prediction"

            store.add(protein_acc, {
                feature_acc: {
                    "name": name,
                    "database": database,
                    "evidence": evidence,
                    "locations": [(pos_start, pos_end, seq_feature)]
                }
            })

            if (i + 1) % 100000000 == 0:
                logger.info(f"{i + 1:>15,}")

        logger.info(f"{i + 1:>15,}")
        cur.close()
        con.close()

        store.merge(apply=_sort_features)
        logger.info(f"temporary files: {store.size / 1024 / 1024:.0f} MB")

    logger.info("done")


def _sort_features(values: Sequence[dict]) -> dict:
    protein = {}
    for value in values:
        copy_dict(value, protein)

    for signature in protein.values():
        signature["locations"].sort()

    return protein


def export_isoforms(url: str, dst: str):
    con = cx_Oracle.connect(url)
    cur = con.cursor()
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


def export_proteins(url: str, file: str, **kwargs):
    chunksize = kwargs.get("chunksize", 10000)
    tempdir = kwargs.get("tempdir")

    logger.info("creating temporary store")
    keys = []
    with SimpleStore(tempdir=tempdir) as tmpstore:
        con = cx_Oracle.connect(url)
        cur = con.cursor()
        cur.execute(
            """
            SELECT 
              PROTEIN_AC, NAME, DBCODE, CRC64, LEN, 
              TO_CHAR(TIMESTAMP, 'YYYY-MM-DD'), FRAGMENT, TO_CHAR(TAX_ID) 
            FROM INTERPRO.PROTEIN
            ORDER BY PROTEIN_AC
            """
        )

        for i, rec in enumerate(cur):
            if i % chunksize == 0:
                accession = rec[0]
                keys.append(accession)

            tmpstore.add(rec)

        cur.close()
        con.close()

        logger.info("creating final store")

        with Store(file, "w", keys=keys, tempdir=tempdir) as store:
            for i, rec in enumerate(tmpstore):
                store.add(rec[0], {
                    "identifier": rec[1],
                    "reviewed": rec[2] == 'S',
                    "crc64": rec[3],
                    "length": rec[4],
                    "date": rec[5],
                    "fragment": rec[6] == 'Y',
                    "taxid": rec[7]
                })

                if (i + 1) % 10000000 == 0:
                    logger.info(f"{i + 1:>15,}")

            logger.info(f"{i + 1:>15,}")
            store.merge(apply=store.get_first)

            size = tmpstore.size + store.size

    logger.info(f"temporary files: {size / 1024 / 1024:.0f} MB")
    logger.info("done")


def export_matches(url: str, src: str, dst: str, **kwargs):
    tempdir = kwargs.get("tempdir")

    logger.info("starting")
    with Store(src, "r") as store:
        keys = store.file_keys

    with Store(dst, "w", keys=keys, tempdir=tempdir) as store:
        con = cx_Oracle.connect(url)
        cur = con.cursor()
        cur.execute(
            """
            SELECT M.PROTEIN_AC, M.METHOD_AC, M.MODEL_AC, M.POS_FROM, 
                   M.POS_TO, M.FRAGMENTS, M.SCORE, E.ENTRY_AC
            FROM INTERPRO.MATCH M
            LEFT OUTER JOIN (
              SELECT E.ENTRY_AC, EM.METHOD_AC
              FROM INTERPRO.ENTRY E
              INNER JOIN INTERPRO.ENTRY2METHOD EM
                ON E.ENTRY_AC = EM.ENTRY_AC
              WHERE E.CHECKED = 'Y'
            ) E ON M.METHOD_AC = E.METHOD_AC
            """
        )

        for i, rec in enumerate(cur):
            protein_acc = rec[0]
            signature_acc = rec[1]
            model_acc = rec[2]
            pos_start = rec[3]
            pos_end = rec[4]
            fragments = rec[5]
            score = rec[6]
            entry_acc = rec[7]

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

            store.add(protein_acc, (signature_acc, model_acc, score,
                                    fragments, entry_acc))

            if (i + 1) % 100000000 == 0:
                logger.info(f"{i + 1:>15,}")

        logger.info(f"{i + 1:>15,}")
        cur.close()
        con.close()

        store.merge(apply=_merge_matches)
        logger.info(f"temporary files: {store.size / 1024 / 1024:.0f} MB")

    logger.info("done")


def _iter_matches(uri: str):
    con = cx_Oracle.connect(uri)
    cur = con.cursor()
    cur.execute(
        """
        SELECT EM.METHOD_AC, E.ENTRY_AC
        FROM INTERPRO.ENTRY E
        INNER JOIN INTERPRO.ENTRY2METHOD EM
            ON E.ENTRY_AC = EM.ENTRY_AC
        WHERE E.CHECKED = 'Y'        
        """
    )
    integrated = dict(cur.fetchall())

    cur.execute(
        """
        SELECT PROTEIN_AC, METHOD_AC, MODEL_AC, POS_FROM, POS_TO, FRAGMENTS, 
               SCORE
        FROM INTERPRO.MATCH
        ORDER BY PROTEIN_AC
        """
    )

    protein_acc = None
    matches = []
    for prot_acc, sig_acc, model_acc, pos_start, pos_end, frags, score in cur:
        if prot_acc != protein_acc:
            if protein_acc:
                yield protein_acc, _merge_matches(matches)

            protein_acc = prot_acc
            matches = []

        if frags:
            fragments = []
            for f in frags.split(','):
                # Format: START-END-STATUS
                s, e, t = f.split('-')
                fragments.append({
                    "start": int(s),
                    "end": int(e),
                    "dc-status": DC_STATUSES[t]
                })
        else:
            fragments = [{
                "start": pos_start,
                "end": pos_end,
                "dc-status": DC_STATUSES['S']  # Continuous
            }]

        matches.append((sig_acc, model_acc, score, fragments,
                        integrated.get(sig_acc)))

    cur.close()
    con.close()

    if protein_acc:
        yield protein_acc, _merge_matches(matches)


def export_matches2(uri: str, output: str):
    logger.info("starting")
    with StoreDirectLoader(output) as store:
        for i, (protein_acc, matches) in enumerate(_iter_matches(uri)):
            store.add(protein_acc, matches)

            if (i + 1) % 1e7 == 0:
                logger.info(f"{i + 1:>15,}")

        logger.info(f"{i + 1:>15,}")

        store.close()

    logger.info("done")


def condense_locations(locations: Sequence[Sequence[dict]],
                       min_overlap: int = 0.1) -> List[Tuple[int, int]]:
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


def _merge_matches(values: Sequence[tuple]) -> dict:
    entries = {}
    signatures = {}
    for signature_acc, model_acc, score, fragments, entry_acc in values:
        fragments.sort(key=lambda x: (x["start"], x["end"]))

        try:
            s = signatures[signature_acc]
        except KeyError:
            s = signatures[signature_acc] = []

        s.append({
            "fragments": fragments,
            "model": model_acc or signature_acc,
            "score": score
        })

        if entry_acc:
            try:
                entries[entry_acc].append(fragments)
            except KeyError:
                entries[entry_acc] = [fragments]

    # Merge overlapping matches
    for entry_acc, locations in entries.items():
        condensed = []
        for start, end in condense_locations(locations):
            condensed.append({
                "fragments": [{
                    "start": start,
                    "end": end,
                    "dc-status": DC_STATUSES['S']
                }],
                "model": None,
                "score": None
            })

        entries[entry_acc] = condensed

    # Add signatures
    for signature_acc, locations in signatures.items():
        # Sort locations using their leftmost fragment (fragments are sorted)
        locations.sort(key=lambda l: (l["fragments"][0]["start"],
                                      l["fragments"][0]["end"]))
        entries[signature_acc] = locations

    return entries


def _iter_residues(uri: str):
    con = cx_Oracle.connect(uri)
    cur = con.cursor()

    cur.execute("SELECT DBCODE, LOWER(DBSHORT) FROM INTERPRO.CV_DATABASE")
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


def export_residues(uri: str, output: str):
    logger.info("starting")

    with SimpleStore(file=output) as store:
        for i, item in enumerate(_iter_residues(uri)):
            store.add(item)

            if (i + 1) % 1e7 == 0:
                logger.info(f"{i + 1:>15,}")

        logger.info(f"{i + 1:>15,}")

    logger.info("done")


def export_sequences(url: str, src: str, dst: str, **kwargs):
    tempdir = kwargs.get("tempdir")

    logger.info("starting")
    with Store(src, "r") as store:
        keys = store.file_keys

    with Store(dst, "w", keys=keys, tempdir=tempdir) as store:
        con = cx_Oracle.connect(url)
        cur = con.cursor()
        cur.outputtypehandler = lob_as_str
        cur.execute(
            """
            SELECT UX.AC, UP.SEQ_SHORT, UP.SEQ_LONG
            FROM UNIPARC.XREF UX
            INNER JOIN UNIPARC.PROTEIN UP ON UX.UPI = UP.UPI
            WHERE UX.DBID IN (2, 3)  -- Swiss-Prot and TrEMBL
            AND UX.DELETED = 'N'
            """
        )

        for i, (accession, seq_short, seq_long) in enumerate(cur):
            sequence = seq_short or seq_long
            store.add(accession, gzip.compress(sequence.encode("utf-8")))

            if (i + 1) % 10000000 == 0:
                logger.info(f"{i + 1:>15,}")

        logger.info(f"{i + 1:>15,}")
        cur.close()
        con.close()

        store.merge(apply=store.get_first)
        logger.info(f"temporary files: {store.size / 1024 / 1024:.0f} MB")

    logger.info("done")


def export_uniparc(uri: str, entries_file: str, proteins_dst: str, **kwargs):
    chunksize = kwargs.get("chunksize", 10000)
    tempdir = kwargs.get("tempdir")

    con = cx_Oracle.connect(uri)
    cur = con.cursor()

    logger.info("exporting UniParc proteins")
    proteins_tmp = f"{proteins_dst}.proteins.tmp"
    keys = []
    with SimpleStore(proteins_tmp) as store:
        cur.execute(
            """
            SELECT UPI, LEN, CRC64
            FROM UNIPARC.PROTEIN
            ORDER BY UPI
            """
        )

        for i, (upi, length, crc64) in enumerate(cur):
            if i % chunksize == 0:
                keys.append(upi)

            store.add((upi, length, crc64))

            if (i + 1) % 1e8 == 0:
                logger.info(f"{i + 1:>15,}")

        logger.info(f"{i + 1:>15,}")

    logger.info("exporting UniParc matches")
    matches_tmp = f"{proteins_dst}.matches.tmp"
    with Store(matches_tmp, "w", keys=keys, tempdir=tempdir) as store:
        cur.execute(
            """
            SELECT MA.UPI, MA.METHOD_AC, MA.MODEL_AC,
                   MA.SEQ_START, MA.SEQ_END, MA.SCORE, MA.SEQ_FEATURE,
                   MA.FRAGMENTS
            FROM IPRSCAN.MV_IPRSCAN MA
            INNER JOIN INTERPRO.METHOD ME
              ON MA.METHOD_AC = ME.METHOD_AC
            """
        )

        for i, rec in enumerate(cur):
            upi = rec[0]
            store.add(upi, rec[1:])

            if (i + 1) % 1e9 == 0:
                logger.info(f"{i + 1:>15,}")

        logger.info(f"{i + 1:>15,}")

        store.merge(apply=_merge_uniparc_matches)

        logger.info(f"temporary files: {store.size / 1024 / 1024:.0f} MB")

    cur.close()
    con.close()

    logger.info("loading entries")
    entries = loadobj(entries_file)

    logger.info("writing final file")
    with SimpleStore(file=proteins_dst) as store:
        with SimpleStore(proteins_tmp) as st1, Store(matches_tmp, "r") as st2:
            for i, (upi, length, crc64) in enumerate(st1):
                matches = []

                for signature_acc, model_acc, locations in st2.get(upi, []):
                    try:
                        entry = entries[signature_acc]
                    except KeyError:
                        continue

                    if entry.integrated_in:
                        interpro_entry = (
                            entry.integrated_in,
                            entries[entry.integrated_in].name,
                            entries[entry.integrated_in].type,
                            entries[entry.integrated_in].relations[0]
                        )
                    else:
                        interpro_entry = None

                    matches.append((
                        signature_acc,
                        entry.name or signature_acc,
                        entry.source_database,
                        entry.evidence,
                        model_acc,
                        interpro_entry,
                        locations
                    ))

                store.add((upi, length, crc64, matches))

                if (i + 1) % 1e8 == 0:
                    logger.info(f"{i + 1:>15,}")

        logger.info(f"{i + 1:>15,}")

    os.unlink(proteins_tmp)
    os.unlink(matches_tmp)
    logger.info("done")


def _merge_uniparc_matches(matches: Sequence[tuple]) -> list:
    signatures = {}
    for acc, model, start, end, score, seq_feature, fragments in matches:
        try:
            s = signatures[acc]
        except KeyError:
            s = signatures[acc] = {
                "model": model or acc,
                "locations": []
            }
        finally:
            s["locations"].append((start, end, score, seq_feature, fragments))

    matches = []
    for acc in sorted(signatures):
        s = signatures[acc]
        s["locations"].sort()
        matches.append((acc, s["model"], s["locations"]))

    return matches
