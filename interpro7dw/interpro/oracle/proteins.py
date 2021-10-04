from typing import List, Sequence, Tuple

import cx_Oracle

from interpro7dw.utils import logger, SimpleStore, Store


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


def export_features(url: str, src: str, dst: str, **kwargs):
    tempdir = kwargs.get("tempdir")
    workers = kwargs.get("workers", 1)

    logger.info("starting")
    with Store(src, "r") as store:
        keys = store.file_keys

    with Store(dst, "w", keys=keys, tempdir=tempdir) as store:
        con = cx_Oracle.connect(url)
        cur = con.cursor()
        cur.execute(
            """
            SELECT FM.PROTEIN_AC, FM.METHOD_AC, LOWER(DB.DBSHORT),
                   FM.POS_FROM, FM.POS_TO, FM.SEQ_FEATURE
            FROM INTERPRO.FEATURE_MATCH FM
            INNER JOIN INTERPRO.CV_DATABASE DB ON FM.DBCODE = DB.DBCODE
            """
        )

        for i, rec in enumerate(cur):
            protein_acc = rec[0]
            signature_acc = rec[1]
            database = rec[2]
            pos_start = rec[3]
            pos_end = rec[4]
            seq_feature = rec[5]

            if database == "mobidblt" and seq_feature is None:
                seq_feature = "Consensus Disorder Prediction"

            store.add(protein_acc, {
                signature_acc: {
                    "database": database,
                    "locations": [(pos_start, pos_end, seq_feature)]
                }
            })

            if (i + 1) % 10000000 == 0:
                logger.info(f"{i + 1:>15,}")

        logger.info(f"{i + 1:>15,}")
        cur.close()
        con.close()

        store.merge(workers, apply=_sort_features)
        logger.info(f"temporary files: {store.size / 1024 / 1024:.0f} MB")

    logger.info("done")


def _sort_features(values: Sequence[dict]) -> dict:
    protein = Store.merge_dicts(values)

    for signature in protein.values():
        signature["locations"].sort()

    return protein


def export_proteins(url: str, file: str, **kwargs):
    chunksize = kwargs.get("chunksize", 10000)
    tempdir = kwargs.get("tempdir")
    workers = kwargs.get("workers", 1)

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
            store.merge(workers=workers, apply=store.get_first)

            size = tmpstore.size + store.size

    logger.info(f"temporary files: {size / 1024 / 1024:.0f} MB")
    logger.info("done")


def export_matches(url: str, src: str, dst: str, **kwargs):
    tempdir = kwargs.get("tempdir")
    workers = kwargs.get("workers", 1)

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

                for frag in fragments[5].split(','):
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

        store.merge(workers, apply=_merge_matches)
        logger.info(f"temporary files: {store.size / 1024 / 1024:.0f} MB")

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
        locations.sort(key=lambda l: (l[0]["fragments"]["start"],
                                      l[0]["fragments"]["end"]))
        entries[signature_acc] = locations

    return entries


def export_residues(url: str, src: str, dst: str, **kwargs):
    tempdir = kwargs.get("tempdir")
    workers = kwargs.get("workers", 1)

    logger.info("starting")
    with Store(src, "r") as store:
        keys = store.file_keys

    with Store(dst, "w", keys=keys, tempdir=tempdir) as store:
        con = cx_Oracle.connect(url)
        cur = con.cursor()
        cur.execute(
            """
            SELECT S.PROTEIN_AC, S.METHOD_AC, M.NAME, LOWER(D.DBSHORT),
                   S.DESCRIPTION, S.RESIDUE, S.RESIDUE_START, S.RESIDUE_END
            FROM INTERPRO.SITE_MATCH S
            INNER JOIN INTERPRO.CV_DATABASE D ON S.DBCODE = D.DBCODE
            LEFT OUTER JOIN INTERPRO.METHOD M ON S.METHOD_AC = M.METHOD_AC  
            """
        )

        for i, rec in enumerate(cur):
            protein_acc = rec[0]
            signature_acc = rec[1]
            signature_name = rec[2]
            database = rec[3]
            description = rec[4]
            residue = rec[5]
            pos_start = rec[6]
            pos_end = rec[7]

            store.add(protein_acc, {
                signature_acc: {
                    "name": signature_name,
                    "database": database,
                    "descriptions": {
                        description: [(residue, pos_start, pos_end)]
                    }
                }
            })

            if (i + 1) % 100000000 == 0:
                logger.info(f"{i + 1:>15,}")

        logger.info(f"{i + 1:>15,}")
        cur.close()
        con.close()

        store.merge(workers, apply=_sort_residues)
        logger.info(f"temporary files: {store.size / 1024 / 1024:.0f} MB")

    logger.info("done")


def _sort_residues(values: Sequence[dict]) -> dict:
    protein = Store.merge_dicts(values)

    for signature in protein.values():
        for locations in signature["descriptions"].values():
            locations.sort(key=lambda x: (x[1], x[2]))

    return protein


def export_sequences(url: str, src: str, dst: str, **kwargs):
    tempdir = kwargs.get("tempdir")
    workers = kwargs.get("workers", 1)

    logger.info("starting")
    with Store(src, "r") as store:
        keys = store.file_keys

    with Store(dst, "w", keys=keys, tempdir=tempdir) as store:
        con = cx_Oracle.connect(url)
        cur = con.cursor()
        cur.execute(
            """
            SELECT UX.AC, UP.SEQ_SHORT, UP.SEQ_LONG
            FROM UNIPARC.XREF UX
            INNER JOIN UNIPARC.PROTEIN UP ON UX.UPI = UP.UPI
            WHERE UX.DBID IN (2, 3)  -- Swiss-Prot and TrEMBL
            AND UX.DELETED = 'N'
            """
        )

        for i, (accession, bytes_obj, clob) in enumerate(cur):
            store.add(accession, bytes_obj or clob.read())

            if (i + 1) % 10000000 == 0:
                logger.info(f"{i + 1:>15,}")

        logger.info(f"{i + 1:>15,}")
        cur.close()
        con.close()

        store.merge(workers=workers, apply=store.get_first)
        logger.info(f"temporary files: {store.size / 1024 / 1024:.0f} MB")

    logger.info("done")
