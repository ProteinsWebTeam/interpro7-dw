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
