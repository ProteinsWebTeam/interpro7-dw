import gzip

import oracledb

from interpro7dw.utils import logger
from interpro7dw.utils.oracle import lob_as_str
from interpro7dw.utils.store import KVStoreBuilder, KVStore


def export_uniprot_proteins(uri: str, output: str):
    logger.info("starting")
    with KVStoreBuilder(output, keys=[], cachesize=10000) as store:
        con = oracledb.connect(uri)
        cur = con.cursor()
        cur.execute(
            """
            SELECT PROTEIN_AC, NAME, DBCODE, CRC64, LEN, 
                   TO_CHAR(TIMESTAMP, 'YYYY-MM-DD'), FRAGMENT, TO_CHAR(TAX_ID) 
            FROM INTERPRO.PROTEIN
            ORDER BY PROTEIN_AC
            """
        )

        for i, rec in enumerate(cur):
            store.append(rec[0], {
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

        cur.close()
        con.close()
        store.close()
        logger.info(f"{i + 1:>15,}")

    logger.info("done")


def export_uniprot_sequences(uri: str, kvstore: str, output: str,
                             tempdir: str | None = None):
    logger.info("starting")
    with KVStore(kvstore) as s:
        keys = s.get_keys()

    with KVStoreBuilder(output, keys=keys, tempdir=tempdir) as store:
        con = oracledb.connect(uri)
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

            if (i + 1) % 1e7 == 0:
                logger.info(f"{i + 1:>15,}")

        logger.info(f"{i + 1:>15,}")
        cur.close()
        con.close()

        store.build(apply=store.get_first)
        logger.info(f"temporary files: {store.get_size() / 1024 ** 2:.0f} MB")

    logger.info("done")


def export_uniparc_proteins(uri: str, output: str):
    logger.info("starting")
    with KVStoreBuilder(output, keys=[], cachesize=10000) as store:
        con = oracledb.connect(uri)
        cur = con.cursor()
        cur.execute(
            """
            SELECT UPI, LEN, CRC64
            FROM UNIPARC.PROTEIN
            ORDER BY UPI
            """
        )

        for i, (upi, length, crc64) in enumerate(cur):
            store.append(upi, (length, crc64))

            if (i + 1) % 1e8 == 0:
                logger.info(f"{i + 1:>15,}")

        cur.close()
        con.close()

        store.close()
        logger.info(f"{i + 1:>15,}")

    logger.info("done")
