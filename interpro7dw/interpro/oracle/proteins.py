import cx_Oracle

from interpro7dw.utils import logger, SimpleStore, Store


def export(url: str, file: str, **kwargs):
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
              PROTEIN_AC, NAME, DBCODE, LEN, FRAGMENT, 
              TO_CHAR(TAX_ID), CRC64
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
                    "length": rec[3],
                    "fragment": rec[4] == 'Y',
                    "taxid": rec[5],
                    "crc64": rec[6]
                })

                if (i + 1) % 10000000 == 0:
                    logger.info(f"{i + 1:>15,}")

            logger.info(f"{i + 1:>15,}")
            store.merge(workers=workers, apply=store.get_first)

        size = tmpstore.size + store.size

    logger.info(f"temporary files: {size / 1024 / 1024:.0f} MB")
    logger.info("done")
