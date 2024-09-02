import gzip
import os
import pickle
import multiprocessing as mp

import oracledb

from interpro7dw.utils import logger
from interpro7dw.utils.store import KVStore


def drop_table(table_name: str, cur: oracledb.Cursor):
    try:
        cur.execute(f"DROP TABLE {table_name}")
    except oracledb.DatabaseError as exception:
        error_obj, = exception.args

        # ORA-00942: table or view does not exist
        # ORA-08103: object no longer exists
        if error_obj.code not in (942, 8103):
            raise exception


def create_md5_table(uri: str, proteins_file: str):
    con = oracledb.connect(uri)
    cur = con.cursor()

    drop_table('IPRSCAN.LOOKUP_MD5', cur)
    cur.execute(
        """
        CREATE TABLE IPRSCAN.LOOKUP_MD5 NOLOGGING
        AS
        SELECT MD5 FROM UNIPARC.PROTEIN
        WHERE 1 = 0
        """
    )

    with KVStore(proteins_file) as proteins:
        rows = []
        for _, _, md5 in proteins.values():
            rows.append((md5,))

            if len(rows) == 10000:
                cur.executemany(
                    """
                    INSERT /*+ APPEND */ INTO IPRSCAN.LOOKUP_MD5
                    VALUES (:1)
                    """,
                    rows
                )
                rows.clear()
                con.commit()

        if rows:
            cur.executemany(
                """
                INSERT /*+ APPEND */ INTO IPRSCAN.LOOKUP_MD5
                VALUES (:1)
                """,
                rows
            )
            rows.clear()
            con.commit()

    cur.execute(
        """
        CREATE UNIQUE INDEX PK_LOOKUP_MD5
        ON IPRSCAN.LOOKUP_MD5 (MD5)
        TABLESPACE IPRSCAN_IND
        NOLOGGING
        """
    )
    cur.close()
    con.close()


def create_matches_table(uri: str, proteins_file: str, workdir: str,
                         processes: int = 8):
    export_workers = []
    queue1 = mp.Queue()
    queue2 = mp.Queue()

    for _ in range(processes):
        processdir = os.path.join(workdir, f"matches-{i}")
        os.makedirs(processdir, exist_ok=True)
        p = mp.Process(target=export_matches,
                       args=(uri, proteins_file, processdir, queue1, queue2))
        p.start()
        export_workers.append(p)

    with KVStore(proteins_file) as store:
        keys = store.get_keys()

        for i, start in enumerate(keys):
            try:
                stop = keys[i + 1]
            except IndexError:
                stop = None
            finally:
                queue1.put((start, stop))

        for _ in export_workers:
            queue1.put(None)

    con = oracledb.connect(uri)
    cur = con.cursor()
    drop_table("IPRSCAN.LOOKUP_MATCH", cur)
    cur.execute(
        """
        CREATE TABLE IPRSCAN.LOOKUP_MATCH (
            MD5 VARCHAR2(32) NOT NULL,
            SIGNATURE_LIBRARY_NAME VARCHAR2(25),
            SIGNATURE_LIBRARY_RELEASE VARCHAR2(20) NOT NULL,
            SIGNATURE_ACCESSION VARCHAR2(255),
            MODEL_ACCESSION VARCHAR2(255),
            SEQ_START NUMBER(10,0),
            SEQ_END NUMBER(10,0),
            FRAGMENTS VARCHAR2(400),
            SEQUENCE_SCORE BINARY_DOUBLE,
            SEQUENCE_EVALUE BINARY_DOUBLE,
            HMM_BOUNDS VARCHAR2(25),
            HMM_START NUMBER,
            HMM_END NUMBER,
            HMM_LENGTH NUMBER,
            ENVELOPE_START NUMBER,
            ENVELOPE_END NUMBER,
            SCORE BINARY_DOUBLE,
            EVALUE BINARY_DOUBLE,
            SEQ_FEATURE VARCHAR2(4000)            
        ) NOLOGGING
        """
    )

    running = len(export_workers)
    insert_workers = []
    while running:
        obj = queue2.get()
        if obj is not None:
            # A file is ready
            if insert_workers:
                queue1.put(obj)
            else:
                _insert_matches(cur, obj)
        else:
            if len(insert_workers) == 0:
                cur.close()
                con.close()

            running -= 1
            p = mp.Process(target=insert_matches,
                           args=(uri, queue1))
            p.start()
            insert_workers.append(p)

    for _ in insert_workers:
        queue1.put(None)

    for p in insert_workers:
        p.join()

    con = oracledb.connect(uri)
    cur = con.cursor()
    cur.execute(
        """
        CREATE INDEX I_LOOKUP_MATCH
        ON IPRSCAN.LOOKUP_MATCH (MD5)
        TABLESPACE IPRSCAN_IND
        NOLOGGING
        """
    )
    cur.close()
    con.close()


def export_matches(uri: str, proteins_file: str, outdir: str,
                   inqueue: mp.Queue, outqueue: mp.Queue):
    con = oracledb.connect(uri)
    cur = con.cursor()

    appls = get_i5_appls(cur)

    with KVStore(proteins_file) as proteins:
        for start, stop in iter(inqueue.get, None):
            if stop is not None:
                where = "UPI >= :1 AND UPI < :2"
                params = [start, stop]
            else:
                where = "UPI >= :1"
                params = [start]

            cur.execute(
                f"""
                SELECT UPI, ANALYSIS_ID, METHOD_AC, MODEL_AC, SEQ_START, 
                       SEQ_END, FRAGMENTS, SEQSCORE, SEQEVALUE, 
                       HMM_BOUNDS, HMM_START, HMM_END, HMM_LENGTH, 
                       ENVELOPE_START, ENVELOPE_END, SCORE, EVALUE, SEQ_FEATURE
                FROM IPRSCAN.MV_IPRSCAN
                WHERE {where}
                """,
                params
            )

            matches = []
            for row in cur.fetchall():
                try:
                    _, _, md5 = proteins[row[0]]
                except KeyError:
                    continue
                else:
                    dbname, dbversion = appls[row[1]]
                    matches.append((
                        md5,
                        dbname,
                        dbversion,
                        *row[2:]
                    ))

            file = os.path.join(outdir, f"match-{start}")
            with gzip.open(file, "wb") as fh:
                pickle.dump(matches, fh, pickle.HIGHEST_PROTOCOL)

            outqueue.put(file)
            i += 1

    cur.close()
    con.close()
    outqueue.put(None)


def get_i5_appls(cur: oracledb.Cursor) -> dict[int, tuple[str, str]]:
    cur.execute(
        """
        SELECT I2D.IPRSCAN_SIG_LIB_REL_ID,
               DECODE(D.DBNAME, 
                      'CATH-Gene3D', 'GENE3D', 
                      UPPER(REPLACE(D.DBNAME, ' ', '_'))),
               V.VERSION
        FROM INTERPRO.IPRSCAN2DBCODE I2D
        INNER JOIN INTERPRO.CV_DATABASE D ON I2D.DBCODE = D.DBCODE
        INNER JOIN INTERPRO.DB_VERSION V ON D.DBCODE = V.DBCODE
        """
    )
    return {row[0]: row[1:] for row in cur.fetchall()}


def insert_matches(uri: str, queue: mp.Queue):
    con = oracledb.connect(uri)
    cur = con.cursor()

    for file in iter(queue.get, None):
        _insert_matches(cur, file)

    cur.close()
    con.close()


def _insert_matches(cur: oracledb.Cursor, file: str, batchsize: int = 10000):
    with gzip.open(file, "rb") as fh:
        records = pickle.load(fh)

    for i in range(0, len(records), batchsize):
        cur.executemany(
            """
            INSERT /*+ APPEND */ INTO IPRSCAN.LOOKUP_MATCH
            VALUES (:1, :2, :3, :4, :5, :6, :7, :8, :9, :10, :11, 
                    :12, :13, :14, :15, :16, :17, :18, :19)
            """,
            records[i:i+batchsize]
        )
        cur.connection.commit()

    # os.unlink(file)


def create_site_table(uri: str, proteins_file: str, workdir: str,
                      processes: int = 8):
    os.makedirs(workdir, exist_ok=True)

    export_workers = []
    queue1 = mp.Queue()
    queue2 = mp.Queue()

    for _ in range(processes):
        processdir = os.path.join(workdir, f"sites-{i}")
        os.makedirs(processdir, exist_ok=True)
        p = mp.Process(target=export_sites,
                       args=(uri, proteins_file, processdir, queue1, queue2))
        p.start()
        export_workers.append(p)

    with KVStore(proteins_file) as store:
        keys = store.get_keys()

        for i, start in enumerate(keys):
            try:
                stop = keys[i + 1]
            except IndexError:
                stop = None
            finally:
                queue1.put((start, stop))

        for _ in export_workers:
            queue1.put(None)

    con = oracledb.connect(uri)
    cur = con.cursor()
    drop_table("IPRSCAN.LOOKUP_SITE", cur)
    cur.execute(
        """
        CREATE TABLE IPRSCAN.LOOKUP_SITE (
            MD5 VARCHAR2(32) NOT NULL,
            SIGNATURE_LIBRARY_NAME VARCHAR2(25),
            SIGNATURE_LIBRARY_RELEASE VARCHAR2(20) NOT NULL,
            SIGNATURE_ACCESSION VARCHAR2(255),
            LOC_START NUMBER(10,0) NOT NULL,
            LOC_END NUMBER(10,0) NOT NULL,
            NUM_SITES NUMBER,
            RESIDUE VARCHAR2(1000),
            RESIDUE_START NUMBER,
            RESIDUE_END NUMBER,
            DESCRIPTION VARCHAR2(255)
        ) NOLOGGING
        """
    )

    running = len(export_workers)
    insert_workers = []
    while running:
        obj = queue2.get()
        if obj is not None:
            # A file is ready
            if insert_workers:
                queue1.put(obj)
            else:
                _insert_sites(cur, obj)
        else:
            if len(insert_workers) == 0:
                cur.close()
                con.close()

            running -= 1
            p = mp.Process(target=insert_sites,
                           args=(uri, queue1))
            p.start()
            insert_workers.append(p)

    for _ in insert_workers:
        queue1.put(None)

    for p in insert_workers:
        p.join()

    con = oracledb.connect(uri)
    cur = con.cursor()
    cur.execute(
        """
        CREATE INDEX I_LOOKUP_SITE
        ON IPRSCAN.LOOKUP_SITE (MD5)
        TABLESPACE IPRSCAN_IND
        NOLOGGING
        """
    )
    cur.close()
    con.close()


def export_sites(uri: str, proteins_file: str, outdir: str,
                 inqueue: mp.Queue, outqueue: mp.Queue):
    con = oracledb.connect(uri)
    cur = con.cursor()

    appls = get_i5_appls(cur)

    with KVStore(proteins_file) as proteins:
        for start, stop in iter(inqueue.get, None):
            if stop is not None:
                where = "UPI >= :1 AND UPI < :2"
                params = [start, stop]
            else:
                where = "UPI >= :1"
                params = [start]

            cur.execute(
                f"""
                SELECT UPI, ANALYSIS_ID, METHOD_AC, LOC_START, LOC_END, 
                       NUM_SITES, RESIDUE, RESIDUE_START, RESIDUE_END, 
                       DESCRIPTION
                FROM IPRSCAN.SITE
                WHERE {where}
                """,
                params
            )

            sites = []
            for row in cur.fetchall():
                try:
                    _, _, md5 = proteins[row[0]]
                except KeyError:
                    continue
                else:
                    dbname, dbversion = appls[row[1]]

                    sites.append((
                        md5,
                        dbname,
                        dbversion,
                        *row[2:]
                    ))

            file = os.path.join(outdir, f"site-{start}")
            with gzip.open(file, "wb") as fh:
                pickle.dump(sites, fh, pickle.HIGHEST_PROTOCOL)

            outqueue.put(file)

    cur.close()
    con.close()
    outqueue.put(None)


def insert_sites(uri: str, queue: mp.Queue):
    con = oracledb.connect(uri)
    cur = con.cursor()

    for file in iter(queue.get, None):
        _insert_sites(cur, file)

    cur.close()
    con.close()


def _insert_sites(cur: oracledb.Cursor, file: str, batchsize: int = 10000):
    with gzip.open(file, "rb") as fh:
        records = pickle.load(fh)

    for i in range(0, len(records), batchsize):
        cur.executemany(
            """
            INSERT /*+ APPEND */ INTO IPRSCAN.LOOKUP_SITE
            VALUES (:1, :2, :3, :4, :5, :6, :7, :8, :9, :10, :11)
            """,
            records[i:i + batchsize]
        )
        cur.connection.commit()

    # os.unlink(file)
