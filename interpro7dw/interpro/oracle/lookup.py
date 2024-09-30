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


def create_md5_table(uri: str, proteins_file: str, batchsize: int = 10000):
    logger.info("starting")

    con = oracledb.connect(uri)
    cur = con.cursor()

    drop_table('IPRSCAN.LOOKUP_MD5', cur)
    cur.execute(
        """
        CREATE TABLE IPRSCAN.LOOKUP_MD5 (
            MD5 VARCHAR2(32) NOT NULL
        ) COMPRESS NOLOGGING
        """
    )

    with KVStore(proteins_file) as proteins:
        rows = []
        for _, _, md5 in proteins.values():
            rows.append((md5,))

            if len(rows) == batchsize:
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

    logger.info("creating index")
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
    logger.info("done")


def create_matches_table(uri: str, proteins_file: str, processes: int = 8,
                         batchsize: int = 10000):
    logger.info("starting")
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
        ) 
        COMPRESS NOLOGGING
        """
    )
    cur.close()
    con.close()

    workers = []
    queue1 = mp.Queue()
    queue2 = mp.Queue()
    for i in range(processes):
        p = mp.Process(target=insert_matches,
                       args=(uri, proteins_file, queue1, batchsize, queue2))
        p.start()
        workers.append(p)

    with KVStore(proteins_file) as store:
        keys = store.get_keys()
        num_tasks = len(keys)

        for i, start in enumerate(keys):
            try:
                stop = keys[i + 1]
                incl_stop = False
            except IndexError:
                stop = store.max()
                incl_stop = True
            finally:
                queue1.put((start, stop, incl_stop))

        for _ in workers:
            queue1.put(None)

    monitor(len(workers), num_tasks, queue2)
    logger.info("done")


def insert_matches(uri: str, proteins_file: str, inqueue: mp.Queue,
                   batchsize: int, outqueue: mp.Queue):
    con = oracledb.connect(uri)
    cur1 = con.cursor()  # select
    cur2 = con.cursor()  # insert
    appls = get_i5_appls(cur1)

    with KVStore(proteins_file) as proteins:
        for start, stop, incl_stop in iter(inqueue.get, None):
            if incl_stop:
                where = "UPI BETWEEN :1 AND :2"
            else:
                where = "UPI >= :1 AND UPI < :2"

            cur1.execute(
                f"""
                SELECT UPI, ANALYSIS_ID, METHOD_AC, MODEL_AC, SEQ_START, 
                       SEQ_END, FRAGMENTS, SEQSCORE, SEQEVALUE, 
                       HMM_BOUNDS, HMM_START, HMM_END, HMM_LENGTH, 
                       ENVELOPE_START, ENVELOPE_END, SCORE, EVALUE, SEQ_FEATURE
                FROM IPRSCAN.MV_IPRSCAN
                WHERE {where}
                """,
                [start, stop]
            )

            while rows := cur1.fetchmany(batchsize):
                matches = []
                for row in rows:
                    _, _, md5 = proteins[row[0]]
                    dbname, dbversion = appls[row[1]]
                    matches.append((
                        md5,
                        dbname,
                        dbversion,
                        *row[2:]
                    ))

                cur2.executemany(
                    """
                    INSERT /*+ APPEND */ INTO IPRSCAN.LOOKUP_MATCH
                    VALUES (:1, :2, :3, :4, :5, :6, :7, :8, :9, :10, :11, 
                            :12, :13, :14, :15, :16, :17, :18, :19)
                    """,
                    matches
                )
                con.commit()

            outqueue.put(False)

    cur1.close()
    cur2.close()
    con.close()
    outqueue.put(True)


def index_matches(uri: str):
    logger.info("indexing")
    con = oracledb.connect(uri)
    cur = con.cursor()
    cur.execute(
        """
        CREATE INDEX I_LOOKUP_MATCH
        ON IPRSCAN.LOOKUP_MATCH (MD5)
        NOLOGGING
        TABLESPACE IPRSCAN_IND
        """
    )
    cur.close()
    con.close()
    logger.info("done")


def create_sites_table(uri: str, proteins_file: str, processes: int = 8,
                       batchsize: int = 10000):
    logger.info("starting")
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
        )
        COMPRESS NOLOGGING
        """
    )
    cur.close()
    con.close()

    workers = []
    queue1 = mp.Queue()
    queue2 = mp.Queue()
    for i in range(processes):
        p = mp.Process(target=insert_sites,
                       args=(uri, proteins_file, queue1, batchsize, queue2))
        p.start()
        workers.append(p)

    with KVStore(proteins_file) as store:
        keys = store.get_keys()
        num_tasks = len(keys)

        for i, start in enumerate(keys):
            try:
                stop = keys[i + 1]
                incl_stop = False
            except IndexError:
                stop = store.max()
                incl_stop = True
            finally:
                queue1.put((start, stop, incl_stop))

        for _ in workers:
            queue1.put(None)

    monitor(len(workers), num_tasks, queue2)
    logger.info("done")


def insert_sites(uri: str, proteins_file: str, inqueue: mp.Queue,
                 batchsize: int, outqueue: mp.Queue):
    con = oracledb.connect(uri)
    cur1 = con.cursor()  # select
    cur2 = con.cursor()  # insert
    appls = get_i5_appls(cur1)

    with KVStore(proteins_file) as proteins:
        for start, stop, incl_stop in iter(inqueue.get, None):
            if incl_stop:
                where = "UPI BETWEEN :1 AND :2"
            else:
                where = "UPI >= :1 AND UPI < :2"

            cur1.execute(
                f"""
                SELECT UPI, ANALYSIS_ID, METHOD_AC, LOC_START, LOC_END, 
                       NUM_SITES, RESIDUE, RESIDUE_START, RESIDUE_END, 
                       DESCRIPTION
                FROM IPRSCAN.SITE
                WHERE {where}
                """,
                [start, stop]
            )

            while rows := cur1.fetchmany(batchsize):
                sites = []
                for row in rows:
                    _, _, md5 = proteins[row[0]]
                    dbname, dbversion = appls[row[1]]
                    sites.append((
                        md5,
                        dbname,
                        dbversion,
                        *row[2:]
                    ))

                cur2.executemany(
                    """
                    INSERT /*+ APPEND */ INTO IPRSCAN.LOOKUP_SITE
                    VALUES (:1, :2, :3, :4, :5, :6, :7, :8, :9, :10, :11)
                    """,
                    sites
                )
                con.commit()

            outqueue.put(False)

    cur1.close()
    cur2.close()
    con.close()
    outqueue.put(True)


def index_sites(uri: str):
    logger.info("indexing")
    con = oracledb.connect(uri)
    cur = con.cursor()
    cur.execute(
        """
        CREATE INDEX I_LOOKUP_SITE
        ON IPRSCAN.LOOKUP_SITE (MD5)
        NOLOGGING
        TABLESPACE IPRSCAN_IND
        """
    )
    cur.close()
    con.close()
    logger.info("done")


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


def monitor(num_workers: int, num_tasks: int, queue: mp.Queue):
    done = 0
    milestone = step = 5
    while num_workers:
        if queue.get():
            # A worker completed
            num_workers -= 1
        else:
            done += 1
            progress = done / num_tasks * 100
            if progress >= milestone:
                logger.info(f"\t{progress:.0f}%")
                milestone += step
