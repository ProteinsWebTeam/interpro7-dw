import itertools
import string
import multiprocessing as mp

import oracledb

from interpro7dw.utils import logger
from interpro7dw.utils.store import KVStore


INSERT_SIZE = 10000


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

            if len(rows) == INSERT_SIZE:
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


def create_matches_table(uri: str, proteins_file: str, matches_file: str,
                         processes: int = 8):
    logger.info("starting")
    con = oracledb.connect(uri)
    cur = con.cursor()
    drop_table("IPRSCAN.LOOKUP_MATCH", cur)
    cur.execute(
        f"""
        CREATE TABLE IPRSCAN.LOOKUP_MATCH (
            MD5 CHAR(32) NOT NULL,
            MD5_PREFIX CHAR(3) NOT NULL,
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
        PARTITION BY LIST (MD5_PREFIX) (
            {','.join(get_partitions())}
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
                       args=(uri, proteins_file, matches_file, queue1, queue2))
        p.start()
        workers.append(p)

    with KVStore(proteins_file) as store:
        keys = store.get_keys()
        num_tasks = len(keys)

        for i, start in enumerate(keys):
            try:
                stop = keys[i + 1]
            except IndexError:
                stop = None
            finally:
                queue1.put((start, stop))

        for _ in workers:
            queue1.put(None)

    monitor(num_tasks, queue2)
    logger.info("done")


def insert_matches(uri: str, proteins_file: str, matches_file: str,
                   inqueue: mp.Queue, outqueue: mp.Queue):
    con = oracledb.connect(uri)
    cur = con.cursor()
    with KVStore(proteins_file) as proteins, KVStore(matches_file) as matches:
        for start, stop in iter(inqueue.get, None):
            for upi, (_, _, md5) in proteins.range(start, stop):
                records = []
                for match in matches.get(upi, {}).values():
                    for (loc_start, loc_end, hmm_start, hmm_end, hmm_length,
                         hmm_bounds, dom_evalue, dom_score, env_start, env_end,
                         fragments, seq_feature) in match["locations"]:
                        records.append((
                            md5,
                            md5[:3],
                            get_i5_appl(match["database"]["name"]),
                            match["database"]["version"],
                            match["signature"]["accession"],
                            match["model"],
                            loc_start,
                            loc_end,
                            fragments,
                            match["score"],
                            match["evalue"],
                            hmm_bounds,
                            hmm_start,
                            hmm_end,
                            hmm_length,
                            env_start,
                            env_end,
                            dom_score,
                            dom_evalue,
                            seq_feature
                        ))

                for i in range(0, len(records), INSERT_SIZE):
                    cur.executemany(
                        """
                        INSERT /*+ APPEND */ INTO IPRSCAN.LOOKUP_MATCH
                        VALUES (:1, :2, :3, :4, :5, :6, :7, :8, :9, :10, :11, 
                                :12, :13, :14, :15, :16, :17, :18, :19, :20)
                        """,
                        records[i:i + INSERT_SIZE]
                    )
                    con.commit()

            outqueue.put(None)

    cur.close()
    con.close()


def create_sites_table(uri: str, proteins_file: str, sites_file: str,
                       processes: int = 8):
    logger.info("starting")
    con = oracledb.connect(uri)
    cur = con.cursor()
    drop_table("IPRSCAN.LOOKUP_SITE", cur)
    cur.execute(
        f"""
        CREATE TABLE IPRSCAN.LOOKUP_SITE (
            MD5 CHAR(32) NOT NULL,
            MD5_PREFIX CHAR(3) NOT NULL,
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
        PARTITION BY LIST (MD5_PREFIX) (
            {','.join(get_partitions())}
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
                       args=(uri, proteins_file, sites_file, queue1, queue2))
        p.start()
        workers.append(p)

    with KVStore(proteins_file) as store:
        keys = store.get_keys()
        num_tasks = len(keys)

        for i, start in enumerate(keys):
            try:
                stop = keys[i + 1]
            except IndexError:
                stop = None
            finally:
                queue1.put((start, stop))

        for _ in workers:
            queue1.put(None)

    monitor(num_tasks, queue2)
    logger.info("done")


def insert_sites(uri: str, proteins_file: str, sites_file: str,
                 inqueue: mp.Queue, outqueue: mp.Queue):
    con = oracledb.connect(uri)
    cur = con.cursor()

    with KVStore(proteins_file) as proteins, KVStore(sites_file) as sites:
        for start, stop in iter(inqueue.get, None):
            for upi, (_, _, md5) in proteins.range(start, stop):
                records = []
                for sig_acc, signature in sites.get(upi, {}).items():
                    for (loc_start, loc_end), descriptions in signature.items():
                        for description, site_locations in descriptions.items():
                            for site in site_locations:
                                records.append((
                                    md5,
                                    md5[:3],
                                    signature["database"]["name"],
                                    signature["database"]["version"],
                                    sig_acc,
                                    loc_start,
                                    loc_end,
                                    len(site_locations),
                                    site["residue"],
                                    site["start"],
                                    site["end"],
                                    description
                                ))

                for i in range(0, len(records), INSERT_SIZE):
                    cur.executemany(
                        """
                        INSERT /*+ APPEND */ INTO IPRSCAN.LOOKUP_SITE
                        VALUES (:1, :2, :3, :4, :5, :6, :7, :8, :9, :10, :11, :12)
                        """,
                        records[i:i + INSERT_SIZE]
                    )
                    con.commit()

            outqueue.put(False)

    cur.close()
    con.close()


def get_i5_appl(dbname: str) -> str:
    if dbname == "CATH-Gene3D":
        return "GENE3D"

    return dbname.upper().replace(" ", "_")


def monitor(num_tasks: int, queue: mp.Queue):
    milestone = step = 5
    for i in range(num_tasks):
        queue.get()
        progress = (i + 1) / num_tasks * 100
        if progress >= milestone:
            logger.info(f"\t{progress:.0f}%")
            milestone += step


def get_partitions() -> list[str]:
    partitions = []
    # Get hexadecimal characters (upper case only), naturally ordered
    hex_chars = sorted(set(string.hexdigits.upper()))

    # Create 4096 possible combinations with three characters (16*16*16)
    for chars in itertools.product(hex_chars, repeat=3):
        name = str(len(partitions)).zfill(4)
        value = "".join(chars)
        partitions.append(f"PARTITION P{name} VALUES ('{value}')")

    return partitions
