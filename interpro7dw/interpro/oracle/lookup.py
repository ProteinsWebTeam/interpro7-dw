import glob
import itertools
import os
import string
from multiprocessing import Process, Queue

import oracledb

from interpro7dw.utils import logger
from interpro7dw.utils.store import BasicStore


def drop_table(table_name: str, cur: oracledb.Cursor):
    try:
        cur.execute(f"DROP TABLE {table_name} PURGE")
    except oracledb.DatabaseError as exception:
        error_obj, = exception.args

        # ORA-00942: table or view does not exist
        # ORA-08103: object no longer exists
        if error_obj.code not in (942, 8103):
            raise exception


def create_md5_table(uri: str, indir: str, processes: int = 8):
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
    cur.close()
    con.close()

    workers = []
    inqueue = Queue()
    outqueue = Queue()
    for _ in range(processes):
        p = Process(target=insert_md5, args=(uri, inqueue, outqueue))
        p.start()
        workers.append(p)

    task_count = 0
    for filepath in glob.glob(os.path.join(indir, "*.dat")):
        inqueue.put(filepath)
        task_count += 1

    for _ in workers:
        inqueue.put(None)

    monitor(task_count, outqueue)

    logger.info("creating index")
    con = oracledb.connect(uri)
    cur = con.cursor()
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


def insert_md5(uri: str, inqueue: Queue, outqueue: Queue):
    con = oracledb.connect(uri)
    cur = con.cursor()
    statement = """
    INSERT /*+ APPEND */ INTO IPRSCAN.LOOKUP_MD5
    VALUES (:1)
    """
    records = []
    for filepath in iter(inqueue.get, None):
        with BasicStore(filepath) as bs:
            for proteins in bs:
                for protein in proteins.values():
                    records.append((protein["md5"],))

                    if len(records) == 10000:
                        cur.executemany(statement, records)
                        con.commit()
                        records.clear()

        outqueue.put(None)

    if records:
        cur.executemany(statement, records)
        con.commit()

    cur.close()
    con.close()


def create_matches_table(uri: str, indir: str, processes: int = 8):
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
    inqueue = Queue()
    outqueue = Queue()
    for _ in range(processes):
        p = Process(target=insert_matches, args=(uri, inqueue, outqueue))
        p.start()
        workers.append(p)

    task_count = 0
    for filepath in glob.glob(os.path.join(indir, "*.dat")):
        inqueue.put(filepath)
        task_count += 1

    for _ in workers:
        inqueue.put(None)

    monitor(task_count, outqueue)
    logger.info("done")


def insert_matches(uri: str, inqueue: Queue, outqueue: Queue):
    con = oracledb.connect(uri)
    cur = con.cursor()
    statement = """
    INSERT /*+ APPEND */ INTO IPRSCAN.LOOKUP_MATCH
    VALUES (:1, :2, :3, :4, :5, :6, :7, :8, :9, :10, :11, 
            :12, :13, :14, :15, :16, :17, :18, :19, :20)
    """
    records = []
    for filepath in iter(inqueue.get, None):
        with BasicStore(filepath) as bs:
            for proteins in bs:
                for protein in proteins.values():
                    for match in protein["matches"]:
                        signature = match["signature"]
                        library = get_i5_appl(signature["signatureLibraryRelease"]["library"])
                        for location in match["locations"]:
                            seq_evalue = match["evalue"]
                            seq_score = match["score"]
                            dom_evalue = location["evalue"]
                            dom_score = location["score"]

                            if library == "CDD":
                                seq_score = dom_score
                                seq_evalue = dom_evalue
                            elif library in ("HAMAP", "PRINTS", "PROSITE_PROFILES"):
                                seq_score = dom_score
                                dom_score = 0
                            elif library == "SUPERFAMILY":
                                seq_evalue = dom_evalue

                            records.append((
                                protein["md5"],
                                protein["md5"][:3],
                                library,
                                signature["signatureLibraryRelease"]["version"],
                                signature["accession"],
                                match["model-ac"],
                                location["start"],
                                location["end"],
                                location["extra"]["fragments"],
                                seq_score,
                                seq_evalue,
                                location["extra"]["hmm_bounds"],
                                location["hmmStart"],
                                location["hmmEnd"],
                                location["hmmLength"],
                                location["envelopeStart"],
                                location["envelopeEnd"],
                                dom_score,
                                dom_evalue,
                                location["sequence-feature"]
                            ))

                            if len(records) == 10000:
                                cur.executemany(statement, records)
                                con.commit()
                                records.clear()

        outqueue.put(None)

    if records:
        cur.executemany(statement, records)
        con.commit()

    cur.close()
    con.close()


def create_sites_table(uri: str, indir: str, processes: int = 8):
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
    inqueue = Queue()
    outqueue = Queue()
    for _ in range(processes):
        p = Process(target=insert_sites, args=(uri, inqueue, outqueue))
        p.start()
        workers.append(p)

    task_count = 0
    for filepath in glob.glob(os.path.join(indir, "*.dat")):
        inqueue.put(filepath)
        task_count += 1

    for _ in workers:
        inqueue.put(None)

    monitor(task_count, outqueue)
    logger.info("done")


def insert_sites(uri: str, inqueue: Queue, outqueue: Queue):
    con = oracledb.connect(uri)
    cur = con.cursor()
    statement = """
    INSERT /*+ APPEND */ INTO IPRSCAN.LOOKUP_SITE
    VALUES (:1, :2, :3, :4, :5, :6, :7, :8, :9, :10, :11, :12)
    """
    records = []
    for filepath in iter(inqueue.get, None):
        with BasicStore(filepath) as bs:
            for proteins in bs:
                records = []
                for protein in proteins.values():
                    for match in protein["matches"]:
                        signature = match["signature"]
                        for location in match["locations"]:
                            for site in location.get("sites", []):
                                for site_location in site["siteLocations"]:
                                    records.append((
                                        protein["md5"],
                                        protein["md5"][:3],
                                        get_i5_appl(signature["signatureLibraryRelease"]["library"]),
                                        signature["signatureLibraryRelease"]["version"],
                                        signature["accession"],
                                        location["start"],
                                        location["end"],
                                        site["numLocations"],
                                        site_location["residue"],
                                        site_location["start"],
                                        site_location["end"],
                                        site["description"]
                                    ))

                                    if len(records) == 10000:
                                        cur.executemany(statement, records)
                                        con.commit()
                                        records.clear()

        outqueue.put(None)

    if records:
        cur.executemany(statement, records)
        con.commit()

    cur.close()
    con.close()


def get_i5_appl(dbname: str) -> str:
    if dbname == "CATH-FunFam":
        return "FUNFAM"
    elif dbname == "CATH-Gene3D":
        return "GENE3D"

    return dbname.upper().replace(" ", "_")


def monitor(task_count: int, queue: Queue):
    milestone = step = 5
    for i in range(task_count):
        queue.get()
        progress = (i + 1) / task_count * 100
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
