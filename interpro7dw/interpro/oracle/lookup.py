import itertools

from oracledb import connect, Cursor, DatabaseError

from interpro7dw.utils import logger


def drop_table(table_name: str, cur: Cursor):
    try:
        cur.execute(f"DROP TABLE {table_name}")
    except DatabaseError as exception:
        error_obj, = exception.args

        # ORA-00942: table or view does not exist
        # ORA-08103: object no longer exists
        if error_obj.code not in (942, 8103):
            raise exception


def get_partitions():
    filter_chars = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E', 'F']
    filter_chars_product = itertools.product(filter_chars, repeat=3)
    partitions = []
    upi_range_base = 'UPI00'

    for el in filter_chars_product:
        el_str = ''.join(el)
        partition_name = upi_range_base + el_str
        partitions.append(partition_name)

    return partitions


def build_upi_md5_table(ipr_uri: str):
    con = connect(ipr_uri)
    cur = con.cursor()

    drop_table('lookup_tmp_upi_md5', cur)
    cur.execute(
        """
        CREATE TABLE lookup_tmp_upi_md5 NOLOGGING AS
        SELECT upi, md5 FROM uniparc.protein
        WHERE 1 = 0
        """
    )

    cur.execute(
        """
        INSERT /*+ APPEND */ INTO lookup_tmp_upi_md5
        SELECT upi, md5
        FROM uniparc.protein
        """
    )

    cur.execute(
        """
        CREATE INDEX lookup_upi_UPIX 
        ON lookup_tmp_upi_md5(UPI)
        TABLESPACE IPRSCAN_IND
        """
    )

    con.commit()

    cur.close()
    con.close()


def build_matches_table(ipr_uri: str):
    con = connect(ipr_uri)
    cur = con.cursor()

    drop_table('db_versions_tmp_tab', cur)
    # Create new db_versions_tmp_tab table (the list of analysis IDs to be in the Berkeley DB build)
    cur.execute(
        """
        CREATE TABLE db_versions_tmp_tab AS
        SELECT r.iprscan_sig_lib_rel_id, 
            DECODE(DBNAME, 'CATH-Gene3D', 'GENE3D', UPPER(REPLACE(DBNAME, ' ', '_'))) LIBRARY, 
            d.VERSION
        FROM INTERPRO.iprscan2dbcode r
        INNER JOIN INTERPRO.CV_DATABASE c
        ON r.DBCODE = c.DBCODE
        INNER JOIN INTERPRO.DB_VERSION d
        ON c.DBCODE = d.DBCODE
        """
    )

    cur.execute(
        """
        ALTER TABLE db_versions_tmp_tab 
        ADD CONSTRAINT DB_VERSIONS_TMP_TAB_PK PRIMARY KEY(iprscan_sig_lib_rel_id)
        """
    )

    drop_table('lookup_tmp_tab', cur)

    cur.execute(
        """
        SELECT *
        FROM db_versions_tmp_tab
        """
    )
    analyses = cur.fetchall()
    if not analyses:
        raise RuntimeError("No analyses found")

    partitions = get_partitions()

    sql = """
        CREATE TABLE lookup_tmp_tab (
            ID NUMBER,
            UPI_RANGE VARCHAR2(8),
            ANALYSIS_ID NUMBER(19,0),
            PROTEIN_MD5 VARCHAR2(32) NOT NULL,
            SIGNATURE_LIBRARY_NAME VARCHAR2(25),
            SIGNATURE_LIBRARY_RELEASE VARCHAR2(20) NOT NULL,
            SIGNATURE_ACCESSION VARCHAR2(255),
            MODEL_ACCESSION VARCHAR2(255),
            SCORE BINARY_DOUBLE,
            SEQUENCE_SCORE BINARY_DOUBLE,
            SEQUENCE_EVALUE BINARY_DOUBLE,
            EVALUE BINARY_DOUBLE,
            SEQ_START NUMBER(10,0),
            SEQ_END NUMBER(10,0),
            HMM_START NUMBER,
            HMM_END NUMBER,
            HMM_LENGTH NUMBER,
            HMM_BOUNDS VARCHAR2(25),
            ENVELOPE_START NUMBER,
            ENVELOPE_END NUMBER,
            SEQ_FEATURE VARCHAR2(4000),
            FRAGMENTS VARCHAR2(400)
        ) PARTITION BY LIST (UPI_RANGE) (
        """

    for i, value in enumerate(sorted(partitions)):
        if i:
            sql += ', '
        sql += f"partition {value} VALUES('{value}')"
    sql += ", partition OTHER VALUES(default))"

    cur.execute(sql)

    for progress_count, analysis in enumerate(analyses, start=1):
        cur.execute(
            """
            INSERT /*+ APPEND */ INTO lookup_tmp_tab nologging
            SELECT rownum as id,
                substr(m.upi,0,8) upi_range,
                m.analysis_id,
                p.md5 as protein_md5,
                cast(:name as VARCHAR2(255 CHAR)) as signature_library_name,
                cast(:version as VARCHAR2(255 CHAR)) as signature_library_release,
                m.method_ac as signature_accession,
                m.model_ac as model_accession,
                m.score as score,
                m.seqscore as sequence_score,
                m.seqevalue as sequence_evalue,
                m.evalue,
                m.seq_start,
                m.seq_end,
                m.hmm_start,
                m.hmm_end,
                m.hmm_length,
                m.hmm_bounds,
                m.envelope_start,
                m.envelope_end,
                m.seq_feature,
                m.fragments
            FROM lookup_tmp_upi_md5 p, mv_iprscan m
            WHERE m.upi = p.upi
            AND m.analysis_id = :id
            """, [analysis[1], analysis[2], analysis[0]]
        )
        con.commit()

        logger.info(f"Processed {progress_count} of {len(analyses)}")

    cur.execute(
        """
        CREATE INDEX LKP_RANGE_MD5X
        ON lookup_tmp_tab(upi_range,protein_md5)
        TABLESPACE IPRSCAN_IND
        """
    )

    con.commit()

    cur.close()
    con.close()


def build_site_table(ipr_uri: str):
    con = connect(ipr_uri)
    cur = con.cursor()

    drop_table('db_versions_site_tmp_tab', cur)

    # Create new db_versions_tmp_tab table (the list of analysis IDs to be in the Berkeley DB build)
    cur.execute(
        """
        CREATE TABLE db_versions_site_tmp_tab AS
        SELECT iprscan_sig_lib_rel_id, library, version FROM (
            SELECT r.iprscan_sig_lib_rel_id, 
                DECODE(DBNAME, 'CATH-Gene3D', 'GENE3D', UPPER(REPLACE(DBNAME, ' ', '_'))) LIBRARY, 
                d.VERSION
            FROM INTERPRO.iprscan2dbcode r
            INNER JOIN INTERPRO.CV_DATABASE c
            ON r.DBCODE = c.DBCODE
            INNER JOIN INTERPRO.DB_VERSION d
            ON c.DBCODE = d.DBCODE
        ) WHERE LIBRARY in ('SFLD', 'CDD', 'PIRSR')
        """
    )

    cur.execute(
        """
        ALTER TABLE db_versions_site_tmp_tab
        ADD CONSTRAINT DB_VERSIONS_STMP_TAB_PK PRIMARY KEY(iprscan_sig_lib_rel_id)
        """
    )

    drop_table('lookup_site_tmp_tab', cur)

    cur.execute(
        """
        SELECT *
        FROM db_versions_site_tmp_tab
        """
    )
    analyses = cur.fetchall()
    if not analyses:
        raise RuntimeError("No analyses found")

    partitions = get_partitions()

    sql = """
        CREATE TABLE lookup_site_tmp_tab (
            ID NUMBER,
            UPI_RANGE VARCHAR2(8),
            ANALYSIS_ID NUMBER(19,0),
            PROTEIN_MD5 VARCHAR2(32),
            SIGNATURE_LIBRARY_NAME VARCHAR2(25),
            SIGNATURE_LIBRARY_RELEASE VARCHAR2(20),
            SIGNATURE_ACCESSION VARCHAR2(255) NOT NULL,
            LOC_START NUMBER(10,0) NOT NULL,
            LOC_END NUMBER(10,0) NOT NULL,
            NUM_SITES NUMBER,
            RESIDUE VARCHAR2(1000),
            RESIDUE_START NUMBER,
            RESIDUE_END NUMBER,
            DESCRIPTION VARCHAR2(255)
        ) PARTITION BY LIST (UPI_RANGE) (
        """

    for i, value in enumerate(sorted(partitions)):
        if i:
            sql += ', '
        sql += f"partition {value} VALUES('{value}')"
    sql += ", partition OTHER VALUES(default))"

    cur.execute(sql)

    for progress_count, analysis in enumerate(analyses, start=1):
        for upi_range_partition in sorted(partitions):
            cur.execute(
                """
                INSERT /*+ APPEND */ INTO lookup_site_tmp_tab nologging
                SELECT rownum as id,
                    s.upi_range,
                    s.analysis_id,
                    p.md5 as protein_md5,
                    cast(:name as VARCHAR2(255 CHAR)) as signature_library_name,
                    cast(:version as VARCHAR2(255 CHAR)) as signature_library_release,
                    s.method_ac as signature_accession,
                    s.loc_start,
                    s.loc_end,
                    s.num_sites,
                    s.residue,
                    s.residue_start,
                    s.residue_end,
                    s.description
                FROM lookup_tmp_upi_md5 p, site s
                WHERE s.upi = p.upi and s.upi_range = :partition
                AND s.analysis_id = :id
                """,
                [analysis[1], analysis[2], upi_range_partition, analysis[0]]
            )
            con.commit()

        logger.info(f"Processed {progress_count} of {len(analyses)}")

    cur.execute(
        """
        CREATE INDEX LKP_SITE_RANGE_MD5X
        ON lookup_site_tmp_tab(upi_range,protein_md5)
        TABLESPACE IPRSCAN_IND
        """
    )

    con.commit()
    cur.close()
    con.close()
