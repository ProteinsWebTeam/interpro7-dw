import itertools

import cx_Oracle

from interpro7dw.utils import logger


def drop_table(table_name, cur):
    try:
        cur.execute(f"DROP TABLE {table_name}")
        logger.info(f"Dropped table: {table_name}")
    except cx_Oracle.DatabaseError as exception:
        error_obj, = exception.args
        err_code = error_obj.code

        if err_code == 942:
            logger.warning(f"Could not drop table {table_name} as it did not exist, continuing anyway")
        else:
            logger.error(str(exception))
            raise Exception(f"Failed to drop table {table_name}")


def get_partitions():
    filter_chars = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E', 'F']
    filter_chars_product = itertools.product(filter_chars, repeat=2)
    partitions = []
    upi_range_base = 'UPI000'
    upi_range_base_2ndtier = 'UPI001'
    for el in filter_chars_product:
        el_str = ''.join(el)
        partition_name = upi_range_base + el_str
        partitions.append(partition_name)

        if el_str[0] in '0123456789':
            partition_name = upi_range_base_2ndtier + el_str
            partitions.append(partition_name)

    return partitions


def build_upi_md5_tbl(ipr_uri: str):
    con = cx_Oracle.connect(ipr_uri)
    cur = con.cursor()

    cur.execute(
        """
        SELECT MAX(UPI)
        FROM UNIPARC.PROTEIN
        """
    )

    row = cur.fetchone()
    maxupi = row[0]
    logger.info(f"MAXUPI: {maxupi}")

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
        WHERE upi <= :maxupi
        """, [maxupi]
    )

    cur.execute(
        """
        CREATE INDEX lookup_upi_UPIX ON lookup_tmp_upi_md5(UPI)
        """
    )

    con.commit()

    cur.close()
    con.close()

    logger.info("lookup_tmp_upi_md5 table built.")


def build_lookup_tmp_tab(ipr_uri: str):
    con = cx_Oracle.connect(ipr_uri)
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
        ALTER TABLE db_versions_tmp_tab add constraint DB_VERSIONS_TMP_TAB_PK primary key(iprscan_sig_lib_rel_id)
        """
    )

    drop_table('lookup_tmp_tab', cur)
    # Get list of analysis IDs that will be included in the Berkley DB build
    cur.execute(
        """
        SELECT *
        FROM db_versions_tmp_tab
        """
    )

    analyses = cur.fetchall()
    logger.info(f"analysis: {str(analyses)}")

    if not analyses:
        raise Exception("No analyses found in the iprscan.db_versions_tmp_tab table")

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
        ) partition BY list (upi_range) (
        """

    first_partition = True
    for partition_value in sorted(partitions):
        if first_partition:
            first_partition = False
        else:
            sql += ', '
        sql += f"partition {partition_value} VALUES('{partition_value}')"
    sql += ", partition OTHER VALUES(default))"

    cur.execute(sql)

    analysis_count = len(analyses)
    logger.info(f"{str(analysis_count)} analyses to process")

    for progress_count, analysis in enumerate(analyses, start=1):

        logger.info(f"Prepare data for analysis: {str(analysis[0])} - {analysis[1]} ({analysis[2]})")

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

        logger.info(f"Processed {str(progress_count)} of {str(analysis_count)} ...")

    con.commit()

    cur.close()
    con.close()

    logger.info("lookup_tmp_tab table data has been populated.")


def build_lookup_tmp_tab_idx(ipr_uri: str):
    con = cx_Oracle.connect(ipr_uri)
    cur = con.cursor()

    cur.execute(
        """
        CREATE INDEX LKP_RANGE_MD5X
        ON lookup_tmp_tab(upi_range,protein_md5)
        """
    )

    con.commit()

    cur.close()
    con.close()

    logger.info("lookup_tmp_tab table index have been created.")


def build_site_lookup_tmp_tab(ipr_uri: str):
    con = cx_Oracle.connect(ipr_uri)
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

    # Add indexes to db_versions_tmp_tab table
    cur.execute(
        """
        ALTER TABLE db_versions_site_tmp_tab
        ADD constraint DB_VERSIONS_STMP_TAB_PK primary key(iprscan_sig_lib_rel_id)
        """
    )

    drop_table('lookup_site_tmp_tab', cur)

    # Get list of analysis IDs that will be included in the Berkley DB build
    cur.execute(
        """
        SELECT *
        FROM db_versions_site_tmp_tab
        """
    )

    analyses = cur.fetchall()

    logger.info(f"analysis: {str(analyses)}")

    if not analyses:
        raise Exception("No analyses found in the iprscan.db_versions_tmp_tab table")

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
        ) partition BY list (upi_range) (
        """

    first_partition = True
    for partition_value in sorted(partitions):
        if first_partition:
            first_partition = False
        else:
            sql += ', '
        sql += f"partition {partition_value} VALUES('{partition_value}')"
    sql += ", partition OTHER VALUES(default))"

    cur.execute(sql)

    analysis_count = len(analyses)
    logger.info(f"{str(analysis_count)} analyses to process")

    for progress_count, analysis in enumerate(analyses, start=1):

        logger.info(f"Prepare data for analysis: {str(analysis[0])} - {analysis[1]} ({analysis[2]})")

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
                """, [analysis[1], analysis[2], upi_range_partition, analysis[0]]
            )
            con.commit()
            logger.debug(f"Processed range {upi_range_partition} : {str(progress_count)} of {str(analysis_count)} ...")

        logger.info(f"Processed {str(progress_count)} of {str(analysis_count)} ...")

    con.commit()

    cur.close()
    con.close()

    logger.info("lookup_site_tmp_tab table data has been populated.")


def build_site_lookup_tmp_tab_idx(ipr_uri: str):
    con = cx_Oracle.connect(ipr_uri)
    cur = con.cursor()

    cur.execute(
        """
        CREATE INDEX LKP_SITE_RANGE_MD5X
        ON lookup_site_tmp_tab(upi_range,protein_md5)
        """
    )

    con.commit()

    cur.close()
    con.close()

    logger.info("lookup_site_tmp_tab table index have been created.")
