from interpro7dw.utils import logger, oracle
import re
import cx_Oracle
import sys
import traceback
import itertools



def get_maxupi(ipr_uri: str):

    logger.info("retrieving maxupi")
    con = cx_Oracle.connect(ipr_uri)
    cur = con.cursor()


    # upi_md5 table creation

    cur.execute(
        """
        SELECT MAX(UPI)
        FROM UNIPARC.PROTEIN
        """
    )

    row = cur.fetchone()
    maxupi = row[0]

    cur.close()
    con.close()

    return maxupi


def drop_table(table_name, cur):
    try:
        cur.execute( "DROP TABLE {0}".format(table_name) )
        logger.debug('Dropped table: ' + table_name)
    except cx_Oracle.DatabaseError as exception:
        if 'table or view does not exist' in str(exception):
            logger.warning('Could not drop table ' + table_name + ' as it did not exist, continuing anyway')
        else:
            logger.error('Failed to execute: ' + "DROP TABLE {0}".format(table_name))
            logger.debug(exception)



def get_partitions():

    # Create partitions array

    # First character of the protein_md5 to filter on ('G' is OK as it's partitioned  with less than)
    md5_filter_chars = ['1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E', 'F', 'G']

    #get the partition_list
    filter_chars = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E', 'F']
    filter_chars_product = itertools.product(filter_chars, repeat=2)
    partitions = []
    upi_range_base = 'UPI000'
    upi_range_base_2ndtier = 'UPI001'
    for el in filter_chars_product:
      el_str = ''.join(el)
      partition_name = upi_range_base + el_str
      partitions.append(partition_name)
      #if el_str.startswith('0') or el_str.startswith('1') or el_str.startswith('2') or el_str.startswith('3'):
      if el_str[0] in '0123456789':
        partition_name = upi_range_base_2ndtier + el_str
        partitions.append(partition_name)

    return partitions



def build_upi_md5_tbl(ipr_uri: str, maxupi: str):
    logger.setLevel('DEBUG')

    logger.info("Preparing to create upi_md5 table with max_upi: " + maxupi)
    con = cx_Oracle.connect(ipr_uri)
    cur = con.cursor()

    logger.info("Connected to " + ipr_uri)

    # upi_md5 table creation

    cur.execute(
        """
        ALTER SESSION enable parallel dml
        """
    )

    tmp_upi_md5_table = 'lookup_tmp_upi_md5'
    drop_table(tmp_upi_md5_table, cur)

    cur.execute(
        """
        CREATE TABLE {0} NOLOGGING AS
        SELECT /*+ PARALLEL */ upi, md5
        FROM iprscan.uniparc_protein
        WHERE upi <= '{1}'
        """.format(tmp_upi_md5_table, maxupi)
    )

    cur.execute(
        """
        CREATE INDEX lookup_upi_UPIX ON {0}(UPI) TABLESPACE IPRSCAN_IND NOLOGGING PARALLEL 5
        """.format(tmp_upi_md5_table)
    )

    con.commit()

    cur.close()
    con.close()

    logger.info('Done')






def build_lookup_tmp_tab(ipr_uri: str, maxupi: str):
    logger.setLevel('DEBUG')

    logger.info("Preparing to build lookup tmp tables")
    con = cx_Oracle.connect(ipr_uri)
    cur = con.cursor()

    logger.info("Connected to " + ipr_uri)

    # Build lookup_tmp_tab table

    lookup_table = 'lookup_tmp_tab'
    logger.info("preparing to build table " + lookup_table)

    # Drop db_versions_tmp_tab table
    db_versions_table = 'db_versions_tmp_tab'
    drop_table(db_versions_table, cur)

    # Create new db_versions_tmp_tab table (the list of analysis IDs to be in the Berkeley DB build)
    cur.execute(
        """
        CREATE TABLE {0} AS
        SELECT r.iprscan_sig_lib_rel_id, l.library, l.version
        FROM INTERPRO.iprscan2dbcode r
        INNER JOIN iprscan.ipm_signature_library_release@ispro l
        ON (r.iprscan_sig_lib_rel_id=l.id)
        """.format(db_versions_table)
    )

    # Add indexes to db_versions_tmp_tab table
    cur.execute(
        """
        ALTER TABLE {0} add constraint DB_VERSIONS_TMP_TAB_PK primary key(iprscan_sig_lib_rel_id)
        """.format(db_versions_table)
    )

    drop_table(lookup_table, cur)

    # Get list of analysis IDs that will be included in the Berkley DB build
    cur.execute(
        """
        SELECT *
        FROM {0}
        """.format(db_versions_table)
    )

    analyses = cur.fetchall()

    logger.info("analysis: " + str(analyses))


    if not analyses:
        logger.error("No analyses found in the iprscan.db_versions_tmp_tab table")
        exit (1)


    # # Create partitions array
    partitions = get_partitions()


    # Create empty lookup_tmp_tab table
    sql = 'CREATE TABLE ' + lookup_table + ' partition BY list (upi_range) ('
    first_partition = True
    for partition_value in sorted(partitions):
        if first_partition:
            first_partition = False
        else:
            sql += ', '
        sql += 'partition ' + partition_value + " VALUES ('" + partition_value + "')"
    sql += ", partition OTHER VALUES(default)"


    sql += """)
           AS (SELECT  /* PARALLEL */  rownum AS id,
                                       substr(m.upi,0,8) upi_range,
                                       m.analysis_id,
                                       p.md5 AS protein_md5,
                                       v.library AS signature_library_name,
                                       v.version AS signature_library_release,
                                       m.method_ac AS signature_accession,
                                       m.model_ac AS model_accession,
                                       m.score AS score,
                                       m.seqscore AS sequence_score,
                                       m.seqevalue AS sequence_evalue,
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
                                  FROM lookup_tmp_upi_md5 p,
                                       mv_iprscan m,
                                       db_versions_tmp_tab v
                                 WHERE 1=0)"""

    cur.execute(sql)


    analysis_count_dict = {}

    analysis_count = len(analyses)

    progress_count = 0
    logger.info(str(analysis_count) + ' analyses to process')
    for analysis in analyses:
        progress_count = progress_count + 1

        logger.info('Prepare data for analysis: ' + str(analysis[0]) + ' - ' + analysis[1] + ' (' + analysis[2] + ')')


        sql = 'insert into ' + lookup_table + ' nologging '
        sql += """ SELECT  /*+ PARALLEL */  rownum as id,
                     substr(m.upi,0,8) upi_range,
                     m.analysis_id,
                     p.md5 as protein_md5,
                     cast('{0}' as VARCHAR2(255 CHAR)) as signature_library_name,
                     cast('{1}' as VARCHAR2(255 CHAR)) as signature_library_release,
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
                FROM lookup_tmp_upi_md5 p,
                     mv_iprscan m
               WHERE m.upi = p.upi
                 AND m.analysis_id={2}""".format(analysis[1], analysis[2], analysis[0])

        cur.execute(sql)
        con.commit()

        logger.info("Processed " + str(progress_count) + " of " + str(analysis_count) + " ...")

    logger.info(lookup_table + ' table data has been populated.')

    con.commit()

    cur.close()
    con.close()

    logger.info('Done')



def build_lookup_tmp_tab_idx(ipr_uri: str, maxupi: str):
    logger.setLevel('DEBUG')

    logger.info("Preparing to build lookup tmp tables index")
    con = cx_Oracle.connect(ipr_uri)
    cur = con.cursor()

    logger.info("Connected to " + ipr_uri)

    lookup_table = 'lookup_tmp_tab'

    partitions = get_partitions()

    #create the local indices

    logger.info('Creating indices for ' + lookup_table)

    lookup_main_index = 'LKP_RANGE_MD5X'
    sql_idx = 'CREATE INDEX ' + lookup_main_index + ' ON ' + lookup_table + '(upi_range,protein_md5) LOCAL TABLESPACE IPRSCAN_IND NOLOGGING unusable'
    cur.execute(sql_idx)

    #now rebuild each index in each partition

    for partition_value in partitions:
        sql_rebuild_idx = 'alter index ' + lookup_main_index + ' rebuild partition ' + partition_value + ' PARALLEL 5'
        cur.execute(sql_rebuild_idx)

    logger.info('Completed partition ' + partition_value + '...')

    logger.info(lookup_table + ' table indices have been created.')


    con.commit()

    cur.close()
    con.close()

    logger.info('Done')





def build_site_lookup_tmp_tab(ipr_uri: str, maxupi: str):
    logger.setLevel('DEBUG')

    logger.info("preparing to built site lookup tmp tables")
    con = cx_Oracle.connect(ipr_uri)
    cur = con.cursor()

    logger.info("Connected to " + ipr_uri)


    # Do the site matches now

    lookup_table = 'lookup_site_tmp_tab'

    logger.info("Preparing to build table " + lookup_table)

    # Drop db_versions_site_tmp_tab table
    db_versions_table = 'db_versions_site_tmp_tab'

    # Drop db_versions_tmp_tab table
    drop_table(db_versions_table, cur)

    # Create new db_versions_tmp_tab table (the list of analysis IDs to be in the Berkeley DB build)
    # TODO Remove SFLD and CDD exclusion once their sites are also included in the Berkeley DB - until then they will always need to be calculated locally anyway, therefore there is no point including them in lookup service!
    sql = """
        CREATE TABLE {0} AS
        SELECT r.iprscan_sig_lib_rel_id, l.library, l.version
        FROM INTERPRO.iprscan2dbcode r inner join iprscan.mv_signature_library_release l on (r.iprscan_sig_lib_rel_id=l.id)
        AND library in ('SFLD', 'CDD')
        """.format(db_versions_table)

    cur.execute(sql)

    # Add indexes to db_versions_tmp_tab table
    sql = """
        ALTER TABLE {0}
        ADD constraint DB_VERSIONS_STMP_TAB_PK primary key(iprscan_sig_lib_rel_id)
        """.format(db_versions_table)

    cur.execute(sql)

    # Drop lookup_tmp_tab table
    drop_table(lookup_table, cur)

    # Get list of analysis IDs that will be included in the Berkley DB build
    cur.execute(
        """
        SELECT *
        FROM {0}
        """.format(db_versions_table)
    )

    analyses = cur.fetchall()

    logger.info("analysis: " + str(analyses))


    #analyses = [[37,'PROSITE_PROFILES','20.119'], [50,'SFLD','2']] # Small example test case commented out!
    if not analyses:
        logger.error("No analyses found in the iprscan.db_versions_tmp_tab table")
        exit (1)


    # Create partitions array
    partitions = get_partitions()

    # Create empty lookup_tmp_tab table
    sql = 'create table ' + lookup_table + ' partition BY list (upi_range) ('
    first_partition = True
    for partition_value in sorted(partitions):
        if first_partition:
            first_partition = False
        else:
            sql += ', '
        sql += 'partition ' + partition_value + " VALUES('" + partition_value + "')"
    sql += ", partition OTHER VALUES(default)"


    sql += """)
           as (SELECT  /* PARALLEL */  rownum as id,
                           s.upi_range,
                           s.analysis_id,
                           p.md5 as protein_md5,
                           v.library as signature_library_name,
                           v.version as signature_library_release,
                           s.method_ac as signature_accession,
                           s.loc_start,
                           s.loc_end,
                           s.num_sites,
                           s.residue,
                           s.residue_start,
                           s.residue_end,
                           s.description
                      FROM lookup_tmp_upi_md5 p,
                           site s,
                           db_versions_site_tmp_tab v
                     WHERE 1=0)"""
    cur.execute(sql)


    # Note: lookup_tmp_tab and lookup_tmp_tab_partition tables will need to have identical column structures for
    # partition exchange to succeed.
    # These two columns should not really be null but this is done for speed (because constraints on columns won't be
    # transferred to lookup_tmp_tab_partition table when using "cast('STRING' as VARCHAR2(255 CHAR))" in the CTAS as we
    # do below).
    # We could add the "not null" constraints to lookup_tmp_tab table at the end if we want. TODO review later?
    sql = 'alter table ' + lookup_table + ' modify (signature_library_name null)'
    cur.execute(sql)
    sql = 'alter table ' + lookup_table + ' modify (signature_library_release null)'
    cur.execute(sql)
    sql = 'alter table ' + lookup_table + ' modify (protein_md5 null)'
    cur.execute(sql)

    analysis_count_dict = {}

    analysis_count = len(analyses)

    progress_count = 0
    logger.info(str(analysis_count) + ' analyses to process')
    #where s.upi = p.upi and substr(s.upi,0,8) = '{3}' AND s.analysis_id={2}""".format(analysis[1], analysis[2], analysis[0], upi_range_partition)
    for analysis in analyses:
        progress_count = progress_count + 1

        logger.info('Prepare data for analysis: ' + str(analysis[0]) + ' - ' + analysis[1] + ' (' + analysis[2] + ')')


        for upi_range_partition in  sorted(partitions):
           sql = 'INSERT INTO ' + lookup_table + ' nologging '
           sql += """ SELECT  /*+ PARALLEL (2) */  rownum as id,
                    s.upi_range,
                    s.analysis_id,
                    p.md5 as protein_md5,
                    cast('{0}' as VARCHAR2(255 CHAR)) as signature_library_name,
                    cast('{1}' as VARCHAR2(255 CHAR)) as signature_library_release,
                    s.method_ac as signature_accession,
                    s.loc_start,
                    s.loc_end,
                    s.num_sites,
                    s.residue,
                    s.residue_start,
                    s.residue_end,
                    s.description
                FROM lookup_tmp_upi_md5 p,
                     site s
                WHERE s.upi = p.upi and s.upi_range = '{3}'
                AND s.analysis_id={2}""".format(analysis[1], analysis[2], analysis[0], upi_range_partition)

           cur.execute(sql)
           con.commit()
           logger.debug("Processed range " + upi_range_partition + ": " + str(progress_count) + " of " + str(analysis_count) + " ...")

        logger.info("Processed " + str(progress_count) + " of " + str(analysis_count) + " ...")

    logger.info(lookup_table + ' table data has been populated.')


    con.commit()

    cur.close()
    con.close()

    logger.info('Done')



def build_site_lookup_tmp_tab_idx(ipr_uri: str, maxupi: str):
    logger.setLevel('DEBUG')

    logger.info("preparing to built site lookup tmp tables")
    con = cx_Oracle.connect(ipr_uri)
    cur = con.cursor()

    logger.info("Connected to " + ipr_uri)


    # Do the site matches now

    lookup_table = 'lookup_site_tmp_tab'

    partitions = get_partitions()


    #create the local indices

    logger.info('Creating indices for ' + lookup_table)

    lookup_main_index = 'LKP_SITE_RANGE_MD5X'
    sql_idx = "CREATE INDEX {0} ON {1}(upi_range,protein_md5) LOCAL TABLESPACE IPRSCAN_IND NOLOGGING unusable".format(lookup_main_index, lookup_table)
    cur.execute(sql_idx)
    #now rebuild each index in each partition
    for partition_value in partitions:
        sql_rebuild_idx = "ALTER INDEX {0} rebuild partition {1} PARALLEL 5".format(lookup_main_index, partition_value)
        cur.execute(sql_rebuild_idx)
    logger.debug('Completed partition ' + partition_value + '...')



    logger.info(lookup_table + ' table indices have been created.')

    con.commit()

    cur.close()
    con.close()

    logger.info('Done')

