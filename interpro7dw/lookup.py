#import logging
from interpro7dw.utils import logger, oracle
#import time, datetime
import re
import cx_Oracle
import sys
import traceback
import itertools





def build_tmp_tables(ipr_uri: str):


    logger.info("preparing to built lookup tmp tables")
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
    logger.info("upi_md5 table creation, MAXUPI: " + maxupi)

    cur.execute(
        """
        ALTER SESSION enable parallel dml
        """
    )

    try:
        cur.execute(
            """
            DROP TABLE lookup_tmp_upi_md5
            """
        )
    except cx_Oracle.DatabaseError as exception:
            if 'table or view does not exist' in str(exception):
                logger.warning('Could not drop table as it did not exist, continuing anyway')

    cur.execute(
        """
        CREATE TABLE lookup_tmp_upi_md5 NOLOGGING AS
        SELECT /*+ PARALLEL */ upi, md5
        FROM iprscan.uniparc_protein
        WHERE upi <= '" + maxupi + "'
        """
    )

    cur.execute(
        """
        CREATE INDEX lookup_upi_UPIX ON lookup_tmp_upi_md5(UPI) TABLESPACE IPRSCAN_IND NOLOGGING PARALLEL 5
        """
    )

    # con.commit()


    # Build lookup_tmp_tab table

    lookup_table = 'lookup_tmp_tab'
    logger.info("preparing to build table " + lookup_table)

    # Drop db_versions_tmp_tab table
    cur.execute(
        """
        DROP TABLE db_versions_tmp_tab
        """
    )

    # Create new db_versions_tmp_tab table (the list of analysis IDs to be in the Berkeley DB build)
    cur.execute(
        """
        CREATE TABLE db_versions_tmp_tab AS
        SELECT r.iprscan_sig_lib_rel_id, l.library, l.version
        FROM INTERPRO.iprscan2dbcode r
        INNER JOIN iprscan.ipm_signature_library_release@ispro l
        ON (r.iprscan_sig_lib_rel_id=l.id)
        """
    )

    # Add indexes to db_versions_tmp_tab table
    cur.execute(
        """
        ALTER TABLE db_versions_tmp_tab add constraint DB_VERSIONS_TMP_TAB_PK primary key(iprscan_sig_lib_rel_id)
        """
    )

    cur.execute(
        """
        DROP TABLE :1
        """,
        [lookup_table]
    )





    # Get list of analysis IDs that will be included in the Berkley DB build
    cur.execute(
        """
        SELECT *
        FROM :1
        """,
        [lookup_table]
    )

    # analyses = []
    # for row in cur:
    #     analyses.append(row)

    analyses = cur.fetchall()

    print(analyses)


    # Get list of analysis IDs that will be included in the Berkley DB build
    analyses  = data_exchanger.get_list_from_table('db_versions_tmp_tab')
    #analyses = [[37,'PROSITE_PROFILES','20.119'], [50,'SFLD','2']] # Small example test case commented out!
    if not analyses:
        logger.error("No analyses found in the iprscan.db_versions_tmp_tab table")
        exit (1)











print('HEllo')

build_tmp_tables('iprscan/cebula@ora-dlvm-118.ebi.ac.uk:1521/IPTST')

