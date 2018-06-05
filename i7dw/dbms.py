#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re

import cx_Oracle
import MySQLdb
import MySQLdb.cursors


def connect(uri, sscursor=False):
    # driver:user/password@host:port/dbname
    m = re.match(r'(\w+):([^/]+)/([^@]+)@([^:]+):(\d+)/(.+)', uri)

    if m is None:
        raise RuntimeError('invalid connection string: {}'.format(uri))

    driver = m.group(1).lower()
    if driver == 'oracle':
        dsn = cx_Oracle.makedsn(m.group(4), m.group(5), m.group(6))
        con = cx_Oracle.connect(m.group(2), m.group(3), dsn, encoding='utf-8', nencoding='utf-8')
        return con, con.cursor()
    elif driver == 'mysql':
        con = MySQLdb.connect(**{
            'user': m.group(2),
            'passwd': m.group(3),
            'host': m.group(4),
            'port': int(m.group(5)),
            'db': m.group(6),
            'use_unicode': True,
            'charset': 'utf8'
        })

        if sscursor:
            return con, MySQLdb.cursors.SSCursor(con)
        else:
            return con, con.cursor()
    else:
        raise RuntimeError('driver not supported: {}'.format(driver))
