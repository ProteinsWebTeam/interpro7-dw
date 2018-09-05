#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re

import cx_Oracle
import MySQLdb
import MySQLdb.cursors


def connect(uri, sscursor=False, encoding='utf-8'):
    # driver:user/password@host:port/dbname[:encoding]
    # m = re.match(r'(\w+):([^/]+)/([^@]+)@([^:]+):(\d+)/(\w+)(?::([\w-]+))?', uri)
    m = re.match(r'(\w+):([^/]+)/([^@]+)@([^:]+):(\d+)/(\w+)', uri)

    if m is None:
        raise RuntimeError('invalid connection string: {}'.format(uri))

    driver = m.group(1).lower()
    user = m.group(2)
    passwd = m.group(3)
    host = m.group(4)
    port = int(m.group(5))
    db = m.group(6)

    if driver == 'oracle':
        dsn = cx_Oracle.makedsn(host, port, db)
        con = cx_Oracle.connect(user, passwd, dsn, encoding=encoding, nencoding=encoding)
        return con, con.cursor()
    elif driver == 'mysql':
        encoding = encoding.replace('-', '').lower()  # do not support hyphen

        con = MySQLdb.connect(**{
            'user': user,
            'passwd': passwd,
            'host': host,
            'port': port,
            'db': db,
            'use_unicode': encoding == 'utf8',
            'charset': encoding
        })

        if sscursor:
            return con, MySQLdb.cursors.SSCursor(con)
        else:
            return con, con.cursor()
    else:
        raise RuntimeError('driver not supported: {}'.format(driver))
