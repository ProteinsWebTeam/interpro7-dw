#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re

import cx_Oracle
import MySQLdb
import MySQLdb.cursors


def connect(uri, sscursor=False, encoding='utf-8'):
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
        con = cx_Oracle.connect(user, passwd, dsn,
                                encoding=encoding, nencoding=encoding)
        return con, con.cursor()
    elif driver == 'mysql':
        # supports 'utf8', not 'utf-8'
        encoding = encoding.replace('-', '').lower()

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


class Populator(object):
    def __init__(self, uri: str, query: str, autocommit: bool=False,
                 chunk_size: int=100000):
        self.con, self.cur  = connect(uri)
        self.query = query
        self.autocommit = autocommit
        self.chunk_size = chunk_size
        self.rows = []
        self.count = 0

    def __del__(self):
        self.close()

    def insert(self, row:tuple):
        self.rows.append(row)

        if len(self.rows) == self.chunk_size:
            self.flush(self.autocommit)

    def flush(self, commit: bool=False):
        self.cur.executemany(self.query, self.rows)

        if commit:
            self.con.commit()

        self.count += len(self.rows)
        self.rows = []

    def close(self):
        if self.con is not None:
            if self.rows:
                self.flush(True)
            self.cur.close()
            self.con.close()
            self.con = None
