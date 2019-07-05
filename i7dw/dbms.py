import re
from typing import Union

import cx_Oracle
import MySQLdb
import MySQLdb.cursors


def connect(uri: str, sscursor: bool=False, encoding: str='utf-8') -> tuple:
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
    def __init__(self, con: Union[cx_Oracle.Connection, MySQLdb.Connection],
                 query: str, autocommit: bool=False, buffer_size: int=100000):
        self.con = con
        self.cur = con.cursor()
        self.query = query
        self.autocommit = autocommit
        self.buffer_size = buffer_size
        self.rows = []
        self.count = 0

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __del__(self):
        self.close()

    def _execute(self, record: Union[dict, tuple]):
        self.rows.append(record)
        self.count += 1

        if len(self.rows) == self.buffer_size:
            self.flush()

    def insert(self, record: Union[dict, tuple]):
        self._execute(record)

    def update(self, record: Union[dict, tuple]):
        self._execute(record)

    def delete(self, record: Union[dict, tuple]):
        self._execute(record)

    def flush(self):
        if not self.rows:
            return

        self.cur.executemany(self.query, self.rows)
        self.count += len(self.rows)
        self.rows = []

        if self.autocommit:
            self.con.commit()

    def close(self):
        if self.cur is not None:
            self.flush()
            self.cur.close()
            self.con = None
