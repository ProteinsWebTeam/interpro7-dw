from cx_Oracle import DB_TYPE_BLOB, DB_TYPE_LONG
from cx_Oracle import DB_TYPE_CLOB, DB_TYPE_LONG_RAW


def lob_as_str(cursor, name, default_type, size, precision, scale):
    if default_type == DB_TYPE_BLOB:
        return cursor.var(DB_TYPE_LONG_RAW, arraysize=cursor.arraysize)
    if default_type == DB_TYPE_CLOB:
        return cursor.var(DB_TYPE_LONG, arraysize=cursor.arraysize)
