# -*- coding: utf-8 -*-

import json

import MySQLdb

from interpro7dw.utils import url2dict


def reduce(src: dict) -> dict:
    dst = {}
    for key, value in src.items():
        if isinstance(value, dict):
            dst[key] = reduce(value)
        elif isinstance(value, (list, set, tuple)):
            dst[key] = len(value)
        else:
            dst[key] = value

    return dst


def jsonify(obj, nullable: bool = True):
    if isinstance(obj, (float, int)):  # todo: remove before merging
        raise TypeError(f"{type(obj)}")

    if obj:
        if isinstance(obj, set):
            return json.dumps(list(obj))
        else:
            return json.dumps(obj)
    elif nullable:
        return None
    else:
        return json.dumps(None)


def drop_database(url: str):
    con = MySQLdb.connect(**url2dict(url))
    cur = con.cursor()

    try:
        cur.execute("DROP DATABASE interpro")
    except MySQLdb.OperationalError as exc:
        code = exc.args[0]
        if code == 1008:
            # Can't drop database '<name>'; database doesn't exist
            pass
        else:
            raise exc
    finally:
        cur.close()
        con.close()
