import json

import MySQLdb

from interpro7dw.utils.mysql import url2dict


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
    if obj or not nullable:
        return json.dumps(obj)
    else:
        return None


def drop_database(uri: str):
    con = MySQLdb.connect(**url2dict(uri))
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
