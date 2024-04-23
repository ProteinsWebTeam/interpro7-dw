import json

from MySQLdb import OperationalError
from MySQLdb.cursors import Cursor


def create_index(cur: Cursor, statement: str):
    try:
        cur.execute(statement)
    except OperationalError as exc:
        code, message = exc.args
        if code != 1061:
            # If 1061, key already exists
            raise exc


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
