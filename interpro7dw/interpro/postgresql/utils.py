import json
import re

import psycopg


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


def connect(uri: str) -> psycopg.Connection:
    return psycopg.connect(**parse_uri(uri))


def parse_uri(uri: str) -> dict:
    m = re.match(r'([^/]+)/([^@]+)@([^:]+):(\d+)/(\w+)', uri)

    if m is None:
        raise RuntimeError(f"invalid connection string: {uri}")

    return {
        "user": m.group(1),
        "password": m.group(2),
        "host": m.group(3),
        "port": int(m.group(4)),
        "dbname": m.group(5),
    }
