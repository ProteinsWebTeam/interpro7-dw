# -*- coding: utf-8 -*-

import json


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


def jsonify(obj, nullable=True):
    if obj is None:
        return None if nullable else json.dumps(None)
    elif isinstance(obj, set):
        return json.dumps(list(obj))
    else:
        return json.dumps(obj)
