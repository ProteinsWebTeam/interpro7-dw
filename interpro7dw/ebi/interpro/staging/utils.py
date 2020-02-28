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
    if obj or not nullable:
        return json.dumps(obj)
    else:
        return None
