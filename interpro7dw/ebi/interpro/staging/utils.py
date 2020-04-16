# -*- coding: utf-8 -*-

import heapq
import json
import os
from typing import Sequence

from interpro7dw.utils import DataDump, deepupdate


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


def merge_xrefs(files: Sequence[str]):
    iterables = [DataDump(path) for path in files]
    _key = None
    _xrefs = None

    for key, xrefs in heapq.merge(*iterables, key=lambda x: x[0]):
        if key != _key:
            if _key is not None:
                yield _key, _xrefs

            _key = key
            _xrefs = xrefs

        deepupdate(xrefs, _xrefs, replace=False)

    if _key is not None:
        yield _key, _xrefs

    for datadump, path in zip(iterables, files):
        datadump.close()
        os.remove(path)
