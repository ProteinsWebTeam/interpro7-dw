# -*- coding: utf-8 -*-

import hashlib
from typing import Dict, List, Optional, Tuple, Union

import cx_Oracle
import MySQLdb

from i7dw import io, logger
# from .mysql.entry import get_entries  # TODO fix import


MIN_OVERLAP = 0.1


def extract_frag(frag: dict) -> tuple:
    return frag["start"], frag["end"]


def is_overlapping(x1: int, x2: int, y1: int, y2: int) -> bool:
    return x1 <= y2 and y1 <= x2


def condense(to_condense: Dict[str, List]) -> Dict[str, List]:
    condensed = {}
    for entry_ac, locations in to_condense.items():
        start = end = None
        _locations = []

        # We assume `locations` and `fragments` are sorted
        for fragments in locations:
            # We do not consider fragmented matches:
            s = fragments[0]["start"]
            e = fragments[-1]["end"]

            if start is None:
                start = s
                end = e
                continue
            elif s <= end:
                # matches are overlapping (at least one residue)
                overlap = min(end, e) - max(start, s) + 1
                shortest = min(end - start, e - s) + 1

                if overlap >= shortest * MIN_OVERLAP:
                    # Merge
                    end = e
                    continue

            _locations.append((start, end))
            start = s
            end = e

        _locations.append((start, end))

        condensed[entry_ac] = []
        for start, end in _locations:
            condensed[entry_ac].append({
                "fragments": [{
                    "start": start,
                    "end": end,
                    "dc-status": "CONTINUOUS"
                }],
                "model_acc": None
            })

    return condensed


def export_ida(url: str, src_matches: str, dst_ida: str,
               tmpdir: Optional[str]=None, processes: int=1,
               sync_frequency: int=1000000):

    logger.info("starting")
    pfam_entries = {}
    for e in get_entries(url).values():
        if e["database"] == "pfam":
            pfam_ac = e["accession"]
            interpro_ac = e["integrated"]
            pfam_entries[pfam_ac] = interpro_ac

    with io.Store(src_matches) as src, io.Store(dst_ida, src.keys, tmpdir) as dst:
        i = 0

        for acc, matches in src:
            dom_arch = []
            for m in matches:
                method_ac = m["method_ac"]

                if method_ac in pfam_entries:

                    interpro_ac = pfam_entries[method_ac]
                    if interpro_ac:
                        dom_arch.append("{}:{}".format(method_ac, interpro_ac))
                    else:
                        dom_arch.append("{}".format(method_ac))

            if dom_arch:
                ida = '-'.join(dom_arch)
                ida_id = hashlib.sha1(ida.encode("utf-8")).hexdigest()
                dst[acc] = (ida, ida_id)

            i += 1
            if sync_frequency and not i % sync_frequency:
                dst.sync()

            if not i % 10000000:
                logger.info("{:>12,}".format(i))

        logger.info("{:>12,}".format(i))
        dst.merge(processes=processes)
        logger.info("temporary files: {:.0f} MB".format(dst.size/1024/1024))


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
        if self.con is not None:
            self.flush()
            self.cur.close()
            self.con = None
