# -*- coding: utf-8 -*-

import hashlib
from typing import Dict, List, Optional, Tuple, Union

import cx_Oracle
import MySQLdb

from i7dw import io, logger


MIN_OVERLAP = 0.1


def extract_frag(frag: dict) -> Tuple[int, int]:
    return frag["start"], frag["end"]


def is_overlapping(x1: int, x2: int, y1: int, y2: int) -> bool:
    return x1 <= y2 and y1 <= x2


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
        self.rows = []

        if self.autocommit:
            self.con.commit()

    def close(self):
        if self.con is not None:
            self.flush()
            self.cur.close()
            self.con = None


class DomainArchitecture(object):
    def __init__(self, entries: Dict[str, Dict]):
        self.entries = entries
        self.domains = []

    @property
    def identifier(self) -> str:
        blobs = []
        for pfam_acc, interpro_acc in self.domains:
            if interpro_acc:
                blobs.append(f"{pfam_acc}:{interpro_acc}")
            else:
                blobs.append(pfam_acc)

        return '-'.join(blobs)

    @property
    def hash(self) -> str:
        return hashlib.sha1(self.identifier.encode("utf-8")).hexdigest()

    def find(self, other_interpro_acc) -> List[str]:
        accessions = set()
        for pfam_acc, interpro_acc in self.domains:
            if interpro_acc == other_interpro_acc:
                accessions.add(pfam_acc)

        return list(accessions)

    def update(self, entries: Dict[str, List[Dict]]):
        locations = []

        # Merge all Pfam locations
        for signature_acc in entries:
            signature = self.entries[signature_acc]
            if signature["database"] == "pfam":
                entry_acc = signature["integrated"]

                for loc in entries[signature_acc]:
                    # We do not consider fragmented matches
                    locations.append({
                        "pfam": signature_acc,
                        "interpro": entry_acc,
                        "start": loc["fragments"][0]["start"],
                        "end": loc["fragments"][-1]["end"]
                    })

        self.domains = [(loc["pfam"], loc["interpro"])
                        for loc in sorted(locations, key=extract_frag)]
