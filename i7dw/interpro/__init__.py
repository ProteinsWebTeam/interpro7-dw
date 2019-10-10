# -*- coding: utf-8 -*-

import hashlib
from typing import Dict, List, Tuple, Union

import cx_Oracle
import MySQLdb


MIN_OVERLAP = 0.1


def extract_frag(frag: dict) -> Tuple[int, int]:
    return frag["start"], frag["end"]


def is_overlapping(x1: int, x2: int, y1: int, y2: int) -> bool:
    return x1 <= y2 and y1 <= x2


def condense_locations(locations: List[List[Dict]]) -> List[Tuple[int, int]]:
    start = end = None
    condensed = []

    # Sort locations using their leftmost fragment
    for fragments in sorted(locations, key=lambda l: extract_frag(l[0])):
        """
        1) We do not consider fragmented matches
        2) Fragments are sorted by (start, end):
            * `start` of the first frag is guaranteed to be the leftmost one
            * `end` of the last frag is NOT guaranteed to be the rightmost one
                (e.g. [(5, 100), (6, 80)])
        """
        s = fragments[0]["start"]
        e = max([f["end"] for f in fragments])

        if start is None:
            # First location
            start, end = s, e
            continue
        elif e <= end:
            # Current location within "merged" one: nothing to do
            continue
        elif s <= end:
            # Locations are overlapping (at least one residue)
            overlap = min(end, e) - max(start, s) + 1
            shortest = min(end - start, e - s) + 1

            if overlap >= shortest * MIN_OVERLAP:
                # Merge
                end = e
                continue

        condensed.append((start, end))
        start, end = s, e

    # Adding last location
    condensed.append((start, end))

    return condensed


class Table(object):
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
        self.entries = {}
        self.domains = []
        for acc, e in entries.items():
            if e["database"] == "pfam":
                self.entries[acc] = e["integrated"]

    @property
    def accession(self) -> str:
        blobs = []
        for pfam_acc, interpro_acc in self.domains:
            if interpro_acc:
                blobs.append(f"{pfam_acc}:{interpro_acc}")
            else:
                blobs.append(pfam_acc)

        return '-'.join(blobs)

    @property
    def identifier(self) -> str:
        return hashlib.sha1(self.accession.encode("utf-8")).hexdigest()

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
            try:
                entry_acc = self.entries[signature_acc]
            except KeyError:
                # Not a Pfam signature
                continue

            for loc in entries[signature_acc]:
                # We do not consider fragmented matches
                locations.append({
                    "pfam": signature_acc,
                    "interpro": entry_acc,
                    "start": loc["fragments"][0]["start"],
                    "end": max([f["end"] for f in loc["fragments"]])
                })

        self.domains = [(loc["pfam"], loc["interpro"])
                        for loc in sorted(locations, key=extract_frag)]
