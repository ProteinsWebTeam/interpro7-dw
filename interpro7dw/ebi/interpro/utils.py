# -*- coding: utf-8 -*-

from typing import List, Mapping, Sequence, Tuple, Union
SoM = Sequence[Mapping]


DC_STATUSES = {
    # Continuous single chain domain
    "S": "CONTINUOUS",
    # N terminus discontinuous
    "N": "N_TERMINAL_DISC",
    # C terminus discontinuous
    "C": "C_TERMINAL_DISC",
    # N and C terminus discontinuous
    "NC": "NC_TERMINAL_DISC"
}
MIN_OVERLAP = 0.1


def condense_locations(locations: Sequence[Sequence[dict]]) -> List[Tuple[int, int]]:
    start = end = None
    condensed = []

    # Sort locations using their leftmost fragment
    for fragments in sorted(locations, key=lambda l: repr_fragment(l[0])):
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


def repr_fragment(fragment: dict) -> Tuple[int, int]:
    return fragment["start"], fragment["end"]


def overlaps_pdb_chain(locations: SoM, segments: SoM) -> bool:
    for loc in locations:
        # We do not consider fragmented matches
        loc_start = loc["fragments"][0]["start"]
        loc_end = max([f["end"] for f in loc["fragments"]])

        for segment in segments:
            seg_start = segment["protein_start"]
            seg_end = segment["protein_end"]

            if loc_start <= seg_end and seg_start <= loc_end:
                return True

    return False


class Table(object):
    def __init__(self, con, query: str, autocommit: bool=False,
                 buffer_size: int=100000, depends_on=None):
        self.con = con
        self.cur = con.cursor()
        self.query = query
        self.autocommit = autocommit
        self.buffer_size = buffer_size
        self.depends_on = depends_on
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
        elif self.depends_on:
            self.depends_on.flush()

        self.cur.executemany(self.query, self.rows)
        self.rows = []

        if self.autocommit:
            self.con.commit()

    def close(self):
        if self.con is not None:
            self.flush()
            self.cur.close()
            self.con = None
