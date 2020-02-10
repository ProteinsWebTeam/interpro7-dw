# -*- coding: utf-8 -*-

import hashlib
from typing import List, Mapping, Optional, Sequence, Tuple


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


class DomainArchitecture(object):
    def __init__(self, pfam2interpro: Mapping[str, Optional[str]]):
        self.pfam2interpro = pfam2interpro
        self.domains = []

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

    def update(self, entries: Mapping[str, Sequence[dict]]):
        # Merge all Pfam locations
        pfam_locations = []
        for accession in entries:
            try:
                interpro_acc = self.pfam2interpro[accession]
            except KeyError:
                continue  # Not a Pfam signature

            for loc in entries[accession]:
                # We do not consider fragmented matches
                pfam_locations.append({
                    "pfam": accession,
                    "interpro": interpro_acc,
                    "start": loc["fragments"][0]["start"],
                    "end": max(f["end"] for f in loc["fragments"])
                })

        self.domains = []
        for loc in sorted(pfam_locations, key=repr_fragment):
            self.domains.append((loc["pfam"], loc["interpro"]))
