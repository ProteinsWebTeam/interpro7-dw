from typing import Dict, List, Tuple


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
MIN_OVERLAP = 0.1  # NOT percent


def repr_frag(f: dict) -> Tuple[int, int]:
    return f["start"], f["end"]


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
