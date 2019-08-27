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


def repr_frag(f: dict) -> Tuple[int, int]:
    return f["start"], f["end"]


def is_overlapping(x1: int, x2: int, y1: int, y2: int) -> bool:
    return x1 <= y2 and y1 <= x2


def condense(to_condense: Dict[str, List]) -> Dict[str, List]:
    condensed = {}
    for entry_ac, locations in to_condense.items():
        start = end = None
        _locations = []

        # Sort location by the position of the leftmost fragment
        for frags in sorted(locations, key=lambda l: repr_frag(l[0])):
            """
            We do not consider fragmented matches:
                - `s` is the leftmost start position
                - `e` is the rightmost end position
                (assuming `frags` is sorted by (start, end) keys)
            """
            s = frags[0]["start"]  # leftmost start position
            e = frags[-1]["end"]  # rightmost end position

            if start is None:
                start = s
                end = e
            elif s > end:
                """
                      end
                   [----] [----]
                          s
                -> new location
                """
                _locations.append((start, end))
                start = s
                end = e
            elif e > end:
                """
                        end
                   [----]
                     [------]
                            e
                -> extend
                """
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
                "seq_feature": None,
                "model_acc": None
            })

    return condensed
