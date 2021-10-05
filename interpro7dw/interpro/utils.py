from typing import Sequence


def overlaps_pdb_chain(locations: Sequence[dict],
                       segments: Sequence[dict]) -> bool:
    """Evaluate in protein matches and chain segments overlap.

    :param locations:
    :param segments:
    :return:
    """
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
