import copy


def copy_dict(src: dict, dst: dict, concat_or_incr: bool = False):
    for key, value in src.items():
        if key in dst:
            if isinstance(value, dict):
                copy_dict(value, dst[key], concat_or_incr)
            elif isinstance(value, (list, tuple)):
                dst[key] += value
            elif isinstance(value, set):
                dst[key] |= value
            elif isinstance(value, (int, float, str)) and concat_or_incr:
                dst[key] += value
            else:
                dst[key] = value
        else:
            dst[key] = copy.deepcopy(value)


def overlaps_pdb_chain(locations: list[dict], segments: list[dict]) -> dict[tuple, int]:
    """Evaluate in protein matches and chain segments overlap.

    :param locations:
    :param segments:
    :return:
    """

    matches_count = {}
    for loc in locations:
        # We do not consider fragmented matches
        loc_start = loc["fragments"][0]["start"]
        loc_end = max([f["end"] for f in loc["fragments"]])

        for segment in segments:
            seg_start = segment["protein_start"]
            seg_end = segment["protein_end"]

            if loc_start <= seg_end and seg_start <= loc_end:
                match =
                matches_count.update()

    return matches_count
