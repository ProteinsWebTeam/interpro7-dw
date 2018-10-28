#!/usr/bin/env python
# -*- coding: utf-8 -*-

class Supermatch(object):
    def __init__(self, entry_ac, entry_root, start, end):
        self.entries = {(entry_ac, entry_root)}
        self.start = start
        self.end = end

    def in_same_hierarchy(self, other):
        for acc_1, root_1 in self.entries:
            for acc_2, root_2 in other.entries:
                if root_1 != root_2:
                    return False

        return True

    def merge_if_overlap(self, other, min_overlap):
        overlap = min(self.end, other.end) - max(self.start, other.start) + 1
        shortest = min(self.end - self.start, other.end - other.start) + 1

        if (overlap / shortest * 100) >= min_overlap:
            self.entries = self.entries | other.entries
            self.start = min(self.start, other.start)
            self.end = max(self.end, other.end)
            return self
        else:
            return None

    def format(self):
        return '{}:{:.0f}-{:.0f}'.format(
            self.format_entries(),
            self.start,
            self.end
        )

    def format_entries(self):
        return '&'.join(sorted([re.sub(r'IPR0*', '', accession) for accession in self.get_entries()]))

    def get_entries(self):
        return [entry_ac for entry_ac, entry_root in self.entries]

    def __eq__(self, other):
        return (
                isinstance(other, Supermatch) and
                self.start == other.start and
                self.end == other.end and
                self.entries == other.entries
        )


class SupermatchSet(object):
    def __init__(self, supermatch):
        self.supermatches = [supermatch]

    def add(self, candidate, min_overlap):
        if not self.supermatches[0].in_same_hierarchy(candidate):
            return False

        merged = None
        for sm in self.supermatches:
            merged = sm.merge_if_overlap(candidate, min_overlap)

            if merged is not None:
                break

        if merged is None:
            self.supermatches.append(candidate)
        else:
            # Merged supermatch: we now need to remove overlaps between the newly merged supermatch and others
            indexes_ok = set()

            while True:
                index = None

                for i, sm in enumerate(self.supermatches):
                    if sm == merged or i in indexes_ok:
                        continue

                    if merged.merge_if_overlap(sm, min_overlap):
                        # Overlap so merged, we now have to remove the merged supermatch (sm)
                        index = i
                        break
                    else:
                        # No overlap, might be skipped during next iteration
                        indexes_ok.add(i)

                if index is None:
                    # No move overlaps
                    break
                else:
                    self.supermatches.pop(index)

        return True


def merge_supermatches(supermatches, min_overlap=20):
    sets = []

    for sm in supermatches:
        in_set = False

        for s in sets:
            in_set = s.add(sm, min_overlap)

            if in_set:
                break

        if not in_set:
            sets.append(SupermatchSet(sm))

    return sets
