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


class MatchPostProcessor:
    def __init__(self, cur):
        cur.execute(
            """
            SELECT E.ENTRY_AC, E.NAME, ET.ABBREV, EE.PARENT_AC
            FROM INTERPRO.ENTRY E
            INNER JOIN INTERPRO.CV_ENTRY_TYPE ET
              ON E.ENTRY_TYPE = ET.CODE
            LEFT OUTER JOIN INTERPRO.ENTRY2ENTRY EE
              ON E.ENTRY_AC = EE.ENTRY_AC AND EE.RELATION = 'TY'
            WHERE E.CHECKED = 'Y'
            """
        )
        self.entries = {rec[0]: rec[1:] for rec in cur}

        self.integrated = {}
        self.signatures = {}
        cur.execute(
            """
            SELECT M.METHOD_AC, M.DESCRIPTION, D.DBSHORT, ET.ABBREV, 
                   I2D.EVIDENCE, EM.ENTRY_AC
            FROM INTERPRO.METHOD M
            INNER JOIN INTERPRO.CV_DATABASE D
              ON M.DBCODE = D.DBCODE
            INNER JOIN INTERPRO.CV_ENTRY_TYPE ET
              ON M.SIG_TYPE = ET.CODE
            INNER JOIN INTERPRO.IPRSCAN2DBCODE I2D 
              ON M.DBCODE = I2D.DBCODE          
            LEFT OUTER JOIN INTERPRO.ENTRY2METHOD EM
              ON M.METHOD_AC = EM.METHOD_AC
            """
        )
        for rec in cur:
            self.signatures[rec[0]] = rec[1:5]
            if rec[5]:
                self.integrated[rec[0]] = rec[5]

    def digest(self, matches: list[tuple]) -> tuple[dict, dict]:
        entries = {}
        signatures = {}
        for signature_acc, model_acc, score, fragments in matches:
            fragments.sort(key=lambda x: (x["start"], x["end"]))

            try:
                s = signatures[signature_acc]
            except KeyError:
                sig = self.signatures[signature_acc]
                s = signatures[signature_acc] = {
                    "accession": signature_acc,
                    "name": sig[0],
                    "database": sig[1],
                    "type": sig[2],
                    "evidence": sig[3],
                    "entry": sig[4],
                    "locations": []
                }

            s["locations"].append({
                "fragments": fragments,
                "model": model_acc or signature_acc,
                "score": score
            })

            if s["entry"]:
                e_acc = s["entry"]

                try:
                    e = entries[e_acc]
                except KeyError:
                    e_name, e_type, e_parent = self.entries[e_acc]
                    e = entries[e_acc] = {
                        "accession": e_acc,
                        "name": e_name,
                        "database": "InterPro",
                        "type": e_type,
                        "parent": e_parent,
                        "locations": []
                    }

                e["locations"].append(fragments)

        # Sort signature locations using their leftmost fragment
        for sig in signatures.values():
            sig["locations"].sort(key=lambda l: (l["fragments"][0]["start"],
                                                 l["fragments"][0]["end"]))

        # Merge overlapping matches
        for entry in entries.values():
            condensed = []
            for start, end in self.condense_locations(entry["locations"]):
                condensed.append({
                    "fragments": [{
                        "start": start,
                        "end": end,
                        "dc-status": DC_STATUSES['S']
                    }],
                    "model": None,
                    "score": None
                })

            entry["locations"] = condensed

        return signatures, entries

    @staticmethod
    def condense_locations(locations: list[list[dict]],
                           min_overlap: int = 0.1) -> list[tuple[int, int]]:
        start = end = None
        condensed = []

        """
        Sort locations using their leftmost fragment
        (assume that fragments are sorted in individual locations)
        """
        for fragments in sorted(locations,
                                key=lambda l: (l[0]["start"], l[0]["end"])):
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

                if overlap >= shortest * min_overlap:
                    # Merge
                    end = e
                    continue

            condensed.append((start, end))
            start, end = s, e

        # Adding last location
        condensed.append((start, end))
        return condensed
