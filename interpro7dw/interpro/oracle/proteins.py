import gzip
import os
from typing import Optional

import cx_Oracle

from interpro7dw.utils import logger
from interpro7dw.utils.oracle import lob_as_str
from interpro7dw.utils.store import BasicStore, KVStoreBuilder, KVStore


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
        self.entries = {}
        for rec in cur:
            self.entries[rec[0]] = {
                "name": rec[1],
                "type": rec[2],
                "parent": rec[3]
            }

        cur.execute(
            """
            SELECT M.METHOD_AC, M.NAME, M.DESCRIPTION, D.DBSHORT, ET.ABBREV, 
                   EVI.ABBREV, EM.ENTRY_AC
            FROM INTERPRO.METHOD M
            INNER JOIN INTERPRO.CV_DATABASE D
              ON M.DBCODE = D.DBCODE
            INNER JOIN INTERPRO.CV_ENTRY_TYPE ET
              ON M.SIG_TYPE = ET.CODE
            INNER JOIN INTERPRO.IPRSCAN2DBCODE I2D 
              ON M.DBCODE = I2D.DBCODE
            INNER JOIN INTERPRO.CV_EVIDENCE EVI
              ON I2D.EVIDENCE = EVI.CODE
            LEFT OUTER JOIN (
                SELECT E.ENTRY_AC, EM.METHOD_AC
                FROM INTERPRO.ENTRY E
                INNER JOIN INTERPRO.ENTRY2METHOD EM
                  ON E.ENTRY_AC = EM.ENTRY_AC
                WHERE E.CHECKED = 'Y'
            ) EM ON M.METHOD_AC = EM.METHOD_AC
            """
        )
        self.signatures = {}
        for rec in cur:
            self.signatures[rec[0]] = {
                "short_name": rec[1],
                "name": rec[2],
                "database": rec[3],
                "type": rec[4],
                "evidence": rec[5],
                "entry": rec[6]
            }

    def digest(self, matches: list[tuple]) -> tuple[dict, dict]:
        entries = {}
        signatures = {}
        for signature_acc, model_acc, score, fragments in matches:
            try:
                s = signatures[signature_acc]
            except KeyError:
                s = signatures[signature_acc] = {
                    "name": self.signatures[signature_acc]["name"],
                    "database": self.signatures[signature_acc]["database"],
                    # "type": self.signatures[signature_acc]["type"],
                    "evidence": self.signatures[signature_acc]["evidence"],
                    "entry": self.signatures[signature_acc]["entry"],
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
                    e = entries[e_acc] = {
                        "name": self.entries[e_acc]["name"],
                        "database": "InterPro",
                        "type": self.entries[e_acc]["type"],
                        "parent": self.entries[e_acc]["parent"],
                        "locations": []
                    }

                e["locations"].append(fragments)

        """
        Sort signature locations using their leftmost fragment
        (expects individual locations to be sorted by fragment)
        """
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

    def digest_uniparc(self, matches: list[tuple]) -> dict:
        signatures = {}
        for sig_acc, mod_acc, start, end, score, aln, frags in matches:
            try:
                s = signatures[sig_acc]
            except KeyError:
                if self.signatures[sig_acc]["entry"]:
                    entry_acc = self.signatures[sig_acc]["entry"]
                    entry = {
                        "accession": entry_acc,
                        "name": self.entries[entry_acc]["name"],
                        "type": self.entries[entry_acc]["type"],
                        "parent": self.entries[entry_acc]["parent"]
                    }
                else:
                    entry = None

                s = signatures[sig_acc] = {
                    "name": self.signatures[sig_acc]["short_name"],
                    "database": self.signatures[sig_acc]["database"],
                    "evidence": self.signatures[sig_acc]["evidence"],
                    "entry": entry,
                    "model": mod_acc,
                    "locations": []
                }

            s["locations"].append((start, end, score, aln, frags))

        for s in signatures.values():
            s["locations"].sort()

        return signatures

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


def _iter_features(uri: str):
    con = cx_Oracle.connect(uri)
    cur = con.cursor()
    cur.execute(
        """
        SELECT M.METHOD_AC, M.NAME, D.DBSHORT, EVI.ABBREV
        FROM INTERPRO.FEATURE_METHOD M
        INNER JOIN INTERPRO.CV_DATABASE D
          ON M.DBCODE = D.DBCODE
        INNER JOIN INTERPRO.IPRSCAN2DBCODE I2D
          ON D.DBCODE = I2D.DBCODE
        INNER JOIN INTERPRO.CV_EVIDENCE EVI
          ON I2D.EVIDENCE = EVI.CODE
        """
    )
    features_info = {row[0]: row[1:] for row in cur}

    cur.execute(
        """
        SELECT PROTEIN_AC, METHOD_AC, POS_FROM, POS_TO, SEQ_FEATURE
        FROM INTERPRO.FEATURE_MATCH
        ORDER BY PROTEIN_AC
        """
    )

    protein_acc = None
    features = {}
    for prot_acc, feat_acc, pos_start, pos_end, seq_feature in cur:
        if prot_acc != protein_acc:
            if protein_acc:
                yield protein_acc, _sort_features(features)

            protein_acc = prot_acc
            features = {}

        try:
            feature = features[feat_acc]
        except KeyError:
            name, database, evidence = features_info[feat_acc]
            feature = features[feat_acc] = {
                "name": name,
                "database": database,
                "evidence": evidence,
                "locations": [(pos_start, pos_end, seq_feature)]
            }

        if seq_feature is None and feature["database"].lower() == "mobidblt":
            seq_feature = "Consensus Disorder Prediction"

        feature["locations"].append((pos_start, pos_end, seq_feature))

    cur.close()
    con.close()

    if protein_acc:
        yield protein_acc, _sort_features(features)


def _sort_features(features: dict) -> dict:
    for feature in features.values():
        feature["locations"].sort()

    return features


def export_features(uri: str, output: str):
    logger.info("starting")

    with BasicStore(output, mode="w") as store:
        for i, item in enumerate(_iter_features(uri)):
            store.write(item)

            if (i + 1) % 1e7 == 0:
                logger.info(f"{i + 1:>15,}")

        store.close()
        logger.info(f"{i + 1:>15,}")

    logger.info("done")


def _get_fragments(pos_start: int, pos_end: int, fragments: str) -> list[dict]:
    if fragments:
        result = []
        for frag in fragments.split(','):
            # Format: START-END-STATUS
            s, e, t = frag.split('-')
            result.append({
                "start": int(s),
                "end": int(e),
                "dc-status": DC_STATUSES[t]
            })

        result.sort(key=lambda x: (x["start"], x["end"]))
    else:
        result = [{
            "start": pos_start,
            "end": pos_end,
            "dc-status": DC_STATUSES['S']  # Continuous
        }]

    return result


def export_isoforms(uri: str, output: str):
    con = cx_Oracle.connect(uri)
    cur = con.cursor()

    processor = MatchPostProcessor(cur)

    cur.outputtypehandler = lob_as_str
    cur.execute(
        """
        SELECT V.PROTEIN_AC, V.VARIANT, V.LENGTH, V.CRC64, 
               P.SEQ_SHORT, P.SEQ_LONG
        FROM INTERPRO.VARSPLIC_MASTER V
        INNER JOIN UNIPARC.PROTEIN P ON V.CRC64 = P.CRC64
        """
    )

    isoforms = {}
    for rec in cur:
        variant_acc = rec[0] + '-' + str(rec[1])
        isoforms[variant_acc] = {
            "accession": variant_acc,
            "protein": rec[0],
            "length": rec[2],
            "crc64": rec[3],
            "sequence": rec[4] or rec[5],
            "matches": []
        }

    # PROTEIN_AC is actually PROTEIN-VARIANT (e.g. Q13733-1)
    cur.execute(
        """
        SELECT PROTEIN_AC, METHOD_AC, MODEL_AC, SCORE, POS_FROM, POS_TO, 
               FRAGMENTS
        FROM INTERPRO.VARSPLIC_MATCH V
        """
    )

    for var_acc, sig_acc, model_acc, score, pos_start, pos_end, frags in cur:
        try:
            isoform = isoforms[var_acc]
        except KeyError:
            continue

        fragments = _get_fragments(pos_start, pos_end, frags)
        isoform["matches"].append((sig_acc, model_acc, score, fragments))

    cur.close()
    con.close()

    with BasicStore(output, "w") as store:
        for isoform in isoforms.values():
            matches = isoform.pop("matches")
            signatures, entries = processor.digest(matches)
            isoform["matches"] = signatures, entries
            store.write(isoform)


def _iter_matches(uri: str):
    con = cx_Oracle.connect(uri)
    cur = con.cursor()
    processor = MatchPostProcessor(cur)

    cur.execute(
        """
        SELECT PROTEIN_AC, METHOD_AC, MODEL_AC, POS_FROM, POS_TO, FRAGMENTS, 
               SCORE
        FROM INTERPRO.MATCH
        ORDER BY PROTEIN_AC
        """
    )

    protein_acc = None
    matches = []
    for prot_acc, sig_acc, model_acc, pos_start, pos_end, frags, score in cur:
        if prot_acc != protein_acc:
            if protein_acc:
                yield protein_acc, *processor.digest(matches)

            protein_acc = prot_acc
            matches = []

        fragments = _get_fragments(pos_start, pos_end, frags)
        matches.append((sig_acc, model_acc, score, fragments))

    cur.close()
    con.close()

    if protein_acc:
        yield protein_acc, *processor.digest(matches)


def export_matches(uri: str, output: str):
    logger.info("starting")

    with KVStoreBuilder(output, keys=[], cachesize=10000) as store:
        for i, (acc, signatures, entries) in enumerate(_iter_matches(uri)):
            store.append(acc, (signatures, entries))

            if (i + 1) % 1e7 == 0:
                logger.info(f"{i + 1:>15,}")

        store.close()
        logger.info(f"{i + 1:>15,}")

    logger.info("done")


def export_proteins(uri: str, output: str):
    logger.info("starting")
    with KVStoreBuilder(output, keys=[], cachesize=10000) as store:
        con = cx_Oracle.connect(uri)
        cur = con.cursor()
        cur.execute(
            """
            SELECT PROTEIN_AC, NAME, DBCODE, CRC64, LEN, 
                   TO_CHAR(TIMESTAMP, 'YYYY-MM-DD'), FRAGMENT, TO_CHAR(TAX_ID) 
            FROM INTERPRO.PROTEIN
            ORDER BY PROTEIN_AC
            """
        )

        for i, rec in enumerate(cur):
            store.append(rec[0], {
                "identifier": rec[1],
                "reviewed": rec[2] == 'S',
                "crc64": rec[3],
                "length": rec[4],
                "date": rec[5],
                "fragment": rec[6] == 'Y',
                "taxid": rec[7]
            })

            if (i + 1) % 10000000 == 0:
                logger.info(f"{i + 1:>15,}")

        cur.close()
        con.close()
        store.close()
        logger.info(f"{i + 1:>15,}")

    logger.info("done")


def _iter_residues(uri: str):
    con = cx_Oracle.connect(uri)
    cur = con.cursor()

    cur.execute("SELECT DBCODE, DBSHORT FROM INTERPRO.CV_DATABASE")
    databases = dict(cur.fetchall())

    cur.execute("SELECT METHOD_AC, NAME FROM INTERPRO.METHOD")
    signatures = dict(cur.fetchall())

    cur.execute(
        """
        SELECT PROTEIN_AC, METHOD_AC, DBCODE, DESCRIPTION, RESIDUE, 
               RESIDUE_START, RESIDUE_END
        FROM INTERPRO.SITE_MATCH
        ORDER BY PROTEIN_AC
        """
    )

    protein_acc = None
    matches = {}
    for prot_acc, sig_acc, dbcode, descr, res, pos_start, pos_end in cur:
        if prot_acc != protein_acc:
            if protein_acc:
                yield protein_acc, _sort_residues(matches)

            protein_acc = prot_acc
            matches = {}

        try:
            signature = matches[sig_acc]
        except KeyError:
            signature = matches[sig_acc] = {
                "name": signatures.get(sig_acc),
                "database": databases[dbcode],
                "descriptions": {}
            }

        try:
            signature["descriptions"][descr].append((res, pos_start, pos_end))
        except KeyError:
            signature["descriptions"][descr] = [(res, pos_start, pos_end)]

    cur.close()
    con.close()

    if protein_acc:
        yield protein_acc, _sort_residues(matches)


def _sort_residues(matches: dict) -> dict:
    for signature in matches.values():
        for locations in signature["descriptions"].values():
            locations.sort(key=lambda x: (x[1], x[2]))

    return matches


def export_residues(uri: str, output: str):
    logger.info("starting")

    with BasicStore(output, mode="w") as store:
        for i, item in enumerate(_iter_residues(uri)):
            store.write(item)

            if (i + 1) % 1e7 == 0:
                logger.info(f"{i + 1:>15,}")

        logger.info(f"{i + 1:>15,}")

    logger.info("done")


def export_sequences(uri: str, kvstore: str, output: str,
                     tempdir: Optional[str] = None):
    logger.info("starting")
    with KVStore(kvstore) as s:
        keys = s.get_keys()

    with KVStoreBuilder(output, keys=keys, tempdir=tempdir) as store:
        con = cx_Oracle.connect(uri)
        cur = con.cursor()
        cur.outputtypehandler = lob_as_str
        cur.execute(
            """
            SELECT UX.AC, UP.SEQ_SHORT, UP.SEQ_LONG
            FROM UNIPARC.XREF UX
            INNER JOIN UNIPARC.PROTEIN UP ON UX.UPI = UP.UPI
            WHERE UX.DBID IN (2, 3)  -- Swiss-Prot and TrEMBL
            AND UX.DELETED = 'N'
            """
        )

        for i, (accession, seq_short, seq_long) in enumerate(cur):
            sequence = seq_short or seq_long
            store.add(accession, gzip.compress(sequence.encode("utf-8")))

            if (i + 1) % 1e7 == 0:
                logger.info(f"{i + 1:>15,}")

        logger.info(f"{i + 1:>15,}")
        cur.close()
        con.close()

        store.build(apply=store.get_first)
        logger.info(f"temporary files: {store.get_size() / 1024 ** 2:.0f} MB")

    logger.info("done")


def _iter_uniparc_matches(cur: cx_Oracle.Cursor):
    processor = MatchPostProcessor(cur)

    # SEQ_FEATURE -> contains the alignment for ProSite and HAMAP
    cur.execute(
        """
        SELECT UPI, METHOD_AC, MODEL_AC, SEQ_START, SEQ_END, SCORE, 
               SEQ_FEATURE, FRAGMENTS
        FROM IPRSCAN.MV_IPRSCAN
        ORDER BY UPI
        """
    )

    upi = None
    matches = []
    for _upi, sig_acc, mod_acc, start, end, score, alignment, frags in cur:
        if _upi != upi:
            if upi:
                yield upi, processor.digest_uniparc(matches)

            upi = _upi
            matches = []

        fragments = _get_fragments(start, end, frags)
        matches.append((sig_acc, mod_acc, start, end, score, alignment,
                        fragments))

    if upi:
        yield upi, processor.digest_uniparc(matches)


def export_uniparc(uri: str, output: str):
    con = cx_Oracle.connect(uri)
    cur = con.cursor()

    logger.info("exporting UniParc proteins")
    proteins_tmp = f"{output}.tmp"
    with KVStoreBuilder(proteins_tmp, keys=[], cachesize=10000) as store:
        cur.execute(
            """
            SELECT UPI, LEN, CRC64
            FROM UNIPARC.PROTEIN
            ORDER BY UPI
            """
        )

        for i, (upi, length, crc64) in enumerate(cur):
            store.append(upi, (length, crc64))

            if (i + 1) % 1e8 == 0:
                logger.info(f"{i + 1:>15,}")

        store.close()
        logger.info(f"{i + 1:>15,}")

    logger.info("exporting UniParc matches")
    with KVStore(proteins_tmp) as st1, BasicStore(output, mode="w") as st2:
        for i, (upi, signatures) in enumerate(_iter_uniparc_matches(cur)):
            if (i + 1) % 1e7 == 0:
                logger.info(f"{i + 1:>15,}")

            try:
                length, crc64 = st1[upi]
            except KeyError:
                """
                This may append if matches are calculated against sequences
                in UAPRO instead of UAREAD
                """
                continue

            st2.write((upi, length, crc64, signatures))

        logger.info(f"{i + 1:>15,}")

    os.unlink(proteins_tmp)
    logger.info("done")
