import re

import oracledb

from interpro7dw.utils import logger
from interpro7dw.utils.oracle import lob_as_str
from interpro7dw.utils.store import BasicStore, KVStoreBuilder, KVStore
from .databases import get_databases_codes
from .entries import load_entries, load_signatures


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
REPR_DOMAINS_DATABASES = {
    # Order: Pfam first
    "pfam": 0,
    # No priority for other databases
    "cdd": 1,
    "ncbifam": 1,
    "profile": 1,
    "smart": 1
}
REPR_DOMAINS_TYPES = {"domain", "repeat"}


def get_fragments(pos_start: int, pos_end: int, fragments: str) -> list[dict]:
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


def select_representative_domains(signatures: dict[str, dict]):
    domains = []
    for accession, match in signatures.items():
        try:
            rank = REPR_DOMAINS_DATABASES[match["database"].lower()]
        except KeyError:
            continue

        for i, loc in enumerate(match["locations"]):
            for j, frag in enumerate(loc["fragments"]):
                if not frag["representative"]:
                    continue

                domains.append({
                    "offset": (accession, i, j),
                    "start": frag["start"],
                    "end": frag["end"],
                    "length": frag["end"] - frag["start"] + 1,
                    "rank": rank
                })

    domains.sort(key=lambda x: (-x["length"], x["rank"]))
    for i, dom1 in enumerate(domains):
        start1 = dom1["start"]
        end1 = dom1["end"]

        for dom2 in domains[i+1:]:
            start2 = dom2["start"]
            end2 = dom2["end"]
            overlap = min(end1, end2) - max(start1, start2) + 1

            if overlap >= 0.7 * dom2["length"]:
                # Shorter domain significantly overlapped: discard it
                k, x, y = dom2["offset"]
                fragment = signatures[k]["locations"][x]["fragments"][y]
                fragment["representative"] = False


def export_uniprot_matches(uri: str, proteins_file: str, output: str,
                           processes: int = 1, tempdir: str | None = None):
    logger.info("starting")

    with KVStore(proteins_file) as store:
        keys = store.get_keys()

    with KVStoreBuilder(output, keys=keys, tempdir=tempdir) as store:
        con = oracledb.connect(uri)
        cur = con.cursor()
        entries = load_entries(cur)
        signatures = load_signatures(cur)

        cur.execute(
            """
            SELECT PROTEIN_AC, METHOD_AC, MODEL_AC, FEATURE, 
                   POS_FROM, POS_TO, FRAGMENTS, SCORE
            FROM INTERPRO.MATCH
            UNION ALL
            SELECT PROTEIN_AC, METHOD_AC, NULL, NULL,
                   POS_FROM, POS_TO, NULL, NULL
            FROM INTERPRO.FEATURE_MATCH PARTITION (ANTIFAM)
            """
        )
        i = 0
        for rec in cur:
            store.add(rec[0], (
                rec[1],  # signature acc
                rec[2],  # model acc or subfamily acc (PANTHER)
                rec[3],  # ancestral node ID (PANTHER)
                rec[7],  # score
                get_fragments(rec[4], rec[5], rec[6])
            ))

            i += 1
            if i % 1e8 == 0:
                logger.info(f"{i:>15,}")

        logger.info(f"{i:>15,}")
        cur.close()
        con.close()

        size = store.get_size()
        store.build(apply=merge_uniprot_matches,
                    processes=processes,
                    extraargs=[signatures, entries])

        size = max(size, store.get_size())
        logger.info(f"temporary files: {size / 1024 ** 2:.0f} MB")

    logger.info("done")


def merge_uniprot_matches(matches: list[tuple], signatures: dict,
                          entries: dict) -> tuple[dict, dict]:
    entry_matches = {}
    signature_matches = {}

    panther_subfamily = re.compile(r"PTHR\d+:SF\d+")

    for signature_acc, model_acc, feature, score, fragments in matches:
        if signature_acc in signature_matches:
            match = signature_matches[signature_acc]
        else:
            signature = signatures[signature_acc]
            match = signature_matches[signature_acc] = {
                "name": signature["name"],
                "database": signature["database"],
                "type": signature["type"],
                "evidence": signature["evidence"],
                "entry": signature["entry"],
                "locations": []
            }

            if match["entry"] and match["entry"] not in entry_matches:
                entry = entries[match["entry"]]
                entry_matches[match["entry"]] = {
                    "name": entry["name"],
                    "database": "INTERPRO",
                    "type": entry["type"],
                    "parent": entry["parent"],
                    "locations": []
                }

        for f in fragments:
            f["representative"] = (
                    match["type"].lower() in REPR_DOMAINS_TYPES and
                    match["database"].lower() in REPR_DOMAINS_DATABASES
            )

        location = {
            "fragments": fragments,
            "model": model_acc or signature_acc,
            "score": score
        }

        if model_acc and panther_subfamily.fullmatch(model_acc):
            location["subfamily"] = {
                "accession": model_acc,
                "name": signatures[model_acc]["name"],
                "node": feature
            }

        match["locations"].append(location)

        if match["entry"]:
            entry_matches[match["entry"]]["locations"].append(fragments)

    """
    Sort signature locations using their leftmost fragment
    (expects individual locations to be sorted by fragment)
    """
    for match in signature_matches.values():
        match["locations"].sort(key=lambda l: (l["fragments"][0]["start"],
                                               l["fragments"][0]["end"]))

    select_representative_domains(signature_matches)

    # Merge overlapping matches
    for match in entry_matches.values():
        condensed = []
        for start, end in condense_locations(match["locations"]):
            condensed.append({
                "fragments": [{
                    "start": start,
                    "end": end,
                    "dc-status": DC_STATUSES['S'],
                    "representative": False
                }],
                "model": None,
                "score": None
            })

        match["locations"] = condensed

    return signature_matches, entry_matches


def export_residues(uri: str, output: str):
    logger.info("starting")

    with BasicStore(output, mode="w") as store:
        for i, item in enumerate(_iter_residues(uri)):
            store.write(item)

            if (i + 1) % 1e7 == 0:
                logger.info(f"{i + 1:>15,}")

        logger.info(f"{i + 1:>15,}")

    logger.info("done")


def _iter_residues(uri: str):
    con = oracledb.connect(uri)
    cur = con.cursor()

    cur.execute("SELECT DBCODE, DBSHORT FROM INTERPRO.CV_DATABASE")
    databases = dict(cur.fetchall())

    cur.execute(
        """
        SELECT METHOD_AC, NAME 
        FROM INTERPRO.METHOD
        WHERE NAME IS NOT NULL
        """
    )
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


def export_features(uri: str, proteins_file: str, output: str,
                    processes: int = 1, tempdir: str | None = None):
    logger.info("starting")

    with KVStore(proteins_file) as store:
        keys = store.get_keys()

    with KVStoreBuilder(output, keys=keys, tempdir=tempdir) as store:
        con = oracledb.connect(uri)
        cur = con.cursor()
        cur.execute(
            """
            SELECT PROTEIN_AC, DBCODE, METHOD_AC, POS_FROM, POS_TO, SEQ_FEATURE
            FROM INTERPRO.FEATURE_MATCH
            """
        )

        i = 0
        for rec in cur:
            store.add(rec[0], rec[1:])

            i += 1
            if i % 1e8 == 0:
                logger.info(f"{i:>15,}")

        logger.info(f"{i:>15,}")

        cur.execute(
            """
            SELECT A.METHOD_AC, A.NAME, A.DESCRIPTION, B.DBCODE, B.DBSHORT, 
                   D.ABBREV 
            FROM INTERPRO.FEATURE_METHOD A
            INNER JOIN INTERPRO.CV_DATABASE B 
                ON A.DBCODE = B.DBCODE
            LEFT OUTER JOIN INTERPRO.IPRSCAN2DBCODE C
              ON B.DBCODE = C.DBCODE
            LEFT OUTER JOIN INTERPRO.CV_EVIDENCE D 
                ON C.EVIDENCE = D.CODE
            """
        )

        features = {}
        for acc, name, descr, dbcode, dbname, evidence in cur:
            try:
                db = features[dbcode]
            except KeyError:
                db = features[dbcode] = {}

            db[acc] = (name, descr, dbname, evidence)

        cur.close()
        con.close()

        size = store.get_size()
        store.build(apply=_merge_feature_matches,
                    processes=processes,
                    extraargs=[features])

        size = max(size, store.get_size())
        logger.info(f"temporary files: {size / 1024 ** 2:.0f} MB")

    logger.info("done")


def _merge_feature_matches(matches: list[tuple], features: dict) -> list[dict]:
    databases = {}
    for dbcode, accession, pos_start, pos_end, seq_feature in matches:
        try:
            db = databases[dbcode]
        except KeyError:
            db = databases[dbcode] = {}

        try:
            feature = db[accession]
        except KeyError:
            name, descr, dbname, evidence = features[dbcode][accession]
            feature = db[accession] = {
                "accession": accession,
                "name": name,
                "description": descr,
                "database": dbname,
                "evidence": evidence,
                "locations": []
            }

        if seq_feature is None and feature["database"].lower() == "mobidblt":
            seq_feature = "Consensus Disorder Prediction"

        feature["locations"].append((pos_start, pos_end, seq_feature))

    results = []
    for db in databases.values():
        for feature in db.values():
            feature["locations"].sort(key=lambda x: (x[0], x[1]))
            results.append(feature)

    # Sort features by the leftmost domain
    results.sort(key=lambda x: (x["locations"][0][0], x["locations"][0][1]))

    return results


def export_isoforms(uri: str, output: str):
    con = oracledb.connect(uri)
    cur = con.cursor()

    entries = load_entries(cur)
    signatures = load_signatures(cur)

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

        fragments = get_fragments(pos_start, pos_end, frags)
        isoform["matches"].append((sig_acc, model_acc, None, score, fragments))

    cur.close()
    con.close()

    with BasicStore(output, "w") as store:
        for isoform in isoforms.values():
            matches = isoform.pop("matches")
            isoform["matches"] = merge_uniprot_matches(matches, signatures,
                                                       entries)
            store.write(isoform)


def export_uniparc_matches(uri: str, proteins_file: str, output: str,
                           processes: int = 8, tempdir: str | None = None):
    logger.info("starting")

    with KVStore(proteins_file) as store:
        keys = store.get_keys()

    with KVStoreBuilder(output, keys=keys, tempdir=tempdir) as store:
        con = oracledb.connect(uri)
        cur = con.cursor()

        entries = load_entries(cur)
        signatures = load_signatures(cur)
        dbcodes, _ = get_databases_codes(cur)
        params = [f":{i}" for i in range(len(dbcodes))]

        cur.execute("SELECT MAX(UPI) FROM UNIPARC.PROTEIN")
        max_upi, = cur.fetchone()

        # SEQ_FEATURE -> contains the alignment for ProSite, HAMAP, FunFam
        cur.execute(
            f"""
            WITH ANALYSES AS (
                SELECT IPRSCAN_SIG_LIB_REL_ID AS ID
                FROM INTERPRO.IPRSCAN2DBCODE
                WHERE DBCODE IN ({','.join(params)})
            )
            SELECT UPI, METHOD_AC, MODEL_AC, SEQ_START, SEQ_END, SCORE, 
                   SEQ_FEATURE, FRAGMENTS
            FROM IPRSCAN.MV_IPRSCAN
            WHERE ANALYSIS_ID IN (SELECT ID FROM ANALYSES)
              AND UPI <= :{len(params)} 
            """,
            dbcodes + [max_upi]
        )

        i = 0
        for rec in cur:
            if rec[1] in signatures:
                store.add(rec[0], rec[1:])

            i += 1
            if i % 1e9 == 0:
                logger.info(f"{i:>15,}")

        logger.info(f"{i:>15,}")
        cur.close()
        con.close()

        size = store.get_size()
        store.build(apply=_merge_uniparc_matches,
                    processes=processes,
                    extraargs=[signatures, entries])

        size = max(size, store.get_size())
        logger.info(f"temporary files: {size / 1024 ** 2:.0f} MB")

    logger.info("done")


def _merge_uniparc_matches(matches: list[tuple], signatures: dict,
                           entries: dict) -> dict:
    signature_matches = {}
    for sig_acc, mod_acc, start, end, score, aln, fragments in matches:
        if sig_acc in signature_matches:
            match = signature_matches[sig_acc]
        else:
            signature = signatures[sig_acc]
            entry_acc = signature["entry"]
            if entry_acc:
                entry = {
                    "accession": entry_acc,
                    "name": entries[entry_acc]["name"],
                    "type": entries[entry_acc]["type"],
                    "parent": entries[entry_acc]["parent"],
                }
            else:
                entry = None

            match = signature_matches[sig_acc] = {
                "name": signature["short_name"],
                "database": signature["database"],
                "evidence": signature["evidence"],
                "entry": entry,
                "model": mod_acc,
                "locations": []
            }

        match["locations"].append((start, end, score, aln, fragments))

    # Sort locations
    for match in signature_matches.values():
        match["locations"].sort(key=lambda x: (x[0], x[1]))

    return signature_matches
