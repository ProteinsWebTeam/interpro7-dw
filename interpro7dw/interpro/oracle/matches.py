import re

import oracledb

from interpro7dw.utils import logger
from interpro7dw.utils.oracle import lob_as_str
from interpro7dw.utils.store import BasicStore, KVStoreBuilder, KVStore
from .databases import get_databases_codes
from .entries import (load_entries, load_signatures,
                      REPR_DOM_DATABASES, REPR_DOM_TYPES,
                      REPR_FAM_DATABASES, REPR_FAM_TYPES)


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
MAX_DOM_BY_GROUP = 20
DOM_OVERLAP_THRESHOLD = 0.3


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


def select_repr_domains(domains: list[dict]):
    # Sort by boundaries
    domains.sort(key=lambda d: (d["fragments"][0]["start"],
                                d["fragments"][-1]["end"]))

    # Group overlapping domains together
    domain = domains[0]
    domain["residues"] = calc_coverage(domain)
    stop = domain["fragments"][-1]["end"]
    group = [domain]
    groups = []

    for domain in domains[1:]:
        domain["residues"] = calc_coverage(domain)
        start = domain["fragments"][0]["start"]

        if start <= stop:
            group.append(domain)
            stop = max(stop, domain["fragments"][-1]["end"])
        else:
            groups.append(group)
            group = [domain]
            stop = domain["fragments"][-1]["end"]

    groups.append(group)

    # Select representative domain in each group
    for group in groups:
        """
        Only consider the "best" N domains of the group, 
        otherwise the number of possible combinations/sets is too high 
        (if M domains, max number of combinations is `2 ^ M`)
        """
        group = sorted(group,
                       key=lambda d: (-len(d["residues"]), d["rank"])
                       )[:MAX_DOM_BY_GROUP]

        nodes = set(range(len(group)))
        graph = {i: nodes - {i} for i in nodes}

        for i, dom_a in enumerate(group):
            for j in range(i + 1, len(group)):
                dom_b = group[j]
                if eval_overlap(dom_a, dom_b, DOM_OVERLAP_THRESHOLD):
                    graph[i].remove(j)
                    graph[j].remove(i)

        # Find possible domains combinations
        subgroups = resolve_domains(graph)

        # Find the best combination
        max_coverage = 0
        max_pfams = 0
        best_subgroup = None
        for subgroup in subgroups:
            coverage = set()
            pfams = 0
            _subgroup = []

            for i in subgroup:
                domain = group[i]
                coverage |= domain["residues"]
                if domain["rank"] == 0:
                    pfams += 1

                _subgroup.append(domain)

            coverage = len(coverage)
            if coverage < max_coverage:
                continue
            elif coverage > max_coverage or pfams > max_pfams:
                max_coverage = coverage
                max_pfams = pfams
                best_subgroup = _subgroup

        # Flag selected representative domains
        for domain in best_subgroup:
            domain["representative"] = True


def calc_coverage(domain: dict) -> set[int]:
    residues = set()
    for f in domain["fragments"]:
        residues |= set(range(f["start"], f["end"] + 1))

    return residues


def eval_overlap(dom_a: dict, dom_b: dict, threshold: float) -> bool:
    overlap = dom_a["residues"] & dom_b["residues"]
    if overlap:
        len_a = len(dom_a["residues"])
        len_b = len(dom_b["residues"])
        return len(overlap) / min(len_a, len_b) >= threshold

    return False


def resolve_domains(graph: dict[int, set[int]]) -> list[set[int]]:
    def is_valid(candidate: list[int]) -> bool:
        for node_a in candidate:
            for node_b in candidate:
                if node_a != node_b and node_a not in graph[node_b]:
                    return False

        return True

    def make_sets(current_set: list[int], remaining_nodes: list[int]):
        if is_valid(current_set):
            if not remaining_nodes:
                all_sets.append(set(current_set))
                return True
        else:
            return False

        current_node = remaining_nodes[0]
        remaining_nodes = remaining_nodes[1:]

        # Explore two possibilities at each step of the recursion
        # 1) current node is added to the set under consideration
        make_sets(current_set + [current_node], remaining_nodes)
        # 2) current node is not added to the set
        make_sets(current_set, remaining_nodes)

    all_sets = []
    make_sets([], list(graph.keys()))
    return all_sets


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
            SELECT PROTEIN_AC, METHOD_AC, MODEL_AC, 
                   POS_FROM, POS_TO, FRAGMENTS, SCORE
            FROM INTERPRO.MATCH
            UNION ALL
            SELECT PROTEIN_AC, METHOD_AC, NULL,
                   POS_FROM, POS_TO, NULL, NULL
            FROM INTERPRO.FEATURE_MATCH PARTITION (ANTIFAM)
            """
        )
        i = 0
        for rec in cur:
            store.add(rec[0], (
                rec[1],  # signature acc
                rec[2],  # model acc or subfamily acc (PANTHER)
                rec[6],  # score
                get_fragments(rec[3], rec[4], rec[5])
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
    domains = []
    families = []
    regions = []
    for signature_acc, model_acc, score, fragments in matches:
        signature = signatures[signature_acc]

        database = signature["database"].lower()
        sig_type = signature["type"].lower()
        match = {
            "signature": signature_acc,
            "model": model_acc or signature_acc,
            "score": score,
            "fragments": fragments,
        }

        if database in REPR_DOM_DATABASES and sig_type in REPR_DOM_TYPES:
            match["rank"] = REPR_DOM_DATABASES.index(database)
            domains.append(match)
        elif database in REPR_FAM_DATABASES and sig_type in REPR_FAM_TYPES:
            match["rank"] = REPR_FAM_DATABASES.index(database)
            families.append(match)
        else:
            regions.append(match)

    if domains:
        select_repr_domains(domains)

    if families:
        select_repr_domains(families)

    entry_matches = {}
    signature_matches = {}
    panther_subfamily = re.compile(r"PTHR\d+:SF\d+")
    for domain in domains + families + regions:
        signature_acc = domain["signature"]
        if signature_acc in signature_matches:
            match = signature_matches[signature_acc]
        else:
            signature = signatures[signature_acc]
            match = signature_matches[signature_acc] = {
                "name": signature["name"],
                "short_name": signature["short_name"],
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
                    "short_name": entry["short_name"],
                    "database": "INTERPRO",
                    "type": entry["type"],
                    "parent": entry["parent"],
                    "locations": []
                }

        location = {
            "fragments": domain["fragments"],
            "representative": domain.get("representative", False),
            "model": domain["model"],
            "score": domain["score"]
        }

        if panther_subfamily.fullmatch(domain["model"]):
            location["subfamily"] = {
                "accession": domain["model"],
                "name": signatures[domain["model"]]["name"],
            }

        match["locations"].append(location)

        if match["entry"]:
            entry_match = entry_matches[match["entry"]]
            entry_match["locations"].append(domain["fragments"])

    # Sort signature locations using the leftmost fragment
    for match in signature_matches.values():
        match["locations"].sort(key=lambda l: (l["fragments"][0]["start"],
                                               l["fragments"][0]["end"]))

    # Merge overlapping matches
    for match in entry_matches.values():
        condensed = []
        for start, end in condense_locations(match["locations"]):
            condensed.append({
                "fragments": [{
                    "start": start,
                    "end": end,
                    "dc-status": DC_STATUSES['S'],
                }],
                "representative": False,
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

            if dbname.lower() == "pfam-n":
                # Pfam-N not in IPRSCAN2DBCODE
                evidence = "Maskformer"

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
        isoform["matches"].append((sig_acc, model_acc, score, fragments))

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
