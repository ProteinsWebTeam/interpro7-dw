import oracledb

from interpro7dw.utils import logger
from interpro7dw.utils.store import KVStoreBuilder, KVStore
from .entries import load_entries, load_signatures


HMM_BOUNDS = {
    "[]": "COMPLETE",
    "[.": "N_TERMINAL_COMPLETE",
    ".]": "C_TERMINAL_COMPLETE",
    "..": "INCOMPLETE",
}


def export_proteins(uri: str, output: str):
    logger.info("starting")
    with KVStoreBuilder(output, keys=[], cachesize=10000) as store:
        con = oracledb.connect(uri)
        cur = con.cursor()
        cur.execute(
            """
            SELECT UPI, LEN, CRC64, MD5
            FROM UNIPARC.PROTEIN
            ORDER BY UPI
            """
        )

        for i, (upi, length, crc64, md5) in enumerate(cur):
            store.append(upi, (length, crc64, md5))

            if (i + 1) % 1e8 == 0:
                logger.info(f"{i + 1:>15,}")

        cur.close()
        con.close()

        store.close()
        logger.info(f"{i + 1:>15,}")

    logger.info("done")


def export_matches(uri: str, proteins_file: str, output: str,
                   processes: int = 8, tempdir: str | None = None):
    logger.info("starting")

    with KVStore(proteins_file) as store:
        keys = store.get_keys()

    with KVStoreBuilder(output, keys=keys, tempdir=tempdir) as store:
        con = oracledb.connect(uri)
        cur = con.cursor()

        cur.execute(
            f"""
            SELECT UPI, METHOD_AC, MODEL_AC,
                   SEQSCORE, SEQEVALUE, SEQ_START, SEQ_END, SCORE, EVALUE, 
                   HMM_START, HMM_END, HMM_LENGTH, HMM_BOUNDS, 
                   ENVELOPE_START, ENVELOPE_END, SEQ_FEATURE, FRAGMENTS
            FROM IPRSCAN.MV_IPRSCAN
            """,
        )

        i = 0
        for row in cur:
            store.add(row[0], row[1:])

            i += 1
            if i % 1e9 == 0:
                logger.info(f"{i:>15,}")

        logger.info(f"{i:>15,}")
        signatures = load_signatures(cur, include_features=True)
        entries = load_entries(cur)
        cur.close()
        con.close()

        size = store.get_size()
        store.build(apply=_merge_matches,
                    processes=processes,
                    extraargs=[signatures, entries])

        size = max(size, store.get_size())
        logger.info(f"temporary files: {size / 1024 ** 2:.0f} MB")

    logger.info("done")


def _merge_matches(matches: list[tuple],
                   signatures: dict,
                   entries: dict) -> dict:
    results = {}
    for (signature_acc, model_acc, seq_score,
         seq_evalue, loc_start, loc_end, dom_score, dom_evalue, hmm_start,
         hmm_end, hmm_length, hmm_bound, env_start, env_end, seq_feature,
         fragments) in matches:
        match_key = model_acc or signature_acc
        try:
            match = results[match_key]
        except KeyError:
            signature = signatures[signature_acc]
            entry_acc = signature["entry"]
            if entry_acc:
                entry = {
                    "accession": entry_acc,
                    "name": entries[entry_acc]["short_name"],
                    "description": entries[entry_acc]["name"],
                    "type": entries[entry_acc]["type"],
                    "parent": entries[entry_acc]["parent"],
                }
            else:
                entry = None

            match = results[match_key] = {
                "signature": {
                    "accession": signature_acc,
                    "name": signature["name"],
                    "description": signature["description"],
                    "database": signature["database"],
                    "evidence": signature["evidence"],
                    "entry": entry,
                },
                "model": model_acc,
                "score": seq_score,
                "evalue": seq_evalue,
                "locations": [],
            }

        match["locations"].append((
            loc_start,
            loc_end,
            hmm_start,
            hmm_end,
            hmm_length,
            hmm_bound,
            dom_evalue,
            dom_score,
            env_start,
            env_end,
            fragments,
            seq_feature,
        ))

    # Sort locations
    for match in results.values():
        match["locations"].sort(key=lambda x: (x[0], x[1]))

    return results


def export_sites(uri: str, proteins_file: str, output: str,
                 processes: int = 8, tempdir: str | None = None):
    logger.info("starting")

    with KVStore(proteins_file) as store:
        keys = store.get_keys()

    with KVStoreBuilder(output, keys=keys, tempdir=tempdir) as store:
        con = oracledb.connect(uri)
        cur = con.cursor()

        cur.execute(
            f"""
            SELECT UPI, METHOD_AC, LOC_START, LOC_END, RESIDUE, RESIDUE_START, 
                   RESIDUE_END, DESCRIPTION
            FROM IPRSCAN.SITE
                """,
        )

        i = 0
        for row in cur:
            store.add(row[0], row[1:])

            i += 1
            if i % 1e9 == 0:
                logger.info(f"{i:>15,}")

        logger.info(f"{i:>15,}")
        cur.close()
        con.close()

        size = store.get_size()
        store.build(apply=_merge_sites, processes=processes)

        size = max(size, store.get_size())
        logger.info(f"temporary files: {size / 1024 ** 2:.0f} MB")

    logger.info("done")


def _merge_sites(sites: list[tuple]) -> dict:
    results = {}
    for (signature_acc, loc_start, loc_end, residues, res_start,
         res_end, descr) in sites:
        try:
            locations = results[signature_acc]
        except KeyError:
            locations = results[signature_acc] = {}

        loc_key = (loc_start, loc_end)
        try:
            descriptions = locations[loc_key]
        except KeyError:
            descriptions = locations[loc_key] = {}

        try:
            site_locations = descriptions[descr]
        except KeyError:
            site_locations = descriptions[descr] = []

        site_locations.append((res_start, res_end, residues))

    return results


def merge_matches_sites(matches: dict, sites: dict) -> dict[str, list[dict]]:
    results = {}
    for upi, protein_matches in matches.items():
        obj = results[upi] = []
        for match in protein_matches.values():
            for location in match["locations"]:
                location["sites"] = format_sites(
                    sites
                    .get(upi, {})
                    .get(match["signature"]["accession"], {})
                    .get((location["start"], location["end"]), {})
                )

            obj.append(match)

    return results


def format_sites(sites: dict) -> list[dict]:
    results = []
    for descr, site_locations in sites.items():
        results.append({
            "description": descr,
            "numLocations": len(site_locations),
            "siteLocations": site_locations
        })

    return results
