import gzip
import os
import pickle
from multiprocessing import Process, Queue

import oracledb

from interpro7dw.utils import logger
from interpro7dw.utils.store import Directory, KVStoreBuilder, KVStore
from .entries import load_entries, load_signatures
from .matches import get_fragments


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

    inqueue = Queue()
    outqueue = Queue()
    workers = []
    for _ in range(max(processes - 1, 1)):
        p = Process(target=export_matches_in_range,
                    args=(uri, inqueue, outqueue))
        p.start()
        workers.append(p)

    directory = Directory(tempdir=tempdir)

    files = []
    ready = []
    with KVStore(proteins_file) as store:
        keys = store.get_keys()

        for i, start in enumerate(keys):
            try:
                stop = keys[i + 1]
                include_stop = False
            except IndexError:
                stop = store.max()
                include_stop = True
            finally:
                filepath = directory.mktemp(createfile=False)
                inqueue.put((start, stop, include_stop, filepath))
                files.append(filepath)
                ready.append(False)

    for _ in workers:
        inqueue.put(None)

    with KVStoreBuilder(output, keys=[], cachesize=10000) as store:
        i = 0
        milestone = step = 5
        for _ in range(len(files)):
            filepath = outqueue.get()

            # Flag returned file as read
            j = files.index(filepath)
            ready[j] = True

            # Load files that are ready in order
            for j in range(i, len(files)):
                if not ready[j]:
                    break

                filepath = files[j]
                with gzip.open(filepath, "rb") as fh:
                    matches = pickle.load(fh)

                for upi in sorted(matches):
                    store.append(upi, matches[upi])

                i = j
                os.unlink(filepath)

            progress = (i + 1) * 100 / len(files)
            if progress >= milestone:
                logger.info(f"{progress:.0f}%")
                while milestone <= progress:
                    milestone += step

    for p in workers:
        p.join()

    directory.remove()
    logger.info("done")


def export_matches_in_range(
    uri: str,
    inqueue: Queue,
    outqueue: Queue
):
    con = oracledb.connect(uri)
    cur = con.cursor()
    signatures = load_signatures(cur, include_features=True)
    entries = load_entries(cur)

    for from_upi, to_upi, include_stop, filepath in iter(inqueue.get, None):
        matches = get_matches(cur, from_upi, to_upi, include_stop,
                              signatures, entries)
        sites = get_sites(cur, from_upi, to_upi, include_stop)
        matches = merge_matches_sites(matches, sites)

        with gzip.open(filepath, "wb", compresslevel=6) as fh:
            pickle.dump(matches, fh, pickle.HIGHEST_PROTOCOL)

        outqueue.put(filepath)

    cur.close()
    con.close()


def get_matches(
    cur: oracledb.Cursor,
    from_upi: str,
    to_upi: str,
    include_stop: bool,
    signatures: dict[str, dict],
    entries: dict[str, dict]
) -> dict[str, dict[str, dict]]:
    if include_stop:
        where_upi = "M.UPI >= :1 AND M.UPI <= :2"
    else:
        where_upi = "M.UPI >= :1 AND M.UPI < :2"

    cur.execute(
        f"""
        SELECT M.UPI, D.DBNAME, V.VERSION, M.METHOD_AC, M.MODEL_AC,
               M.SEQSCORE, M.SEQEVALUE, M.SEQ_START, M.SEQ_END,
               M.SCORE, M.EVALUE, M.HMM_START, M.HMM_END, M.HMM_LENGTH, 
               M.HMM_BOUNDS, M.ENVELOPE_START, M.ENVELOPE_END,
               M.SEQ_FEATURE, M.FRAGMENTS
        FROM IPRSCAN.MV_IPRSCAN M
        INNER JOIN INTERPRO.IPRSCAN2DBCODE I2D 
            ON M.ANALYSIS_ID = I2D.IPRSCAN_SIG_LIB_REL_ID
        INNER JOIN INTERPRO.CV_DATABASE D 
            ON I2D.DBCODE = D.DBCODE
        INNER JOIN INTERPRO.DB_VERSION V 
            ON D.DBCODE = V.DBCODE
        WHERE {where_upi}
        """,
        [from_upi, to_upi]
    )

    matches = {}
    while rows := cur.fetchmany(size=100000):
        for (upi, db_name, db_version, signature_acc, model_acc, seq_score,
             seq_evalue, loc_start, loc_end, dom_score, dom_evalue, hmm_start,
             hmm_end, hmm_length, hmm_bound, env_start, env_end, seq_feature,
             fragments) in rows:
            try:
                seq_matches = matches[upi]
            except KeyError:
                seq_matches = matches[upi] = {}

            match_key = model_acc or signature_acc
            try:
                match = seq_matches[match_key]
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

                match = seq_matches[match_key] = {
                    "signature": {
                        "accession": signature_acc,
                        "name": signature["short_name"],
                        "description": signature["name"],
                        "signatureLibraryRelease": {
                            "library": db_name,
                            "version": db_version,
                        },
                        "entry": entry,
                    },
                    "evidence": signature["evidence"],
                    "model-ac": model_acc,
                    "score": seq_score,
                    "evalue": seq_evalue,
                    "locations": [],
                }

            match["locations"].append({
                "start": loc_start,
                "end": loc_end,
                "hmmStart": hmm_start,
                "hmmEnd": hmm_end,
                "hmmLength": hmm_length,
                "hmmBounds": HMM_BOUNDS.get(hmm_bound),
                "evalue": dom_evalue,
                "score": dom_score,
                "envelopeStart": env_start,
                "envelopeEnd": env_end,
                "location-fragments": get_fragments(loc_start, loc_end,
                                                    fragments),
                "sequence-feature": seq_feature,
            })

    # Sort locations
    for seq_matches in matches.values():
        for match in seq_matches.values():
            match["locations"].sort(key=lambda x: (x["start"], x["end"]))

    return matches


def get_sites(
    cur: oracledb.Cursor,
    from_upi: str,
    to_upi: str,
    include_stop: bool,
) -> dict[
        str,                                    # UniParc ID
        dict[
            str,                                # Model accession
            dict[
                tuple[int, int],                # Location (start, end)
                dict[
                    str,                        # Description
                    list[dict]                  # Residues
                ]
            ]
        ]
]:
    if include_stop:
        where_upi = "UPI >= :1 AND UPI <= :2"
    else:
        where_upi = "UPI >= :1 AND UPI < :2"

    cur.execute(
        f"""
        SELECT UPI, METHOD_AC, LOC_START, LOC_END, RESIDUE, RESIDUE_START, 
               RESIDUE_END, DESCRIPTION
        FROM IPRSCAN.SITE
        WHERE {where_upi}
        """,
        [from_upi, to_upi]
    )

    sites = {}
    while rows := cur.fetchmany(size=100000):
        for (upi, signature_acc, loc_start, loc_end, residues, res_start,
             res_end, descr) in rows:
            try:
                locations = sites[signature_acc]
            except KeyError:
                locations = sites[signature_acc] = {}

            loc_key = (loc_start, loc_end)
            try:
                descriptions = locations[loc_key]
            except KeyError:
                descriptions = locations[loc_key] = {}

            try:
                site_locations = descriptions[descr]
            except KeyError:
                site_locations = descriptions[descr] = []

            site_locations.append({
                "start": res_start,
                "end": res_end,
                "residue": residues
            })

    return sites


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
