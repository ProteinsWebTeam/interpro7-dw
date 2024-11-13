import os
import pickle
import shutil
from multiprocessing import Process, Queue

import oracledb

from interpro7dw.utils import logger
from interpro7dw.utils.store import BasicStore
from .entries import load_entries, load_signatures


HMM_BOUNDS = {
    "[]": "COMPLETE",
    "[.": "N_TERMINAL_COMPLETE",
    ".]": "C_TERMINAL_COMPLETE",
    "..": "INCOMPLETE",
}


def export(uri: str, outdir: str, processes: int = 8):
    logger.info("starting")

    try:
        shutil.rmtree(outdir)
    except FileNotFoundError:
        pass

    os.makedirs(outdir, mode=0o775)
    inqueue = Queue()    # parent -> workers
    outqueue = Queue()   # workers -> parent
    workers = []

    # Start workers
    for _ in range(max(1, processes - 1)):
        p = Process(target=export_matches,
                    args=(uri, inqueue, outqueue))
        p.start()
        workers.append(p)

    # Export proteins and send tasks to workers
    task_count = protein_count = 0
    for proteins in iter_proteins(uri):
        filepath = os.path.join(outdir, f"{task_count:06d}")
        task_count += 1
        protein_count += len(proteins)

        with open(filepath, "wb") as fh:
            pickle.dump(proteins, fh, pickle.HIGHEST_PROTOCOL)

        inqueue.put(filepath)

    # Poison pill
    for _ in workers:
        inqueue.put(None)

    logger.info(f"{protein_count:,} proteins exported")

    milestone = step = 5
    for i in range(task_count):
        outqueue.get()
        progress = (i + 1) / task_count * 100
        if progress >= milestone:
            logger.info(f"\t{progress:.0f}%")
            milestone += step


def iter_proteins(uri: str):
    con = oracledb.connect(uri)
    cur = con.cursor()
    try:
        # TODO: remove limit after tests
        cur.execute(
            """
            SELECT UPI, LEN, CRC64, MD5
            FROM UNIPARC.PROTEIN
            WHERE UPI <= 'UPI00004C4B40'
            ORDER BY UPI
            """
        )
        while proteins := cur.fetchmany(size=1000000):
            yield proteins
    finally:
        cur.close()
        con.close()


def export_matches(uri: str, inqueue: Queue, outqueue: Queue):
    con = oracledb.connect(uri)
    cur = con.cursor()
    entries = load_entries(cur)
    signatures = load_signatures(cur, include_features=True)
    cur.close()
    con.close()

    for filepath in iter(inqueue.get, None):
        with open(filepath, "rb") as fh:
            all_proteins = pickle.load(fh)

        con = oracledb.connect(uri)
        cur = con.cursor()

        with BasicStore(filepath, "w") as bs:
            for i in range(0, len(all_proteins), 10000):
                batch_proteins = {}
                for upi, length, crc64, md5 in all_proteins[i:i + 10000]:
                    batch_proteins[upi] = {
                        "length": length,
                        "crc64": crc64,
                        "md5": md5
                    }

                start = min(batch_proteins.keys())
                stop = max(batch_proteins.keys())

                # Get matches for proteins with UPI between `start` and `stop`
                matches = get_matches(cur, start, stop, entries, signatures)

                # Get sites
                sites = get_sites(cur, start, stop)

                # Merge sites in matches
                matches = merge_matches_sites(matches, sites)

                for upi, protein in batch_proteins.items():
                    protein["matches"] = matches.pop(upi, [])

                bs.write(batch_proteins)

        cur.close()
        con.close()

        outqueue.put(None)


def get_matches(cur: oracledb.Cursor, start: str, stop: str,
                entries: dict, signatures: dict) -> dict[str, dict[str, dict]]:
    proteins = {}
    cur.execute(
        """
        SELECT UPI, METHOD_AC, MODEL_AC,
               SEQ_START, SEQ_END, HMM_START, HMM_END, HMM_LENGTH, HMM_BOUNDS,
               ENVELOPE_START, ENVELOPE_END, SEQSCORE, SEQEVALUE, SCORE, EVALUE,
               SEQ_FEATURE, FRAGMENTS
        FROM IPRSCAN.MV_IPRSCAN
        WHERE UPI BETWEEN :1 AND :2
        """,
        [start, stop]
    )
    for (upi, signature_acc, model_acc, seq_start, seq_end, hmm_start, hmm_end,
         hmm_length, hmm_bounds, env_start, env_end, seq_score, seq_evalue,
         dom_score, dom_evalue, seq_feature, fragments) in cur.fetchall():
        try:
            matches = proteins[upi]
        except KeyError:
            matches = proteins[upi] = {}

        key = model_acc or signature_acc
        try:
            match = matches[key]
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

            match = matches[key] = {
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

        match["locations"].append({
            "seq_start": seq_start,
            "seq_end": seq_end,
            "hmm_start": hmm_start,
            "hmm_end": hmm_end,
            "hmm_length": hmm_length,
            "hmm_bounds": hmm_bounds,
            "env_start": env_start,
            "env_end": env_end,
            "dom_evalue": dom_evalue,
            "dom_score": dom_score,
            "fragments": fragments,
            "feature": seq_feature
        })

    # Sort locations
    for matches in proteins.values():
        for match in matches.values():
            match["locations"].sort(
                key=lambda x: (x["seq_start"], x["seq_end"]))

    return proteins


def get_sites(cur: oracledb.Cursor,
              from_upi: str,
              to_upi: str) -> dict[str, dict]:
    proteins = {}
    cur.execute(
        """
        SELECT S.UPI, D.DBSHORT, D.DBNAME, V.VERSION, S.METHOD_AC, 
               S.LOC_START, S.LOC_END, S.RESIDUE, S.RESIDUE_START, 
               S.RESIDUE_END, S.DESCRIPTION
        FROM IPRSCAN.SITE S
        INNER JOIN INTERPRO.IPRSCAN2DBCODE I2D 
            ON S.ANALYSIS_ID = I2D.IPRSCAN_SIG_LIB_REL_ID
        INNER JOIN INTERPRO.CV_DATABASE D 
            ON I2D.DBCODE = D.DBCODE
        INNER JOIN INTERPRO.DB_VERSION V 
            ON D.DBCODE = V.DBCODE
        WHERE S.UPI BETWEEN :1 AND :2
        """,
        [from_upi, to_upi]
    )
    for (upi, dbshort, dbname, version, signature_acc, loc_start, loc_end,
         residues, res_start, res_end, description) in cur.fetchall():
        try:
            sites = proteins[upi]
        except KeyError:
            sites = proteins[upi] = {}

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
            site_locations = descriptions[description]
        except KeyError:
            site_locations = descriptions[description] = []

        site_locations.append({
            "start": res_start,
            "end": res_end,
            "residue": residues
        })

    return proteins


def merge_matches_sites(matches: dict, sites: dict) -> dict[str, list[dict]]:
    results = {}
    for upi, protein_matches in matches.items():
        seq_sites = sites.pop(upi, {})
        obj = results[upi] = []
        for match in protein_matches.values():
            sig_sites = seq_sites.pop(match["signature"]["accession"], {})
            for loc in match["locations"]:
                loc_key = (loc["seq_start"], loc["seq_end"])
                loc_sites = sig_sites.pop(loc_key, {})
                loc["sites"] = format_sites(loc_sites)

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
