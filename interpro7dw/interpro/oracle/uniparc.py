import os
import pickle
import shutil
from multiprocessing import Process, Queue

import oracledb

from interpro7dw.utils import logger
from interpro7dw.utils.store import BasicStore
from .entries import load_entries, load_signatures
from .matches import get_fragments, get_hmm_boundaries


def export(uri: str, outdir: str, processes: int = 8, chunksize: int = 100000,
           limit: int = 0):
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
    for proteins in iter_proteins(uri, chunksize=chunksize):
        filepath = os.path.join(outdir, f"{task_count:06d}.dat")
        task_count += 1
        protein_count += len(proteins)

        with open(filepath, "wb") as fh:
            pickle.dump(proteins, fh, pickle.HIGHEST_PROTOCOL)

        inqueue.put(filepath)

        if protein_count >= limit > 0:
            break

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

    for p in workers:
        p.join()

    logger.info("done")


def iter_proteins(uri: str, chunksize: int = 100000):
    con = oracledb.connect(uri)
    cur = con.cursor()
    try:
        cur.execute(
            """
            SELECT UPI, LEN, CRC64, MD5
            FROM UNIPARC.PROTEIN
            ORDER BY UPI
            """
        )
        while proteins := cur.fetchmany(size=chunksize):
            yield proteins
    finally:
        cur.close()
        con.close()


def export_matches(uri: str, inqueue: Queue, outqueue: Queue):
    con = oracledb.connect(uri)
    cur = con.cursor()
    entries = load_entries(cur)
    signatures = load_signatures(cur, include_features=True)

    cur.execute(
        """
        SELECT I2D.IPRSCAN_SIG_LIB_REL_ID, D.DBSHORT, D.DBNAME, V.VERSION, 
               EVI.ABBREV
        FROM INTERPRO.IPRSCAN2DBCODE I2D
        INNER JOIN INTERPRO.CV_DATABASE D ON I2D.DBCODE = D.DBCODE
        INNER JOIN INTERPRO.DB_VERSION V ON D.DBCODE = V.DBCODE
        INNER JOIN INTERPRO.CV_EVIDENCE EVI ON I2D.EVIDENCE = EVI.CODE        
        """
    )
    analyses = {}
    for analysis_id, dbshort, dbname, version, evidence in cur.fetchall():
        analyses[analysis_id] = {
            "abbrev": dbshort,
            "name": dbname,
            "version": version,
            "evidence": evidence
        }

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
                matches = get_matches(cur, start, stop, entries, signatures,
                                      analyses)

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


def get_matches(cur: oracledb.Cursor,
                start: str,
                stop: str,
                entries: dict,
                signatures: dict,
                analyses: dict[int, dict]) -> dict[str, dict[str, dict]]:
    proteins = {}
    cur.execute(
        """
        SELECT ANALYSIS_ID ,UPI, METHOD_AC, MODEL_AC,
               SEQ_START, SEQ_END, HMM_START, HMM_END, HMM_LENGTH, HMM_BOUNDS,
               ENVELOPE_START, ENVELOPE_END, SEQSCORE, SEQEVALUE, SCORE, EVALUE,
               SEQ_FEATURE, FRAGMENTS
        FROM IPRSCAN.MV_IPRSCAN
        WHERE UPI BETWEEN :1 AND :2
        """,
        [start, stop]
    )
    for (analysis_id, upi, signature_acc, model_acc, seq_start, seq_end,
         hmm_start, hmm_end, hmm_length, hmm_bounds, env_start, env_end,
         seq_score, seq_evalue, dom_score, dom_evalue, seq_feature, fragments
         ) in cur.fetchall():
        try:
            matches = proteins[upi]
        except KeyError:
            matches = proteins[upi] = {}

        key = (analysis_id, model_acc or signature_acc)
        try:
            match = matches[key]
        except KeyError:
            if signature_acc.startswith("PIRSR"):
                signature = {
                    "short_name": signature_acc,
                    "name": None,
                    "type": "Region",
                    "entry": None,
                }
            else:
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

            analysis = analyses[analysis_id]
            match = matches[key] = {
                "signature": {
                    "accession": signature_acc,
                    "name": signature["short_name"],
                    "description": signature["name"],
                    "type": signature["type"],
                    "signatureLibraryRelease": {
                        "library": analysis["name"],
                        "version": analysis["version"],
                    },
                    "entry": entry,
                },
                "model-ac": model_acc,
                "score": seq_score,
                "evalue": seq_evalue,
                "locations": [],
                # Fields not in the InterProScan JSON output
                "extra": {
                    "dbname": analysis["abbrev"],
                    "evidence": analysis["evidence"],
                }
            }

        """
        Notes
        =====
        
        CATH-Gene3D
          - Different seq_score/seq_evalue for two locations 
                from the same model/signature? Pick the "best" one for now
        
        CDD:
          - seq_score == dom_score and seq_evalue == dom_evalue,
                but we only need the score/evalue for locations
                
        HAMAP:
          - location.score: score of individual domains stored in SEQ_SCORE, 
                so it is overwritten when we have multiple domains 
                from the same profile
          - location.cigarAlignment: stored in seq_feature
          
        PANTHER:
          - seq_score == dom_score and seq_evalue == dom_evalue,
                TODO: investigate why
          - PANTHER subfamily ID in model_acc (family ID if no subfamily found)
          - PANTHER annotation node (AN) ID in seq_feature (can be None)

        PRINTS:
          - location.score: like HAMAP
          - location.motifNumber: stored in hmm_length 
          - location.pvalue: stored in dom_evalue
          - match.graphscan: stored in seq_feature

        PIRSR:
          - PIRSR matches are not stored anyhere and code would break when trying to retrieve
            signature data from the METHOD and FEATURE_METHOD table, so manually build a sig

        PROSITE patterns
          - location.cigarAlignment: stored in seq_feature
          
        PROSITE profiles
          - location.score: like PRINTS
          - location.cigarAlignment: stored in seq_feature
        """

        library = match["signature"]["signatureLibraryRelease"]["library"]
        if library == "CATH-Gene3D":
            if seq_score > match["score"]:
                match["score"] = seq_score
            if seq_evalue < match["evalue"]:
                match["evalue"] = seq_evalue
        elif library == "CDD":
            match["score"] = match["evalue"] = None
        elif library in ("HAMAP", "PRINTS", "PROSITE profiles"):
            dom_score = seq_score
            match["score"] = None

        match["locations"].append({
            "start": seq_start,
            "end": seq_end,
            "hmmStart": hmm_start,
            "hmmEnd": hmm_end,
            "hmmLength": hmm_length,
            "hmmBounds": get_hmm_boundaries(hmm_bounds),
            "envelopeStart": env_start,
            "envelopeEnd": env_end,
            "evalue": dom_evalue,
            "score": dom_score,
            "location-fragments": get_fragments(seq_start, seq_end, fragments),
            "sequence-feature": seq_feature,
            # Fields not in the InterProScan JSON output
            "extra": {
                "fragments": fragments,
                "hmm_bounds": hmm_bounds,
            }
        })

    # Sort locations
    for matches in proteins.values():
        for match in matches.values():
            match["locations"].sort(
                key=lambda x: (x["start"], x["end"]))

    return proteins


def get_sites(cur: oracledb.Cursor,
              from_upi: str,
              to_upi: str) -> dict[str, dict]:
    proteins = {}
    cur.execute(
        """
        SELECT UPI, METHOD_AC, LOC_START, LOC_END, RESIDUE, 
               RESIDUE_START, RESIDUE_END, DESCRIPTION
        FROM IPRSCAN.SITE
        WHERE UPI BETWEEN :1 AND :2
        """,
        [from_upi, to_upi]
    )
    for (upi, signature_acc, loc_start, loc_end, residues,
         res_start, res_end, description) in cur.fetchall():
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
                loc_key = (loc["start"], loc["end"])
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
