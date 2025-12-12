import datetime
import glob
import heapq
import json
import os
import shutil
from concurrent.futures import ProcessPoolExecutor, as_completed

from rocksdict import Rdict, Options, WriteBatch

from interpro7dw.utils import logger
from interpro7dw.utils.store import BasicStore


METADATA = "interpro.json"


def build(
    indir: str,
    workdir: str,
    version: str,
    date: datetime.date,
    processes: int = 8,
    max_files: int = 0,
    max_records: int = 0,
):
    logger.info("sorting by MD5")
    try:
        shutil.rmtree(workdir)
    except FileNotFoundError:
        pass

    outdir = os.path.join(workdir, version)
    tmpdir = os.path.join(workdir, "tmp")
    for dirpath in [outdir, tmpdir]:
        os.makedirs(dirpath, mode=0o775)

    files = sort_by_md5(indir, tmpdir, processes=processes, limit=max_files)

    logger.info("creating RocksDB database")
    opt = Options(raw_mode=True)
    # Increase the size of the write buffer (default: 64MB)
    opt.set_write_buffer_size(256 * 1024 * 1024)
    opt.set_level_zero_file_num_compaction_trigger(4)

    # Increase the size of in-memory write buffers (default: 2)
    opt.set_max_write_buffer_number(3)
    # Increase the base file size for level 1 (default: 64MB)
    opt.set_target_file_size_base(256 * 1024 * 1024)

    # Increase the number of background compaction/flush threads
    opt.set_max_background_jobs(4)

    # https://github.com/facebook/rocksdb/wiki/RocksDB-FAQ#basic-readwrite
    opt.prepare_for_bulk_load()

    db = Rdict(outdir, options=opt)

    iterable = []
    for filepath in files:
        bs = BasicStore(filepath, mode="r", compresslevel=0)
        iterable.append(iter(bs))

    wb = WriteBatch(raw_mode=True)
    i = 0
    all_analyses = set()
    for md5, matches, analyses in heapq.merge(*iterable, key=lambda x: x[0]):
        wb.put(md5, matches)
        all_analyses |= analyses
        i += 1

        if i == max_records:
            break
        elif i % 1e6 == 0:
            db.write(wb)
            wb = WriteBatch(raw_mode=True)

            if i % 1e8 == 0:
                logger.info(f"{i:>20,} records inserted")

    db.write(wb)
    logger.info(f"{i:>20,} records inserted")

    logger.info("compacting")
    db.compact_range(None, None)
    db.close()

    with open(os.path.join(outdir, METADATA), "wt") as fh:
        json.dump(
            {
                "release": version,
                "release_date": date.strftime("%Y-%m-%d"),
                "analyses": [
                    dict(zip(("name", "version"), values)) for values in all_analyses
                ],
            },
            fh,
        )

    shutil.rmtree(tmpdir)
    logger.info("done")


def sort_by_md5(
    indir: str, outdir: str, processes: int = 8, limit: int = 0
) -> list[str]:
    with ProcessPoolExecutor(max_workers=max(1, processes - 1)) as executor:
        fs = {}
        for src in glob.glob(os.path.join(indir, "*.dat")):
            dst = os.path.join(outdir, os.path.basename(src))
            f = executor.submit(sort_file, src, dst)
            fs[f] = dst
            if len(fs) == limit:
                break

        milestone = step = 5
        errors = 0
        for i, f in enumerate(as_completed(fs)):
            try:
                f.result()
            except Exception as exc:
                logger.error(exc)
                errors += 1

            progress = (i + 1) * 100 / len(fs)
            if progress >= milestone:
                logger.info(f"{progress:.0f}%")
                milestone += step

    if errors:
        raise RuntimeError(f"{errors} errors occurred")

    return list(fs.values())


def sort_file(src: str, dst: str):
    # Sort each chunk of `src`
    files = []
    with BasicStore(src, mode="r", compresslevel=0) as bs:
        for i, proteins in enumerate(bs):
            temppath = f"{dst}.{i}.tmp"

            with BasicStore(temppath, mode="w", compresslevel=0) as bs2:
                for p in sorted(proteins.values(), key=lambda x: x["md5"]):
                    matches = []
                    analyses = set()
                    while p["matches"]:
                        match = p["matches"].pop(0)
                        siglib = match["signature"]["signatureLibraryRelease"]

                        match siglib["library"]:
                            case "AntiFam":
                                match = format_default(match, sites=False)
                            case "CATH-FunFam" | "FunFam":
                                siglib["library"] = "CATH-FunFam"
                                match = format_default(match, sites=False)
                            case "CATH-Gene3D":
                                match = format_default(match, sites=False)
                            case "CDD":
                                match = format_cdd(match)
                            case "COILS":
                                match = format_minimal(match)
                            case "HAMAP":
                                match = format_prosite(match)
                            case "MobiDB Lite" | "MobiDB-lite":
                                siglib["library"] = "MobiDB-lite"
                                match = format_mobidblite(match)
                            case "NCBIFAM":
                                match = format_default(match, sites=False)
                            case "PANTHER":
                                match = format_panther(match)
                            case "Pfam":
                                match = format_default(match, sites=False)
                            case "Phobius":
                                match = format_minimal(match)
                            case "PIRSF":
                                match = format_default(match, sites=False)
                            case "PIRSR":
                                match = format_default(match, hmm_bounds=False)
                            case "PRINTS":
                                match = format_prints(match)
                            case "PROSITE patterns":
                                match = format_prosite(match, score=False)
                            case "PROSITE profiles":
                                match = format_prosite(match)
                            case "SFLD":
                                match = format_default(match, hmm_bounds=False)
                            case (
                                "SignalP_Euk"
                                | "SignalP_Gram_positive"
                                | "SignalP_Gram_negative"
                                | "TMHMM"
                            ):
                                match = None
                            case "SMART":
                                match = format_default(
                                    match, envelope=False, sites=False
                                )
                            case "SUPERFAMILY":
                                match = format_superfamily(match)
                            case _:
                                raise ValueError(f"Unsupported database: {siglib}")

                        if match is not None:
                            # TODO: handle InterPro-N
                            match["source"] = siglib["library"]
                            matches.append(match)
                            analyses.add((siglib["library"], siglib["version"]))

                    bs2.write((p["md5"], matches, analyses))

            files.append(temppath)

    # Merge chunks in one single sorted file
    iterable = []
    for filepath in files:
        bs = BasicStore(filepath, mode="r", compresslevel=0)
        iterable.append(iter(bs))

    with BasicStore(dst, mode="w", compresslevel=0) as bs:
        for md5, matches, analyses in heapq.merge(*iterable, key=lambda x: x[0]):
            # Add the key-value pair to add to the RocksDB
            bs.write(
                (md5.encode("utf-8"), json.dumps(matches).encode("utf-8"), analyses)
            )

    for filepath in files:
        os.unlink(filepath)


def format_default(
    match: dict, hmm_bounds: bool = True, envelope: bool = True, sites: bool = True
) -> dict:
    locations = []
    for loc in match["locations"]:
        locations.append(
            {
                "start": loc["start"],
                "end": loc["end"],
                "hmmStart": loc["hmmStart"],
                "hmmEnd": loc["hmmEnd"],
                "hmmLength": loc["hmmLength"],
                "evalue": loc["evalue"],
                "score": loc["score"],
                "location-fragments": loc["location-fragments"],
            }
        )

        if hmm_bounds:
            locations[-1]["hmmBounds"] = loc["hmmBounds"]

        if envelope:
            locations[-1]["envelopeStart"] = loc["envelopeStart"]
            locations[-1]["envelopeEnd"] = loc["envelopeEnd"]

        if sites:
            locations[-1]["sites"] = loc["sites"]

    return {
        "signature": match["signature"],
        "model-ac": match["model-ac"],
        "score": match["score"],
        "evalue": match["evalue"],
        "locations": locations,
    }


def format_cdd(match: dict) -> dict:
    locations = []
    for loc in match["locations"]:
        locations.append(
            {
                "start": loc["start"],
                "end": loc["end"],
                "evalue": loc["evalue"],
                "score": loc["score"],
                "location-fragments": loc["location-fragments"],
                "sites": loc["sites"],
            }
        )

    return {
        "signature": match["signature"],
        "model-ac": match["model-ac"],
        "locations": locations,
    }


def format_minimal(match: dict) -> dict:
    locations = []
    for loc in match["locations"]:
        locations.append(
            {
                "start": loc["start"],
                "end": loc["end"],
                "location-fragments": loc["location-fragments"],
            }
        )

    return {
        "signature": match["signature"],
        "model-ac": match["model-ac"],
        "locations": locations,
    }


def format_mobidblite(match: dict) -> dict:
    locations = []
    for loc in match["locations"]:
        locations.append(
            {
                "start": loc["start"],
                "end": loc["end"],
                "location-fragments": loc["location-fragments"],
                "sequenceFeature": loc["sequence-feature"],
            }
        )

    return {
        "signature": match["signature"],
        "model-ac": match["model-ac"],
        "locations": locations,
    }


def format_panther(match: dict) -> dict:
    locations = []
    for loc in match["locations"]:
        locations.append(
            {
                "start": loc["start"],
                "end": loc["end"],
                # TODO: add when available in DB
                # "evalue": loc["evalue"],
                # "score": loc["score"],
                "hmmStart": loc["hmmStart"],
                "hmmEnd": loc["hmmEnd"],
                "hmmLength": loc["hmmLength"],
                "hmmBounds": loc["hmmBounds"],
                "envelopeStart": loc["envelopeStart"],
                "envelopeEnd": loc["envelopeEnd"],
                "location-fragments": loc["location-fragments"],
            }
        )

    return {
        "signature": match["signature"],
        "model-ac": match["model-ac"],
        "ancestralNode": match["locations"][0]["sequence-feature"],
        "evalue": match["locations"][0]["evalue"],
        "score": match["locations"][0]["score"],
        "locations": locations,
    }


def format_prints(match: dict) -> dict:
    locations = []
    for loc in match["locations"]:
        locations.append(
            {
                "start": loc["start"],
                "end": loc["end"],
                "pvalue": loc["evalue"],
                "score": loc["score"],
                "motifNumber": loc["hmmLength"],
                "location-fragments": loc["location-fragments"],
            }
        )

    return {
        "signature": match["signature"],
        "model-ac": match["model-ac"],
        "evalue": match["evalue"],
        "graphscan": match["locations"][0]["sequence-feature"],
        "locations": locations,
    }


def format_prosite(match: dict, score: bool = True) -> dict:
    locations = []
    for loc in match["locations"]:
        locations.append(
            {
                "start": loc["start"],
                "end": loc["end"],
                "cigarAlignment": loc["sequence-feature"],
                "location-fragments": loc["location-fragments"],
            }
        )

        if score:
            locations[-1]["score"] = loc["score"]

    return {
        "signature": match["signature"],
        "model-ac": match["model-ac"],
        "locations": locations,
    }


def format_superfamily(match: dict) -> dict:
    locations = []
    for loc in match["locations"]:
        locations.append(
            {
                "start": loc["start"],
                "end": loc["end"],
                "evalue": loc["evalue"],
                "hmmLength": loc["hmmLength"],
                "location-fragments": loc["location-fragments"],
            }
        )

    return {
        "signature": match["signature"],
        "model-ac": match["model-ac"],
        "locations": locations,
    }
