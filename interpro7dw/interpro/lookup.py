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


def build(indir: str, outdir: str, version: str, date: datetime.date,
          processes: int = 8, max_files: int = 0, max_records: int = 0):
    logger.info("sorting by MD5")
    tmpdir = os.path.join(os.path.dirname(outdir),
                          f"tmp{os.path.basename(outdir)}")

    for dirpath in [outdir, tmpdir]:
        try:
            shutil.rmtree(dirpath)
        except FileNotFoundError:
            pass

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
    for md5, matches in heapq.merge(*iterable, key=lambda x: x[0]):
        wb.put(md5, matches)
        i += 1

        if i == max_records:
            break
        elif i % 1e6 == 0:
            db.write(wb)
            wb = WriteBatch(raw_mode=True)

            if i % 1e7 == 0:
                logger.info(f"{i:>20,} records inserted")

    db.write(wb)
    logger.info(f"{i:>20,} records inserted")

    logger.info("compacting")
    db.compact_range(None, None)
    db.close()

    with open(os.path.join(outdir, "interpro.json"), "wt") as fh:
        json.dump({
            "resource": "InterPro",
            "service": "Matches API",
            "release": {
                "version": version,
                "date": date.strftime("%Y-%m-%d")
            }
        }, fh)

    shutil.rmtree(tmpdir)
    logger.info("done")


def sort_by_md5(indir: str, outdir: str,
                processes: int = 8, limit: int = 0) -> list[str]:
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
                    # Remove extra fields
                    for match in p["matches"]:
                        siglib = match["signature"]["signatureLibraryRelease"]
                        if siglib["library"] == "FunFam":
                            siglib["library"] = "CATH-FunFam"

                        del match["extra"]
                        for loc in match["locations"]:
                            del loc["extra"]

                    bs2.write((p["md5"], p["matches"]))

            files.append(temppath)

    # Merge chunks in one single sorted file
    iterable = []
    for filepath in files:
        bs = BasicStore(filepath, mode="r", compresslevel=0)
        iterable.append(iter(bs))

    with BasicStore(dst, mode="w", compresslevel=0) as bs:
        for md5, matches in heapq.merge(*iterable,
                                        key=lambda x: x[0]):
            # Add the key-value pair to add to the RocksDB
            bs.write((
                md5.encode("utf-8"),
                json.dumps(matches).encode("utf-8")
            ))

    for filepath in files:
        os.unlink(filepath)
