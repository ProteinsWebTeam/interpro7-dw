import glob
import heapq
import json
import os
import shutil
from concurrent.futures import ProcessPoolExecutor, as_completed
from tempfile import mkstemp

from rocksdict import Rdict, Options, SstFileWriter

from interpro7dw.utils import logger
from interpro7dw.utils.store import BasicStore


def build(indir: str, outdir: str, processes: int = 8,
          tempdir: str | None = None):
    logger.info("starting")

    try:
        shutil.rmtree(outdir)
    except FileNotFoundError:
        pass

    os.makedirs(outdir, mode=0o775)

    if tempdir:
        os.makedirs(tempdir, exist_ok=True)

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

    errors = 0
    with ProcessPoolExecutor(max_workers=max(1, processes - 1)) as executor:
        fs = []
        for filepath in glob.glob(os.path.join(indir, "*.dat")):
            f = executor.submit(create_sst, filepath, tempdir)
            fs.append(f)

        milestone = step = 5
        for i, f in enumerate(as_completed(fs)):
            try:
                path = f.result()
            except Exception as exc:
                logger.error(exc)
                errors += 1
            else:
                db.ingest_external_file([path])
                os.unlink(path)

            progress = (i + 1) * 100 / len(fs)
            if progress >= milestone:
                logger.info(f"{progress:.0f}%")
                milestone += step

    if errors:
        db.close()
        raise RuntimeError(f"{errors} occurred")

    logger.info("compacting")
    db.compact_range(None, None)
    db.close()
    logger.info("done")


def create_sst(filepath: str, tempdir: str | None = None) -> str:
    stores = []
    with BasicStore(filepath, mode="r", compresslevel=0) as bs:
        for i, proteins in enumerate(bs):
            fd, temppath = mkstemp(dir=tempdir)
            os.close(fd)

            bs2 = BasicStore(temppath, mode="w", compresslevel=0)
            for p in sorted(proteins.values(), key=lambda x: x["md5"]):
                # Remove extra fields
                for match in p["matches"]:
                    del match["extra"]
                    for loc in match["locations"]:
                        del loc["extra"]

                bs2.write((p["md5"], p["matches"]))

            bs2.close()
            stores.append(bs2)

    sstfile = f"{filepath}.sst"
    writer = SstFileWriter(Options(raw_mode=True))
    writer.open(sstfile)
    iterable = [iter(bs) for bs in stores]
    for md5, matches in heapq.merge(*iterable, key=lambda x: x[0]):
        key = md5.encode("utf-8")
        value = json.dumps(matches).encode("utf-8")
        writer[key] = value

    writer.finish()

    for bs in stores:
        os.unlink(bs.file)

    return sstfile
