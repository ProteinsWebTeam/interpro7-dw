import json
import os
import shutil

import rocksdict

from interpro7dw.utils import logger
from interpro7dw.utils.store import BasicStore


def build(indir: str, outdir: str):
    logger.info("starting")

    try:
        shutil.rmtree(outdir)
    except FileNotFoundError:
        pass

    os.makedirs(outdir, mode=0o775)

    opt = rocksdict.Options(raw_mode=True)
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

    db = rocksdict.Rdict(outdir, options=opt)

    files = []
    for filename in sorted(os.listdir(indir)):
        filepath = os.path.join(indir, filename)
        files.append(filepath)

    milestone = step = 5
    for i, filepath in enumerate(files):
        with BasicStore(filepath) as bs:
            for proteins in bs:
                wb = rocksdict.WriteBatch(raw_mode=True)
                for protein in proteins.values():
                    md5 = protein["md5"]
                    matches = protein["matches"]

                    # Remove extra fields
                    for match in matches:
                        del match["extra"]
                        for loc in match["locations"]:
                            del loc["extra"]

                    wb.put(md5.encode("utf-8"),
                           json.dumps(matches).encode("utf-8"))

                db.write(wb)

        progress = (i + 1) * 100 / len(files)
        if progress >= milestone:
            logger.info(f"{progress:.0f}%")
            milestone += step

    logger.info("compacting")
    db.compact_range(None, None)
    db.close()
    logger.info("done")
