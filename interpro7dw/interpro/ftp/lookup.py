import os
import tarfile

from interpro7dw.utils import logger

_LOOKUP_TAR = "matches-api-data.tar.gz"


def archive(indir: str, version: str, outdir: str):
    logger.info("starting")
    os.makedirs(outdir, exist_ok=True)
    prefix = f"interpro-{version}"
    lookup_dir = os.path.join(indir, version)
    with tarfile.open(os.path.join(outdir, _LOOKUP_TAR), "w:gz") as fh:
        for filename in os.listdir(lookup_dir):
            filepath = os.path.join(lookup_dir, filename)
            fh.add(filepath, arcname=os.path.join(prefix, filename))

    logger.info("done")
