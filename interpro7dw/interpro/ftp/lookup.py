import json
import os
import tarfile

from interpro7dw.interpro.lookup import METADATA


_LOOKUP_TAR = "matches-api-data.tar.gz"


def archive(indir: str, outdir: str):
    filepath = os.path.join(indir, METADATA)
    with open(filepath, "rt") as fh:
        version = json.load(fh)["release"]

    prefix = f"interpro-{version}"
    with tarfile.open(os.path.join(outdir, _LOOKUP_TAR), "w:gz") as fh:
        for filename in os.listdir(indir):
            filepath = os.path.join(indir, filename)
            fh.add(filepath, arcname=os.path.join(prefix, filename))
