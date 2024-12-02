import gzip
import json
import tarfile

from interpro7dw.utils import logger
from interpro7dw.utils.store import KVStoreBuilder


def index(bfdbpath: str, output: str):
    skipped = 0
    with KVStoreBuilder(output) as kb, tarfile.open(bfdbpath, "r") as tar:
        for name in sorted(tar.getnames()):
            fh = tar.extractfile(name)
            data = json.loads(gzip.decompress(fh.read()))
            uniprot_acc = data["uniprot_entry"]["ac"]
            sequence_length = data["uniprot_entry"]["sequence_length"]
            structures = data.get("structures", [])
            if len(structures) == 1:
                summary = structures[0]["summary"]
                model_id = summary["model_identifier"]
                model_url = summary["model_url"]

                if summary["confidence_type"] == "pLDDT":
                    plddt = summary["confidence_avg_local_score"]
                    kb.append(uniprot_acc,
                              (model_id, model_url, sequence_length, plddt))
                else:
                    skipped += 1
            else:
                skipped += 1

    if skipped:
        logger.warning(f"{skipped} structures ignored")
