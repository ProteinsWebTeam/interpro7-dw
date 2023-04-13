from typing import Optional

from interpro7dw.utils import logger
from interpro7dw.utils.store import KVStore, KVStoreBuilder


def export(alphafold_file: str, proteins_file: str, output: str,
           keep_fragments: bool = False, tempdir: Optional[str] = None):
    """Export proteins with AlphaFold predictions.
    :param alphafold_file: TSV file of AlphaFold predictions.
    :param proteins_file: File to KVStore of proteins.
    :param output: Output KVStore file.
    :param keep_fragments: If False, ignore proteins where a prediction is
    split into overlapping fragments (i.e. multiple predictions but for
    different segments of the protein).
    :param tempdir: Temporary directory
    """
    logger.info("starting")

    with KVStore(proteins_file) as st:
        keys = st.get_keys()

    with KVStore(proteins_file) as protein:
        with KVStoreBuilder(output, keys=keys, tempdir=tempdir) as ash:
            with open(alphafold_file, "rt") as fh:
                for line in fh:
                    """
                    Columns:
                        - UniProt accession, e.g. A8H2R3
                        - AlphaFold DB identifier, e.g. AF-A8H2R3-F1
                        - mean pLDDT of the prediction
                        - alphafold sequence hash
                    """
                    cols = line.rstrip().split()
                    uniprot_acc = cols[0]
                    alphafold_id = cols[1]
                    score = float(cols[2])
                    crc64 = cols[3]
                    try:
                        protein_info = protein[uniprot_acc]
                    except KeyError:
                        continue
                    else:
                        ash.add(uniprot_acc, (alphafold_id, score), protein_info["crc64"] == crc64)

            if keep_fragments:
                ash.build(apply=lambda x: x)
            else:
                ash.build(apply=lambda x: x if len(x) == 1 else [])

            logger.info(f"temporary files: {ash.get_size() / 1024 ** 2:.0f} MB")

    logger.info("done")
