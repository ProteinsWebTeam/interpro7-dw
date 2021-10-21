def get_proteins(file: str, keep_fragments: bool = True) -> set:
    """Returns a set of protein accessions with AlphaFold predictions.

    :param file: CSV file of AlphaFold predictions.
    :param keep_fragments: If False, ignore proteins where a prediction is
    split into overlapping fragments (i.e. multiple predictions but for
    different segments of the protein).
    :return:
    """
    alphafold = {}
    with open(file, "rt") as fh:
        for line in fh:
            protein_acc, start, end, model_id = line.rstrip().split(',')
            try:
                alphafold[protein_acc] += 1
            except KeyError:
                alphafold[protein_acc] = 1

    if keep_fragments:
        return set(alphafold.keys())

    return {protein_acc for protein_acc, cnt in alphafold.items() if cnt == 1}
