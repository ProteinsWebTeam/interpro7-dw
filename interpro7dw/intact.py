"""
Module to parse the PSI-MI TAB file produced by IntAct, and extract interactions
relevant to InterPro.

You can download the file from the EMBL-EBI FTP site:
https://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.txt

The file is basically a TSV, but InterPro accessions can be found as:

- Interaction identifier (column 13; 0-based index)
- Interactor xref (columns 22 and 23)
- Interaction xref (column 24)
- Interactor annotation (columns 25 and 26)
- Interaction annotation (column 27)
- Interaction feature (columns 36 and 37)

Each column can contain several values, separated by a pipe (|).

"""

import re


def get_interpro_interactions(file: str) -> dict[str, list]:
    entries = {}

    with open(file, "rt") as fh:
        keys = next(fh).rstrip().split("\t")
        check_num_columns(keys)

        for line in map(str.rstrip, fh):
            values = line.split("\t")
            check_num_columns(values)

            # Find all InterPro references
            accessions = set(re.findall(r"IPR\d{6}", "\t".join(line)))

            if not accessions:
                continue

            # Interaction identifier(s)
            interaction_id = find_interaction(values[13])

            acc1, name1 = find_interactor(
                values[0],  # ID(s) interactor A
                values[4],  # Alias(es) interactor A
            )

            acc2, name2 = find_interactor(
                values[1],  # ID(s) interactor B
                values[5],  # Alias(es) interactor B
            )

            pmid = find_pmid(values[8])  # Publication Identifier(s)

            if (
                interaction_id is not None
                and acc1 is not None
                and acc2 is not None
                and pmid is not None
            ):
                obj = (interaction_id, acc1, name1, acc2, name2, pmid)
                for entry_acc in accessions:
                    try:
                        entries[entry_acc].add(obj)
                    except KeyError:
                        entries[entry_acc] = {obj}

    return {
        entry_acc: list(interactions) for entry_acc, interactions in entries.items()
    }


def check_num_columns(values: list[str]):
    if len(values) != 42:
        raise ValueError(f"Invalid format: expecting 42 columns, " f"got {len(values)}")


def find_interactor(identifiers: str, aliases: str) -> tuple[str | None, str | None]:
    accession = name = None
    m = re.search(r"uniprotkb:([A-Z0-9]+)", identifiers, flags=re.I)
    if m:
        accession = m.group(1)

    m = re.search(r"uniprotkb:([^(|]+)", aliases, flags=re.I)
    if m:
        name = m.group(1)

    return accession, name


def find_pmid(identifiers: str) -> int | None:
    m = re.search(r"pubmed:(\d+)", identifiers, flags=re.I)
    return int(m.group(1)) if m else None


def find_interaction(identifiers: str) -> str | None:
    m = re.search(r"intact:(EBI-\d+)", identifiers, flags=re.I)
    return m.group(1) if m else None
