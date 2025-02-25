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


def get_interpro_interactions(file: str) -> dict[str, list[dict]]:
    entries = {}

    with open(file, "rt") as fh:
        keys = next(fh).rstrip().split("\t")
        check_num_columns(keys)

        for line in map(str.rstrip, fh):
            values = line.split("\t")
            check_num_columns(values)

            # Find all InterPro references
            accessions = set(re.findall(r"IPR\d{6}", line))

            if not accessions:
                continue

            # Interaction identifier(s)
            interaction_id = find_interaction(values[13])

            acc1, name1, type1 = find_interactor(
                values[0],  # ID(s) interactor A
                values[4],  # Alias(es) interactor A
                values[20],  # Type(s) interactor A
            )

            acc2, name2, type2 = find_interactor(
                values[1],  # ID(s) interactor B
                values[5],  # Alias(es) interactor B
                values[21],  # Type(s) interactor B
            )

            pmid = find_pmid(values[8])  # Publication Identifier(s)

            if (
                interaction_id
                and acc1
                and name1
                and type1
                and acc2
                and name2
                and type2
                and pmid
            ):
                obj = {
                    "intact_id": interaction_id,
                    "pubmed_id": pmid,
                    "molecule_1": {
                        "accession": acc1,
                        "identifier": name1,
                        "type": type1,
                    },
                    "molecule_2": {
                        "accession": acc2,
                        "identifier": name2,
                        "type": type2,
                    },
                }
                for entry_acc in accessions:
                    try:
                        entries[entry_acc].append(obj)
                    except KeyError:
                        entries[entry_acc] = [obj]

    return entries


def check_num_columns(values: list[str]):
    if len(values) != 42:
        raise ValueError(f"Invalid format: expecting 42 columns, " f"got {len(values)}")


def find_interactor(
    identifiers: str, aliases: str, types: str
) -> tuple[str | None, str | None, str | None]:
    accession = name = _type = None
    m = re.search(r"uniprotkb:([A-Z0-9]+)", identifiers, flags=re.I)
    if m:
        accession = m.group(1)

    # Prefer the PSI-MI long name
    m = re.search(r"psi-mi:([^(]+)\(display_long\)", aliases, flags=re.I)
    if m:
        name = m.group(1).upper()
    else:
        # Fallback to whatever comes for UniprotKB (like gene name)
        m = re.search(r"uniprotkb:([^(|]+)", aliases, flags=re.I)
        if m:
            name = m.group(1)

    l = types.split("|")
    if len(l) == 1:
        m = re.search(r'psi-mi:"MI:\d+"\(([^)]+)\)', types, flags=re.I)
        _type = m.group(1) if m else None

    return accession, name, _type


def find_pmid(identifiers: str) -> int | None:
    m = re.search(r"pubmed:(\d+)", identifiers, flags=re.I)
    return int(m.group(1)) if m else None


def find_interaction(identifiers: str) -> str | None:
    m = re.search(r"intact:(EBI-\d+)", identifiers, flags=re.I)
    return m.group(1) if m else None
