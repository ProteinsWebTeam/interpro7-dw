import pickle

import cx_Oracle


MAIN_RANKS = [
    "superkingdom",
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species"
]


def export_taxa(url: str, file: str):
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    cur.execute(
        """
        SELECT TO_CHAR(TAX_ID), TO_CHAR(PARENT_ID), SCIENTIFIC_NAME, 
          FULL_NAME, RANK
        FROM INTERPRO.ETAXI
        """
    )

    taxa = {}
    for row in cur:
        taxon_id = row[0]

        taxa[taxon_id] = {
            "id": taxon_id,
            "parent": row[1],
            "sci_name": row[2],
            "full_name": row[3],
            "rank": row[4],
            "children": set(),
            "lineage": [taxon_id]
        }

    cur.close()
    con.close()

    for taxon_id, taxon in taxa.items():
        node_id = taxon_id
        parent_id = taxon["parent"]

        # Traverse lineage from child to parent
        while parent_id is not None:
            taxon["lineage"].append(parent_id)
            taxa[parent_id]["children"].add(node_id)

            # We move to the parent
            node_id = parent_id
            parent_id = taxa[parent_id]["parent"]

    for info in taxa.values():
        info["children"] = list(info["children"])
        info["lineage"] = list(map(str, reversed(info["lineage"])))

    """
    Define lineage but with major ranks only.
    Some ranks will stay empty (None, rank) because not all clades 
    exist (e.g. no family between an order and a genus).
    """
    for info in taxa.values():
        lineage = [[None, rank] for rank in MAIN_RANKS]

        for node_id in info["lineage"]:
            node = taxa[node_id]

            try:
                i = MAIN_RANKS.index(node["rank"])
            except ValueError:
                pass
            else:
                lineage[i][0] = node_id

        info["main_ranks"] = lineage

    with open(file, "wb") as fh:
        pickle.dump(taxa, fh)
