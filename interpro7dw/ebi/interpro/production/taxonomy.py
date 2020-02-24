# -*- coding: utf-8 -*-

import cx_Oracle

from interpro7dw.utils import datadump


def export_taxonomy(url: str, output: str):
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    cur.execute(
        """
        SELECT TO_CHAR(TAX_ID), TO_CHAR(PARENT_ID), SCIENTIFIC_NAME, 
          FULL_NAME, RANK
        FROM INTERPRO.ETAXI
        """
    )

    taxonomy = {}
    for row in cur:
        taxon_id = row[0]

        taxonomy[taxon_id] = {
            "parent": row[1],
            "sci_name": row[2],
            "full_name": row[3],
            "rank": row[4],
            "children": set(),
            "lineage": [taxon_id]
        }

    cur.close()
    con.close()

    for taxon_id, taxon in taxonomy.items():
        node_id = taxon_id
        parent_id = taxon["parent"]

        # Traverse lineage from child to parent
        while parent_id is not None:
            taxon["lineage"].append(parent_id)
            taxonomy[parent_id]["children"].add(node_id)

            # We move to the parent
            node_id = parent_id
            parent_id = taxonomy[parent_id]["parent"]

    for taxon_id, info in taxonomy.items():
        info["children"] = list(info["children"])
        info["lineage"] = list(map(str, reversed(info["lineage"])))

    datadump(output, taxonomy)
