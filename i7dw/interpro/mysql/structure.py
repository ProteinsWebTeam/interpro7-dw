import json

from . import reduce
from ... import dbms, logger, pdbe
from ...io import Store


def insert_structures(ora_uri, uri):
    structures = pdbe.get_structures(ora_uri)
    sec_structures = pdbe.get_secondary_structures(ora_uri)

    con, cur = dbms.connect(uri)
    cur.close()
    table = dbms.Populator(
        con=con,
        query="""
            INSERT INTO webfront_structure (
                  accession,
                  name,
                  source_database,
                  experiment_type,
                  release_date,
                  resolution,
                  literature,
                  chains,
                  proteins,
                  secondary_structures
              ) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
        """
    )

    for pdbe_id, s in structures.items():
        proteins = {}
        all_chains = set()

        for acc, chains in s["proteins"].items():
            proteins[acc] = []
            for chain_id in chains:
                all_chains.add(chain_id)
                proteins[acc].append(chain_id)

        table.insert((
            pdbe_id,
            s["name"],
            "pdb",
            s["evidence"],
            s["date"],
            s["resolution"],
            json.dumps(s["citations"]),
            json.dumps(sorted(all_chains)),
            json.dumps(proteins),
            json.dumps(sec_structures.get(pdbe_id, []))
        ))
    table.close()
    con.commit()
    con.close()


def get_structures(uri: str) -> dict:
    structures = {}
    con, cur = dbms.connect(uri)
    cur.execute(
        """
        SELECT accession, name, experiment_type, resolution, proteins
        FROM webfront_structure
        """
    )

    for acc, name, _type, resolution, proteins in cur:
        structures[acc] = {
            "accession": acc,
            "name": name,
            "type": _type,
            "resolution": resolution,
            "proteins": json.loads(proteins)
        }

    cur.close()
    con.close()

    return structures


def update_counts(uri: str, src_structures: str):
    con, cur = dbms.connect(uri)
    logger.info("updating webfront_structure")
    structures = set(get_structures(uri))
    with Store(src_structures) as store:
        for pdb_id, xrefs in store:
            try:
                structures.remove(pdb_id)
            except KeyError:
                continue

            counts = reduce(xrefs)
            try:
                counts["entries"]["total"] = sum(counts["entries"].values())
            except KeyError:
                counts["entries"] = {"total": 0}

            counts["proteins"] = xrefs["proteins"].pop()  # set of one item

            cur.execute(
                """
                UPDATE webfront_structure
                SET counts = %s
                WHERE accession = %s
                """,
                (json.dumps(counts), pdb_id)
            )

    for pdb_id in structures:
        cur.execute(
            """
            UPDATE webfront_structure
            SET counts = %s
            WHERE accession = %s
            """,
            (json.dumps({
                "domain_architectures": 0,
                "entries": {"total": 0},
                "proteins": 0,
                "proteomes": 0,
                "sets": 0,
                "taxa": 0,
            }), pdb_id)
        )

    con.commit()
    cur.close()
    con.close()

