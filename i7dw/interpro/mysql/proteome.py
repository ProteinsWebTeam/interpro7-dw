import json

from . import reduce, taxonomy
from ... import dbms, logger, uniprot
from ...io import Store


def insert_proteomes(ora_uri, my_uri, chunk_size=100000):
    proteomes = uniprot.get_proteomes(ora_uri)
    taxa = set(taxonomy.get_taxa(my_uri, lineage=False))

    data = []
    con, cur = dbms.connect(my_uri)
    for p in proteomes.values():
        if p["tax_id"] not in taxa:
            """
            If tax_id not in taxa, it's very likely that INTERPRO.ETAXI
            (source for taxonomy table) is out-of-date
            """
            logger.warning("missing taxon (ID: {})".format(p["tax_id"]))
            continue

        data.append((
            p["accession"],
            p["name"],
            1 if p["is_reference"] else 0,
            p["strain"],
            p["assembly"],
            p["tax_id"]
        ))

    for i in range(0, len(data), chunk_size):
        cur.executemany(
            """
            INSERT INTO webfront_proteome (
              accession,
              name,
              is_reference,
              strain,
              assembly,
              taxonomy_id
            ) VALUES (
              %s, %s, %s, %s, %s, %s
            )
            """,
            data[i:i+chunk_size]
        )

    cur.close()
    con.commit()
    con.close()


def get_proteomes(uri: str) -> dict:
    con, cur = dbms.connect(uri)

    cur.execute(
        """
        SELECT accession, name, is_reference, strain, assembly, taxonomy_id
        FROM webfront_proteome
        """
    )

    proteomes = {}
    for row in cur:
        proteomes[row[0]] = {
            'name': row[1],
            'is_reference': bool(row[2]),
            'strain': row[3],
            'assembly': row[4],
            'taxon': row[5]
        }

    cur.close()
    con.close()

    return proteomes


def update_counts(uri: str, src_proteomes: str):
    con, cur = dbms.connect(uri)

    logger.info("updating webfront_proteome")
    proteomes = set(get_proteomes(uri))
    with Store(src_proteomes) as store:
        for upid, xrefs in store:
            try:
                proteomes.remove(upid)
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
                UPDATE webfront_proteome
                SET counts = %s
                WHERE accession = %s
                """,
                (json.dumps(counts), upid)
            )

    for upid in proteomes:
        cur.execute(
            """
            UPDATE webfront_proteome
            SET counts = %s
            WHERE accession = %s
            """,
            (json.dumps({
                "domain_architectures": 0,
                "entries": {"total": 0},
                "proteins": 0,
                "sets": 0,
                "structures": 0,
                "taxa": 1,
            }), upid)
        )

    con.commit()
    cur.close()
    con.close()
