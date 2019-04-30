import json

from . import reduce
from .. import oracle
from ... import dbms, logger
from ...io import KVdb, Store


def insert_taxa(ora_uri, my_uri, chunk_size=100000):
    taxa = oracle.get_taxa(ora_uri)

    con, cur = dbms.connect(my_uri)

    data = [(
        taxon['id'],
        taxon['sci_name'],
        taxon['full_name'],
        # leading/trailing whitespaces are important from API queries
        ' {} '.format(' '.join(taxon['lineage'])),
        taxon['parent_id'],
        taxon['rank'],
        json.dumps(taxon['children'])
    ) for taxon in taxa]

    for i in range(0, len(data), chunk_size):
        cur.executemany(
            """
            INSERT INTO webfront_taxonomy (
                accession,
                scientific_name,
                full_name,
                lineage,
                parent_id,
                rank,
                children
            ) VALUES (
              %s, %s, %s, %s, %s, %s, %s
            )
            """,
            data[i:i+chunk_size]
        )

    cur.close()
    con.commit()
    con.close()


def get_taxa(uri: str, lineage: bool=False) -> dict:
    con, cur = dbms.connect(uri)
    if lineage:
        cur.execute(
            """
            SELECT accession, scientific_name, full_name, lineage, rank
            FROM webfront_taxonomy
            """
        )
        cols = ("taxId", "scientificName", "fullName", "lineage", "rank")
        taxa = {row[0]: dict(zip(cols, row)) for row in cur}
    else:
        cur.execute(
            """
            SELECT accession, scientific_name, full_name
            FROM webfront_taxonomy
            """
        )
        cols = ("taxId", "scientificName", "fullName")
        taxa = {row[0]: dict(zip(cols, row)) for row in cur}

    cur.close()
    con.close()

    return taxa


def update_counts(uri: str, src_taxa: str, tmpdir: str=None):
    with KVdb(cache_size=10000, dir=tmpdir) as taxa:
        logger.info("loading taxa")
        with Store(src_taxa) as store:
            for tax_id, xrefs in store:
                for e in xrefs["proteins"]:
                    # xrefs["proteins"] is a set of one item
                    xrefs["proteins_total"] = e
                taxa[tax_id] = xrefs

        logger.info("propagating cross-references to taxa lineage")
        all_taxa = set()
        cnt = 0
        for tax_id, t in get_taxa(uri, lineage=True).items():
            all_taxa.add(tax_id)

            try:
                taxon = taxa[tax_id]
            except KeyError:
                continue

            n_proteins = taxon["proteins"].pop()

            # lineage stored as a string in MySQL (string include the taxon)
            # -2: first item to include (second to last; last is current taxon)
            # -1: negative step (reverse list)
            lineage = t["lineage"].strip().split()[-2::-1]
            for parent_id in lineage:
                try:
                    parent = taxa[parent_id]
                except KeyError:
                    parent = {
                        "domain_architectures": set(),
                        "entries": {},
                        "proteomes": set(),
                        "proteins": {0},
                        "proteins_total": 0,
                        "sets": set(),
                        "structures": set()
                    }

                parent["proteins_total"] += n_proteins

                for _type in ("domain_architectures", "proteomes", "sets", "structures"):
                    try:
                        accessions = taxon[_type]
                    except KeyError:
                        # Type absent in taxon (e.g. no cross-refs)
                        accessions = taxon[_type] = set()
                    finally:
                        if _type in parent:
                            parent[_type] |= set(accessions)
                        else:
                            parent[_type] = set(accessions)

                try:
                    entries = taxon["entries"]
                except KeyError:
                    entries = taxon["entries"] = {}
                finally:
                    if "entries" not in parent:
                        parent["entries"] = {}

                    for entry_db, db_entries in entries.items():
                        if entry_db in parent["entries"]:
                            parent["entries"][entry_db] |= set(db_entries)
                        else:
                            parent["entries"][entry_db] = set(db_entries)

                # Write back parent to DB
                taxa[parent_id] = parent

            # Write back taxon to DB
            taxa[tax_id] = taxon

            cnt += 1
            if not cnt % 100000:
                logger.info("{:>12,}".format(cnt))

        logger.info("{:>12,}".format(cnt))

        logger.info("updating webfront_taxonomy")
        con, cur = dbms.connect(uri)
        for tax_id, taxon in taxa:
            all_taxa.remove(tax_id)
            counts = reduce(taxon)
            counts["proteins"] = counts.pop("proteins_total")

            try:
                counts["entries"]["total"] = sum(counts["entries"].values())
            except KeyError:
                counts["entries"] = {"total": 0}

            cur.execute(
                """
                UPDATE webfront_taxonomy
                SET counts = %s
                WHERE accession = %s
                """,
                (json.dumps(counts), tax_id)
            )

        logger.info("database size: {:,}".format(taxa.getsize()))

    for tax_id in all_taxa:
        cur.execute(
            """
            UPDATE webfront_taxonomy
            SET counts = %s
            WHERE accession = %s
            """,
            (json.dumps({
                "domain_architectures": 0,
                "entries": {"total": 0},
                "proteomes": 0,
                "proteins": 0,
                "sets": 0,
                "structures": 0
            }), tax_id)
        )

    con.commit()
    cur.close()
    con.close()
