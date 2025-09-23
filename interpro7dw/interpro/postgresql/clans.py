import pickle

from interpro7dw.utils import logger
from interpro7dw.utils.store import BasicStore
from .utils import connect, jsonify


def populate(uri: str, clans_file: str, clanxrefs_file: str):
    logger.info("loading clans")
    with open(clans_file, "rb") as fh:
        clans = pickle.load(fh)

    logger.info("creating webfront_set")
    con = connect(uri)
    cur = con.cursor()
    cur.execute("DROP TABLE IF EXISTS interpro.webfront_set")
    cur.execute(
        """
        CREATE TABLE interpro.webfront_set
        (
            accession VARCHAR(20) COLLATE "case_insensitive" PRIMARY KEY NOT NULL,
            name VARCHAR(400),
            description TEXT,
            source_database VARCHAR(10) NOT NULL,
            relationships JSONB NOT NULL,
            authors JSONB,
            literature JSONB,
            wikipedia JSONB,
            counts JSONB
        )
        """
    )

    query = """
        INSERT INTO interpro.webfront_set
        VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s)
    """

    with BasicStore(clanxrefs_file, mode="r") as store:
        params = []

        for accession, xrefs in store:
            # TODO: change this xref export
            xrefs["entries"]["total"] = xrefs["entries"].pop("all")

            clan = clans[accession]
            record = (
                accession,
                clan["name"],
                clan["description"],
                clan["database"].lower(),
                jsonify(clan["relationships"], nullable=False),
                jsonify(clan.get("authors", []), nullable=False),     # only Pfam
                jsonify(clan.get("literature", []), nullable=False),  # only Pfam
                jsonify(clan.get("wikipedia", []), nullable=False),   # only Pfam
                jsonify({
                    "domain_architectures": len(xrefs["dom_orgs"]),
                    "entries": {k.lower(): len(v)
                                for k, v in xrefs["entries"].items()},
                    "proteins": len(xrefs["proteins"]),
                    "proteomes": len(xrefs["proteomes"]),
                    "structures": len(xrefs["structures"]),
                    "taxa": len(xrefs["taxa"])
                })
            )

            params.append(record)
            if len(params) == 1000:
                cur.executemany(query, params)
                params.clear()

        if params:
            cur.executemany(query, params)
            params.clear()

    con.commit()
    cur.close()
    con.close()

    logger.info("done")
