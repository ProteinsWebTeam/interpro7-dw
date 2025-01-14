import json
import pickle
from datetime import datetime, timezone

import oracledb

from interpro7dw import wikipedia
from interpro7dw.utils import logger
from interpro7dw.utils.oracle import lob_as_str
from interpro7dw.utils.store import BasicStore


def export_families(uri: str, output: str):
    logger.debug("loading Pfam families")
    con = oracledb.connect(uri)
    entries = {}
    with con.cursor() as cur:
        cur.execute(
            """
            SELECT ACCESSION, SEQ_ONTOLOGY_ID, AUTHORS, BUILD_CMD, SEARCH_CMD,
                   SEQ_GA, DOM_GA, VERSION, WIKIPEDIA
            FROM INTERPRO.PFAM_A
            """
        )

        for row in cur:
            entries[row[0]] = {
                "details": {
                    "curation": {
                        "sequence_ontology": row[1],
                        "authors": json.loads(row[2]),
                    },
                    "hmm": {
                        "commands": {
                            "build": row[3],
                            "search": row[4]
                        },
                        "cutoffs": {
                            "gathering": {
                                "sequence": row[5],
                                "domain": row[6],
                            },
                        },
                        "version": row[7]
                    },
                },
                "wikipedia": json.loads(row[7])
            }

    con.close()

    logger.debug("Fetching Wikipedia pages")
    for info in entries.values():
        pages = []

        for title in sorted(info["wikipedia"]):
            obj = get_wiki(title)
            if obj:
                pages.append(obj)

        info["wikipedia"] = pages

    with open(output, "wb") as fh:
        pickle.dump(entries, fh)

    logger.info("done")


def get_wiki(title: str, min_hours: int = 0) -> dict | None:
    """Get Wikipedia article for a given page title.

    :param title: Title of the Wikipedia page
    :param min_hours: Threshold for the minimum number of hours since
                  the latest edit of Wikipedia articles to be considered
    :return: A dict of the page details or None
    """
    # Canonicalize: replace spaces by underscores
    title = title.replace(" ", "_")
    summary = wikipedia.get_summary(title)
    if not summary:
        logger.error(f"{title}: could not retrieve summary")
        return None

    # e.g. 2020-04-14T10:10:52Z (UTC)
    timestamp = summary["timestamp"]

    last_rev = datetime.strptime(timestamp, "%Y-%m-%dT%H:%M:%SZ")
    last_rev = last_rev.replace(tzinfo=timezone.utc)

    now = datetime.now(timezone.utc)
    hours_since_last_edit = (now - last_rev).total_seconds() / 3600
    if min_hours and hours_since_last_edit < min_hours:
        logger.warning(f"{title}: skipped (edited "
                       f"less than {min_hours} hours ago)")
        return None

    return {
        "title": title,
        # "extract": summary["extract"],
        "extract": summary["extract_html"],
        "thumbnail": wikipedia.get_thumbnail(summary)
    }


def get_clans(uri: str) -> dict:
    logger.debug("loading Pfam clans")
    con = oracledb.connect(uri)
    with con.cursor() as cur:
        cur.execute(
            """
            SELECT ACCESSION, ABSTRACT, AUTHORS, REFERENCES, WIKIPEDIA
            FROM INTERPRO.PFAM_C
            """
        )
        clans = {}
        for accession, abstract, authors, references, wiki_articles in cur:
            clans[accession] = {
                "authors": json.loads(authors),
                "description": abstract,
                "literature": json.loads(references),
                "wikipedia": json.loads(wiki_articles)
            }

    con.close()

    logger.debug("Fetching Wikipedia pages")
    for info in clans.values():
        pages = []

        for title in sorted(info["wikipedia"]):
            obj = get_wiki(title)
            if obj:
                pages.append(obj)

        info["wikipedia"] = pages

    return clans


def export_alignments(uri: str, alignments_file: str):
    con = oracledb.connect(uri)
    with con.cursor() as cur, BasicStore(alignments_file, "w") as bs:
        cur.outputtypehandler = lob_as_str
        cur.execute(
            """
            SELECT ACCESSION, SEED_ALN, SEED_NUM, FULL_ALN, FULL_NUM
            FROM INTERPRO.PFAM_A
            """
        )

        for accession, seed_aln, seed_num, full_aln, full_num in cur:
            bs.write((
                accession,
                f"alignment:seed",
                seed_aln,  # gzip-compressed steam
                seed_num
            ))
            bs.write((
                accession,
                f"alignment:full",
                full_aln,  # gzip-compressed steam
                full_num
            ))

    con.close()
