import json
import re
from datetime import datetime, timezone

import oracledb

from interpro7dw import wikipedia
from interpro7dw.utils import logger
from interpro7dw.utils.oracle import lob_as_str
from interpro7dw.utils.store import BasicStore


def get_details(uri: str) -> dict:
    con = oracledb.connect(uri)
    entries = {}
    with con.cursor() as cur:
        cur.execute(
            """
            SELECT ACCESSION, SEQ_ONTOLOGY_ID, AUTHORS, BUILD_CMD, SEARCH_CMD,
                   SEQ_GA, DOM_GA, VERSION
            FROM INTERPRO.PFAM_A
            """
        )

        for row in cur:
            entries[row[0]] = {
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
                }
            }

    con.close()
    return entries


def get_wiki(uri: str, hours: int = 0) -> dict[str, list[dict]]:
    """Get Wikipedia articles linked by Pfam entries.

    :param uri: InterPro Oracle connection string
    :param hours: Threshold for the minimum number of hours since
                  the latest edit of Wikipedia articles to be considered
    :return: A dictionary where keys are Pfam accessions and values are lists
             of Wikipedia articles (as dict)
    """
    logger.debug("loading Pfam entries")
    pfam2wiki = {}
    pages = {}
    con = oracledb.connect(uri)
    with con.cursor() as cur:
        cur.execute(
            """
            SELECT ACCESSION, WIKIPEDIA
            FROM INTERPRO.PFAM_A
            """
        )
        for pfam_acc, articles in cur.fetchall():
            for title in json.loads(articles):
                # Canonicalize: replace spaces by underscores
                title = title.replace(" ", "_")

                try:
                    pfam2wiki[pfam_acc].add(title)
                except KeyError:
                    pfam2wiki[pfam_acc] = {title}

                pages[title] = None

    con.close()

    logger.debug("Fetching Wikipedia pages")
    now = datetime.now(timezone.utc)
    for title in pages:
        summary = wikipedia.get_summary(title)
        if not summary:
            logger.error(f"{title}: could not retrieve summary")
            continue

        # e.g. 2020-04-14T10:10:52Z (UTC)
        timestamp = summary["timestamp"]

        last_rev = datetime.strptime(timestamp, "%Y-%m-%dT%H:%M:%SZ")
        last_rev = last_rev.replace(tzinfo=timezone.utc)

        hours_since_last_edit = (now - last_rev).total_seconds() / 3600
        if hours and hours_since_last_edit < hours:
            logger.warning(f"{title}: skipped (edited "
                           f"less than {hours} hours ago)")
            continue

        pages[title] = {
            "title": title,
            # "extract": summary["extract"],
            "extract": summary["extract_html"],
            "thumbnail": wikipedia.get_thumbnail(summary)
        }

    for pfam_acc in pfam2wiki:
        _pages = []
        for title in sorted(pfam2wiki[pfam_acc]):
            info = pages.get(title)
            if info:
                _pages.append(info)

        pfam2wiki[pfam_acc] = _pages

    return pfam2wiki


def get_clans(uri: str) -> dict:
    con = oracledb.connect(uri)
    with con.cursor() as cur:
        cur.execute(
            """
            SELECT ACCESSION, ABSTRACT, AUTHORS, REFERENCES
            FROM INTERPRO.PFAM_C
            """
        )
        clans = {}
        for accession, abstract, authors, references in cur:
            clans[accession] = {
                "authors": json.loads(authors),
                "description": abstract,
                "literature": json.loads(references)
            }

    con.close()
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
