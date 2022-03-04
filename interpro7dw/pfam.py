import re
from datetime import datetime, timezone

import MySQLdb
import MySQLdb.cursors

from interpro7dw import wikipedia
from interpro7dw.utils import logger
from interpro7dw.utils.store import BasicStore
from interpro7dw.utils.mysql import uri2dict


def get_details(uri: str) -> dict:
    con = MySQLdb.connect(**uri2dict(uri))
    cur = con.cursor()
    cur.execute(
        """
        SELECT 
            e.pfamA_acc, e.seed_source, e.type, so.so_id, e.num_seed, 
            e.num_full, e.average_length, e.percentage_id, e.average_coverage,
            e.buildMethod, e.searchMethod, e.sequence_GA, e.domain_GA, 
            e.sequence_TC, e.domain_TC, e.sequence_NC, e.domain_NC, 
            e.model_length, e.version
        FROM pfamA e
        INNER JOIN sequence_ontology so on e.type = so.type
        """
    )
    entries = {}
    for row in cur:
        entries[row[0]] = {
            "curation": {
                # "seed_source": row[1],
                # "type": row[2],
                "sequence_ontology": row[3],
                "authors": [],
                # "num_seed": row[4],  # Number in seed
                # "num_full": row[5],  # Number in full
                # "avg_length": row[6],  # Length of the domain
                # "avg_id": row[7],  # Identity of full alignment
                # "avg_coverage": row[8],  # Coverage of the seq by the domain
            },
            "hmm": {
                "commands": {
                    "build": row[9],
                    "search": row[10]
                },
                "cutoffs": {
                    "gathering": {
                        "sequence": row[11],
                        "domain": row[12],
                    },
                    # "trusted": {
                    #     "sequence": row[13],
                    #     "domain": row[14],
                    # },
                    # "noise": {
                    #     "sequence": row[15],
                    #     "domain": row[16],
                    # },
                },
                # "length": row[17],
                # "version": row[18]
            }
        }

    cur.execute(
        """
        SELECT e.pfamA_acc, a.author, a.orcid
        FROM pfamA_author e
        INNER JOIN author a on e.author_id = a.author_id
        ORDER BY e.author_rank        
        """
    )
    for acc, author, orcid in cur:
        entries[acc]["curation"]["authors"].append({
            "author": author,
            "orcid": orcid
        })

    cur.close()
    con.close()
    return entries


def is_family_url(url: str) -> bool:
    """Checks that a given URL points to a Pfam family page, e.g.:
      - http://pfam.xfam.org/family/PF16122
      - http://pfam.xfam.org/family?acc=PF16122
      - http://pfam.xfam.org/family/DIM

    URL with multiple slashes (e.g. http://pfam.xfam.org//family/PF16122)
    are accepted.

    :param url: Pfam URL
    :return: boolean
    """
    pattern = r"http://pfam\.xfam\.org/+family(/+|\?acc=)([a-z0-9_\-]+)/*"
    return re.fullmatch(pattern, url, flags=re.I) is not None


def is_pfam_infobox(name: str, value: str) -> bool:
    return "pfam" in name.lower()


def get_wiki(uri: str, hours: int = 0) -> tuple[list[tuple[str,
                                                           list[str],
                                                           list[str]]],
                                                dict[str, list[dict]]]:
    # Pfam DB in LATIN1, with special characters in Wikipedia title
    logger.debug("loading Pfam entries")
    con = MySQLdb.connect(**uri2dict(uri), use_unicode=False)
    cur = con.cursor()
    cur.execute(
        """
        SELECT p.pfamA_acc, p.pfamA_id, w.title
        FROM pfamA p
        INNER JOIN pfamA_wiki pw ON p.pfamA_acc = pw.pfamA_acc 
        INNER JOIN wikipedia w ON pw.auto_wiki = w.auto_wiki
        """
    )
    rows = cur.fetchall()
    cur.close()
    con.close()

    # Pfam -> Wikipedia, in the Pfam database
    pfam_acc2wiki = {}
    key2acc = {}
    for pfam_acc, pfam_id, title in rows:
        # cursor returns bytes instead of string due to `use_unicode=False`
        pfam_acc = pfam_acc.decode("utf-8")
        pfam_id = pfam_id.decode("utf-8")
        try:
            title.decode("ascii")
        except UnicodeDecodeError:
            """
            Contains special characters
            As of Mar 2022, these characters seem to be utf-8 
            interpreted as cp1252.
            e.g. en dash (–) returned as \xc3\xa2\xe2\x82\xac\xe2\x80\x9c
            >>> s = b"\xc3\xa2\xe2\x82\xac\xe2\x80\x9c"
            >>> s = s.decode('utf-8')
            'â€“'
            
            So we re-encode in cp1252, then decode in utf-8
            >>> s.decode('utf-8').encode('cp1252').decode('utf-8')
            '–'
            """
            title = title.decode("utf-8").encode("cp1252").decode("utf-8")
        else:
            # No special characters
            title = title.decode("utf-8")

        # Canonicalize: replace spaces by underscores
        title = title.replace(" ", "_")

        if pfam_acc in pfam_acc2wiki:
            pfam_acc2wiki[pfam_acc].add(title)
        else:
            pfam_acc2wiki[pfam_acc] = {title}

        key2acc[pfam_acc.lower()] = pfam_acc
        key2acc[pfam_id.lower()] = pfam_acc

    # Pages containing external links to Pfam families
    logger.debug("Finding external links to Pfam in Wikipedia articles")
    pages = wikipedia.get_ext_links("pfam.xfam.org", validate=is_family_url)

    # Pfam -> Wikipedia, from Wikipedia API
    wiki_acc2wiki = {}
    for title in pages:
        # Canonicalize: replace spaces by underscores
        title = title.replace(" ", "_")

        props = wikipedia.parse_infobox(title, validate=is_pfam_infobox)
        for name in props:
            for value in map(str.lower, props[name]):
                if value in key2acc:
                    pfam_acc = key2acc[value]
                    if pfam_acc in wiki_acc2wiki:
                        wiki_acc2wiki[pfam_acc].add(title)
                    else:
                        wiki_acc2wiki[pfam_acc] = {title}

    # Diff
    to_change = []
    to_fetch = set()
    for pfam_acc, pages in wiki_acc2wiki.items():
        new_pages = sorted(pages - pfam_acc2wiki.get(pfam_acc, set()))
        old_pages = sorted(pfam_acc2wiki.get(pfam_acc, set()) - pages)

        if new_pages or old_pages:
            to_change.append((pfam_acc, old_pages, new_pages))

        to_fetch |= pages

    logger.debug("Fetching Wikipedia pages")
    now = datetime.now(timezone.utc)
    parsed_page = {}
    for title in to_fetch:
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

        parsed_page[title] = {
            "title": title,
            # "extract": summary["extract"],
            "extract": summary["extract_html"],
            "thumbnail": wikipedia.get_thumbnail(summary)
        }

    result = {}
    for pfam_acc, pages in wiki_acc2wiki.items():
        pfam_pages = []
        for title in sorted(pages):
            if title in parsed_page:
                pfam_pages.append(parsed_page[title])

        if pfam_pages:
            result[pfam_acc] = pfam_pages
        else:
            logger.info(f"{pfam_acc}: could not retrieve articles")

    return to_change, result


def get_clans(uri: str) -> dict:
    con = MySQLdb.connect(**uri2dict(uri))
    cur = con.cursor()
    cur.execute(
        """
        SELECT clan_acc, clan_author, clan_comment
        FROM clan
        """
    )
    clans = {}
    for acc, authors, comment in cur:
        if authors:
            # Split on commas to make a list
            authors = list(map(str.strip, authors.split(',')))

        clans[acc] = {
            "authors": authors,
            "description": comment,
            "literature": []
        }

    cur.execute(
        """
        SELECT clr.clan_acc, lr.pmid, lr.title, lr.author, lr.journal
        FROM clan_lit_ref clr
        INNER JOIN literature_reference lr on clr.auto_lit = lr.auto_lit
        ORDER BY clr.clan_acc, clr.order_added        
        """
    )

    for acc, pmid, title, authors, journal in cur:
        if authors:
            # Trim trailing semi-colon and spaces
            authors = re.sub(r";\s*$", '', authors)

            # Split on commas to make a list
            authors = list(map(str.strip, authors.split(',')))

        clans[acc]["literature"].append({
            "PMID": pmid,
            "title": title.strip() if title else None,
            "authors": authors,
            "journal": journal.strip() if journal else None
        })

    cur.close()
    con.close()
    return clans


def export_alignments(uri: str, alignments_file: str):
    con = MySQLdb.connect(**uri2dict(uri))
    cur = MySQLdb.cursors.SSCursor(con)
    cur.execute(
        """
        SELECT pfamA_acc, num_seed, num_full, number_rp15, number_rp35, 
               number_rp55, number_rp75, number_uniprot
        FROM pfamA
        """
    )
    counts = {}
    for row in cur:
        counts[row[0]] = {
            "seed": row[1],
            "full": row[2],
            "rp15": row[3],
            "rp35": row[4],
            "rp55": row[5],
            "rp75": row[6],
            "uniprot": row[7]
        }

    with BasicStore(alignments_file, "w") as store:
        cur.execute(
            """
            SELECT pfamA_acc, type, alignment
            FROM alignment_and_tree
            WHERE alignment IS NOT NULL
            """
        )

        for accession, aln_type, aln_bytes in cur:
            try:
                count = counts[accession][aln_type]
            except KeyError:
                continue

            store.write((
                accession,
                f"alignment:{aln_type}",
                aln_bytes,  # gzip-compressed steam
                count
            ))

    cur.close()
    con.close()
