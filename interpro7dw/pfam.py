import re
from datetime import datetime, timezone

import oracledb

from interpro7dw import wikipedia
from interpro7dw.utils import logger
from interpro7dw.utils.store import BasicStore
from interpro7dw.utils.mysql import uri2dict


def get_details(uri: str) -> dict:
    con = oracledb.connect(uri)
    entries = {}
    with con.cursor() as cur:
        cur.execute(
            """
            SELECT 
            accession,
            sequence_ontology,
            hmm_build, hmm_search,
            seq_gathering, domain_gathering,
            version
            FROM INTERPRO.PFAM_DATA
            """
        )
    
        for row in cur:
            entries[row[0]] = {
                "curation": {
                    # "seed_source": row[1],
                    # "type": row[2],
                    "sequence_ontology": row[1],
                    "authors": [],
                    # "num_seed": row[4],  # Number in seed
                    # "num_full": row[5],  # Number in full
                    # "avg_length": row[6],  # Length of the domain
                    # "avg_id": row[7],  # Identity of full alignment
                    # "avg_coverage": row[8],  # Coverage of the seq by the domain
                },
                "hmm": {
                    "commands": {
                        "build": row[2],
                        "search": row[3]
                    },
                    "cutoffs": {
                        "gathering": {
                            "sequence": row[4],
                            "domain": row[5],
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
                    "version": row[6]
                }
            }

        cur.execute(
            """
            SELECT accession, author, orcid
            FROM INTERPRO.PFAM_AUTHOR    
            """
        )
        for acc, author, orcid in cur:
            entries[acc]["curation"]["authors"].append({
                "author": author,
                "orcid": orcid
            })

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
    """
    Return a tuple (to_change, results)
    where to_change is a list of Pfam entries that need an update:
        PF00001, [page_to_remove_1, page_to_remove_2], [page_to_add_1]

    and results is a dictionary of Pfam entries with their Wikipedia pages:
        Key: PF00001
        Value: [
                    {
                        title: page_1
                        extract: summary of page
                        thumbnail: thumbnail of first image
                    },
                    ...
               ]
    """
    # Pfam DB in LATIN1, with special characters in Wikipedia title
    logger.debug("loading Pfam entries")
    logger.error("**WIKI**")
    con = oracledb.connect(uri, use_unicode=False)
    with con.cursor() as cur:
        cur.execute(
            """
            SELECT W.accession, W.title, D.name
            FROM INTERPRO.PFAM_WIKIPEDIA W
            INNER JOIN INTERPRO.PFAM_DATA D ON W.accession = D.accession
            """
        )
        rows = cur.fetchall()
    con.close()

    # Pfam -> Wikipedia, in the Pfam database
    pfam_acc2wiki = {}
    key2acc = {}
    for pfam_acc, pfam_id, title in rows:  # Pfam_id == name in InterPro
        # pfam_acc = pfam_acc.decode("utf-8")
        # pfam_id = pfam_id.decode("utf-8")
        # try:
        #     title = title.decode("utf-8")
        # except UnicodeDecodeError:
        #     logger.critical(f"{pfam_acc}: {title}")
        #     raise

        """
        May contains special characters
        Some of these characters seem to be utf-8 interpreted as cp1252.
        e.g. en dash (–) returned as \xc3\xa2\xe2\x82\xac\xe2\x80\x9c
        >>> s = b"\xc3\xa2\xe2\x82\xac\xe2\x80\x9c"
        >>> s = s.decode('utf-8')
        'â€“'

        So we re-encode in cp1252, then decode in utf-8
        >>> s.decode('utf-8').encode('cp1252').decode('utf-8')
        '–'
        
        And if we have an encoding/decoding error... it was probably not cp1252
        """
        try:
            obj = title.encode("cp1252").decode("utf-8")
        except (UnicodeEncodeError, UnicodeDecodeError):
            pass
        else:
            title = obj

        # Canonicalize: replace spaces by underscores
        title = title.replace(" ", "_")

        if pfam_acc in pfam_acc2wiki:
            pfam_acc2wiki[pfam_acc].add(title)
        else:
            pfam_acc2wiki[pfam_acc] = {title}

        key2acc[pfam_acc.lower()] = pfam_acc
        key2acc[pfam_id.lower()] = pfam_acc

    # TODO: use links to pfam.xfam.org, pfam.org, and www.ebi.ac.uk/interpro
    # Pages containing external links to Pfam families
    # logger.debug("Finding external links to Pfam in Wikipedia articles")
    # pages = wikipedia.get_ext_links("pfam.xfam.org", validate=is_family_url)
    pages = {page for pages in pfam_acc2wiki.values() for page in pages}

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
    conn = oracledb.connect(uri)
    with conn.cursor() as cur:
        cur.execute(
            """
            SELECT D.CLAN_ID, D.DESCRIPTION, A.AUTHOR
            FROM PFAM_CLAN_DATA D
            INNER JOIN PFAM_CLAN_AUTHOR A ON D.CLAN_ID = A.CLAN_ID
            """
        )
        clans = {}
        for acc, desc, author in cur:
            try:
                clans[acc]['authors'].append(author)
            except KeyError:
                clans[acc] = {
                    "authors": [author],
                    "description": desc,
                    "literature": []
                }

        cur.execute(
            """
            SELECT CLAN_ID, PUBMED_ID, TITLE, AUTHOR, JOURNAL
            FROM INTERPRO.PFAM_CLAN_LITERATURE    
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

    conn.close()
    return clans


def export_alignments(uri: str, alignments_file: str):
    con = oracledb.connect(uri)
    with con.cursor() as cur:
        cur.execute(
            """
            SELECT D.ACCESSION,
                D.SEED_NUM, D.FULL_NUM,
                D.RP15_NUM, D.RP35_NUM,
                D.RP55_NUM, D.RP75_NUM,
                D.UNIPROT_NUM, A.TYPE, A.ALIGNMENT
            FROM INTERPRO.PFAM_DATA D
            INNER JOIN INTERPRO.PFAM_ALIGNMENTS A ON D.ACCESSION = A.ACCESSION
            """
        )
        counts = {}
        for row in cur:
            accession = row[0]
            counts[accession] = {
                "seed": row[1],
                "full": row[2],
                "rp15": row[3],
                "rp35": row[4],
                "rp55": row[5],
                "rp75": row[6],
                "uniprot": row[7]
            }
            aln_type = row[8]
            aln_bytes = row[9]

            try:
                count = counts[accession][aln_type]
            except KeyError:
                continue

            with BasicStore(alignments_file, "w") as store:
                store.write((
                    accession,
                    f"alignment:{aln_type}",
                    aln_bytes,  # gzip-compressed steam
                    count
                ))

    con.close()
