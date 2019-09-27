# -*- coding: utf-8 -*-

import json
import os
import time
from base64 import b64encode
from http.client import IncompleteRead
from tempfile import mkstemp
from typing import Dict, List
from urllib.error import HTTPError
from urllib.parse import quote, unquote
from urllib.request import urlopen

import MySQLdb
import MySQLdb.cursors

from i7dw import hmmer, logger
from i7dw.interpro.mysql.utils import parse_url


def get_wiki(url: str, max_retries: int=4) -> Dict[str, Dict]:
    base_url = "https://en.wikipedia.org/api/rest_v1/page/summary/"

    # Pfam DB in LATIN1, with special characters in Wikipedia title
    con = MySQLdb.connect(**parse_url(url), use_unicode=False)
    cur = con.cursor()
    cur.execute(
        """
        SELECT p.pfamA_acc, w.title
        FROM pfamA_wiki p
        INNER JOIN wikipedia w ON p.auto_wiki = w.auto_wiki
        """
    )
    rows = cur.fetchall()
    cur.close()
    con.close()

    entries = {}
    for acc, title in rows:
        # cursor returns bytes instead of string due to `use_unicode=False`
        acc = acc.decode()
        title = title.decode()

        # Some records contains HTML %xx escapes: we need to replace them
        title = unquote(title)

        # default `safe` is '/' but we *want* to replace it
        url = base_url + quote(title, safe='')

        obj = None
        num_retries = 0
        while True:
            try:
                res = urlopen(url)
                data = res.read()
            except HTTPError as e:
                # Content can be retrieved with e.fp.read()
                logger.error(f"{title}: {e.code} ({e.reason})")
                break
            except IncompleteRead:
                if num_retries == max_retries:
                    logger.error(f"{title}: incomplete")
                    break
                else:
                    num_retries += 1
                    time.sleep(2)
            else:
                obj = json.loads(data.decode("utf-8"))
                break

        if obj is None:
            # Failed to get data
            continue

        entries[acc] = {
            "title": title,
            # "extract": obj["extract"],
            "extract": obj["extract_html"],
            "thumbnail": None
        }

        thumbnail = obj.get("thumbnail")
        if thumbnail:
            num_retries = 0
            while True:
                try:
                    res = urlopen(thumbnail["source"])
                    data = res.read()
                except HTTPError as e:
                    logger.error(f"{title} (thumbnail): "
                                 f"{e.code} ({e.reason})")
                    break
                except IncompleteRead:
                    if num_retries == max_retries:
                        logger.error(f"{title} (thumbnail): incomplete")
                        break
                    else:
                        num_retries += 1
                        time.sleep(2)
                else:
                    entries[acc]["thumbnail"] = b64encode(data).decode("utf-8")
                    break

    return entries


def get_annotations(url: str) -> Dict[str, List[Dict]]:
    con = MySQLdb.connect(**parse_url(url))
    cur = MySQLdb.cursors.SSCursor(con)
    cur.execute(
        """
        SELECT a.pfamA_acc, a.alignment, h.hmm
        FROM alignment_and_tree a
        INNER JOIN pfamA_HMM h 
          ON a.pfamA_acc = h.pfamA_acc AND a.type = 'seed'
        """
    )

    annotations = {}
    for acc, aln, hmm in cur:
        try:
            entry = annotations[acc]
        except KeyError:
            entry = annotations[acc] = []

        if aln is not None:
            entry.append({
                "type": "alignment",
                "value": aln,
                "mime_type": "application/octet-stream"
            })

        if hmm is not None:
            entry.append({
                "type": "hmm",
                "value": hmm,
                "mime_type": "application/octet-stream"
            })

            fd, filename = mkstemp()
            os.close(fd)
            with open(filename, "wb") as fh:
                fh.write(hmm)

            # hmm_logo = call_skylign(filename)
            hmm_logo = hmmer.hmm_to_logo(filename, method="info_content_all",
                                         processing="hmm")
            os.unlink(filename)

            entry.append({
                "type": "logo",
                "value": json.dumps(hmm_logo),
                "mime_type": "application/json"
            })

    cur.close()
    con.close()

    return annotations


def get_clans(url) -> Dict[str, Dict]:
    con = MySQLdb.connect(**parse_url(url))
    cur = MySQLdb.cursors.SSCursor(con)
    cur.execute(
        """
        SELECT
          c.clan_acc, c.clan_id, c.clan_description,
          c.number_sequences, m.pfamA_acc, f.num_full
        FROM clan c
        INNER JOIN clan_membership m ON c.clan_acc = m.clan_acc
        INNER JOIN pfamA f ON m.pfamA_acc = f.pfamA_acc
        """
    )

    clans = {}
    for row in cur:
        clan_acc = row[0]

        try:
            clan = clans[clan_acc]
        except KeyError:
            clan = clans[clan_acc] = {
                "accession": clan_acc,
                "name": row[1],
                "description": row[2],
                "relationships": {
                    "nodes": [],
                    "links": {}
                }
            }

        clan["relationships"]["nodes"].append({
            "accession": row[4],
            "type": "entry",
            "score": row[5] / row[3]
        })

    cur.close()
    con.close()

    for clan in clans.values():
        clan["relationships"]["nodes"].sort(key=lambda x: x["accession"])

    return clans
