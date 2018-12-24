#!/usr/bin/env python
# -*- coding: utf-8 -*-

import base64
import json
import os
import urllib.parse
import urllib.error
import urllib.request
from tempfile import mkstemp

from . import dbms, hmmer


def get_wiki(uri):
    base_url = 'https://en.wikipedia.org/api/rest_v1/page/summary/'

    # Pfam DB in LATIN1, with special charachters in Wikipedia title
    con, cur = dbms.connect(uri, encoding='latin1')
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
        # cursor returns bytes instead of string due to latin1
        acc = acc.decode()
        title = title.decode()

        try:
            url = base_url + urllib.parse.quote(title)
            res = urllib.request.urlopen(url)
        except urllib.error.HTTPError as e:
            # Content can be retrieved with e.fp.read()
            continue
        else:
            obj = json.loads(res.read().decode('utf-8'))
            thumbnail = obj.get('thumbnail')

            if thumbnail:
                try:
                    filename, headers = urllib.request.urlretrieve(thumbnail['source'])
                except urllib.error.ContentTooShortError:
                    b64str = None
                else:
                    with open(filename, 'rb') as fh:
                        b64str = base64.b64encode(fh.read()).decode('utf-8')

                    os.unlink(filename)
            else:
                b64str = None

            entries[acc] = {
                'title': title,
                'extract': obj['extract'],
                'extract_html': obj['extract_html'],
                'thumbnail': b64str
            }

    return entries


def get_annotations(uri):
    con, cur = dbms.connect(uri, sscursor=True)
    cur.execute(
        """
        SELECT a.pfamA_acc, a.alignment, h.hmm
        FROM alignment_and_tree a
        INNER JOIN pfamA_HMM h ON a.pfamA_acc = h.pfamA_acc AND a.type = 'seed'
        """
    )

    entries = {}
    for pfam_ac, aln, hmm in cur:
        if pfam_ac not in entries:
            entries[pfam_ac] = []

        if aln is not None:
            entries[pfam_ac].append({
                'type': 'alignment',
                'value': aln,
                'mime_type': 'application/octet-stream'
            })

        if hmm is not None:
            entries[pfam_ac].append({
                'type': 'hmm',
                'value': hmm,
                'mime_type': 'application/octet-stream'
            })

            fd, filename = mkstemp()
            os.close(fd)

            with open(filename, 'wb') as fh:
                fh.write(hmm)

            # hmm_logo = call_skylign(filename)
            hmm_logo = hmmer.hmm_to_logo(filename, method='info_content_all', processing='hmm')
            os.unlink(filename)

            entries[pfam_ac].append({
                'type': 'logo',
                'value': json.dumps(hmm_logo),
                'mime_type': 'application/json'
            })

    cur.close()
    con.close()

    annotations = []
    for pfam_ac in entries:
        for anno in entries[pfam_ac]:
            annotations.append((
                '{}--{}'.format(pfam_ac, anno['type']),
                pfam_ac,
                anno['type'],
                anno['value'],
                anno['mime_type']
            ))

    return annotations


def get_clans(uri) -> dict:
    con, cur = dbms.connect(uri, sscursor=True)
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
        clan_ac = row[0]

        if clan_ac not in clans:
            clans[clan_ac] = {
                "accession": clan_ac,
                "name": row[1],
                "description": row[2],
                "relationships": {
                    "nodes": [],
                    "links": {}
                }
            }

        clans[clan_ac]["relationships"]["nodes"].append({
            "accession": row[4],
            "type": "entry",
            "score": row[5] / row[3]
        })

    # # Not used any more, as we use profile-profile alignments
    # cur.execute(
    #     """
    #     SELECT
    #       m1.clan_acc, el.pfamA_acc_1,
    #       rel.pfamA_acc_2, rel.evalue
    #     FROM pfamA2pfamA_hhsearch rel
    #     INNER JOIN clan_membership m1
    #       ON rel.pfamA_acc_1 = m1.pfamA_acc
    #     INNER JOIN clan_membership m2
    #       ON rel.pfamA_acc_2 = m2.pfamA_acc
    #       AND m1.clan_acc = m2.clan_acc
    #     """
    # )
    #
    # for clan_ac, ac1, ac2, evalue in cur:
    #     links = clans[clan_ac]["relationships"]["links"]
    #
    #     if ac1 > ac2:
    #         ac1, ac2 = ac2, ac1
    #
    #     if ac1 not in links:
    #         links[ac1] = {ac2: evalue}
    #     elif ac2 not in links[ac1] or evalue < links[ac1][ac2]:
    #         links[ac1][ac2] = evalue

    cur.close()
    con.close()

    for clan in clans.values():
        clan["relationships"]["nodes"].sort(key=lambda x: x["accession"])
        # clan["relationships"]["links"] = [
        #     {
        #         "source": ac1,
        #         "target": ac2,
        #         "score": ev
        #     }
        #     for ac1, targets in clan["relationships"]["links"].items()
        #     for ac2, ev in targets.items()
        # ]

    return clans
