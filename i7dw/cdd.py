# -*- coding: utf-8 -*-

import gzip
import os
import re
import urllib.request


_CDDID = "ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/cddid_all.tbl.gz"
_LINKS = "ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/family_superfamily_links"
_CDTRACK = "ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/cdtrack.txt"


def get_superfamilies() -> dict:
    filename, headers = urllib.request.urlretrieve(_CDDID)

    sets = {}
    nodes = {}
    parent_of = {}
    domains = {}
    with gzip.open(filename, "rt") as fh:
        for line in fh:
            obj = line.rstrip().split("\t")
            pssm_id = obj[0]
            accession = obj[1]
            short_name = obj[2]
            descr = obj[3].lstrip("N/A. ")
            pssm_length = obj[4]

            if re.match(r"cl\d+", accession, re.I):
                sets[accession] = {
                    "accession": accession,
                    "name": short_name,
                    "description": descr
                }
                nodes[accession] = set()
            elif re.match(r"cd\d+", accession, re.I):
                domains[accession] = {
                    "name": short_name,
                    "description": descr
                }

    os.remove(filename)

    filename, headers = urllib.request.urlretrieve(_LINKS)
    with open(filename, "rt") as fh:
        for line in fh:
            cd_id, cd_pssm_id, cl_id, cl_pssm_id = line.rstrip().split("\t")

            if re.match(r"cl\d+", cl_id, re.I) and re.match(r"cd\d+", cd_id, re.I):
                nodes[cl_id].add(cd_id)

    os.remove(filename)

    filename, headers = urllib.request.urlretrieve(_CDTRACK)
    with open(filename, "rt") as fh:
        for line in fh:
            line = line.rstrip()

            if not line or line[0] == "#":
                continue

            fields = line.split()
            ac = fields[0]
            parent_ac = fields[3]

            if parent_ac in domains:
                parent_of[ac] = parent_ac

    os.remove(filename)

    final_sets = {}
    for cl_id in sets:
        _nodes = []
        _links = []
        for cd_id in nodes[cl_id]:
            _nodes.append({
                "accession": cd_id,
                "type": "entry",
                "score": 1
            })

            if cd_id in parent_of:
                _links.append({
                    "source": parent_of[cd_id],
                    "target": cd_id,
                    "score": 1
                })

        if _nodes:
            sets[cl_id]["relationships"] = {
                "nodes": _nodes,
                "links": _links
            }

            final_sets[cl_id] = sets[cl_id]

    return final_sets
