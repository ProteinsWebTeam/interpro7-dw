#!/usr/bin/env python
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
            (pssm_id, accession, short_name,
             descr, pssm_length) = line.rstrip().split("\t")
            descr = descr.lstrip("N/A. ")

            if re.match(r"cl\d+", accession):
                sets[accession] = {
                    "accession": accession,
                    "name": short_name,
                    "description": descr
                }
                nodes[accession] = set()
            elif re.match(r"cd\d+", accession):
                domains[accession] = {
                    "name": short_name,
                    "description": descr
                }

    os.remove(filename)

    filename, headers = urllib.request.urlretrieve(_LINKS)
    with open(filename, "rt") as fh:
        for line in fh:
            cd_id, cd_pssm_id, cl_id, cl_pssm_id = line.rstrip().split("\t")
            if re.match(r"cl\d+", cl_id) and re.match(r"cd\d+", cd_id):
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
            cd_id = cd_id
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
