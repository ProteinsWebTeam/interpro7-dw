# -*- coding: utf-8 -*-

import re
from typing import Dict, List
from urllib.request import urlopen


def fetch_enzyme_xrefs():
    url = "http://rest.kegg.jp/link/enzyme/pathway"
    pathways = {}
    with urlopen(url) as res:
        lines = res.read().decode("utf-8").strip().split('\n')
        for line in lines:
            m = re.match("path:map(\d+)\s+(ec:\d+\.\d+\.\d+\.\d+)$", line)
            if not m:
                continue

            pathway_id = m.group(1)
            ecno = m.group(2)

            try:
                pathways[pathway_id].add(ecno)
            except KeyError:
                pathways[pathway_id] = {ecno}

    return pathways


def fetch_pathways():
    url = "http://rest.kegg.jp/list/pathway"
    pathways = {}
    with urlopen(url) as res:
        lines = res.read().decode("utf-8").strip().split('\n')
        for line in lines:
            m = re.match("path:map(\d+)\s+(.*)$", line)
            pathway_id = m.group(1)
            name = m.group(2)
            pathways[pathway_id] = name

    return pathways


def get_ec2pathways() -> Dict[str, List[tuple]]:
    pathways = fetch_pathways()
    xrefs = fetch_enzyme_xrefs()

    ec2pathways = {}
    for pathway_id in xrefs:
        name = pathways[pathway_id]

        for ecno in xrefs[pathway_id]:
            try:
                ec2pathways[ecno].append((pathway_id, name))
            except KeyError:
                ec2pathways[ecno] = [(pathway_id, name)]

    return ec2pathways
