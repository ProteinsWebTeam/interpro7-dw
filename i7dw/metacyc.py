# -*- coding: utf-8 -*-

import os
import re
import shutil
import tarfile
from tempfile import mkdtemp, mkstemp
from typing import Dict, List
from urllib.request import urlopen
from urllib.error import HTTPError


def load_enzyme_xrefs(filepath):
    pathways = {}
    with open(filepath, "rt", encoding="latin-1") as fh:
        reaction_pathways = []
        ecno = None

        for line in fh:
            if line[0] == '#':
                continue

            line = line.rstrip()
            m = re.match("EC-NUMBER - EC-(\d+\.\d+\.\d+(\.\d+)?)", line)
            if m:
                g1, g2 = m.groups()
                if g2 is None:
                    # 3.4.19 -> 3.4.19.-
                    ecno = g1 + ".-"
                else:
                    ecno = g1

            m = re.match("IN-PATHWAY - (PWYG?-\d+)", line)
            if m:
                reaction_pathways.append(m.group(1))

            if line == "//":
                if ecno:
                    for pathway_id in reaction_pathways:
                        try:
                            pathways[pathway_id].add(ecno)
                        except KeyError:
                            pathways[pathway_id] = {ecno}

                reaction_pathways = []
                ecno = None

    return pathways


def load_pathways(filepath):
    pathways = {}
    with open(filepath, "rt", encoding="latin-1") as fh:
        pathway_id = None
        pathway_name = None
        for line in fh:
            if line[0] == '#':
                continue

            line = line.rstrip()
            m = re.match("UNIQUE-ID - (PWYG?-\d+)", line)
            if m:
                pathway_id = m.group(1)

            m = re.match("COMMON-NAME - (.+)", line)
            if m:
                pathway_name = m.group(1)

                # Replacing HTML symbols for greek letters, e.g. &alpha; -> alpha
                pathway_name = re.sub(r"&([a-z]+);", r"\1", pathway_name, flags=re.I)

                # Remove HTML tags (e.g. <em>c</em> -> c)
                pathway_name = re.sub(r"</?.+?>", '', pathway_name)

            if line == "//":
                if pathway_id:
                    pathways[pathway_id] = pathway_name

                pathway_id = None
                pathway_name = None

    return pathways


def get_ec2pathways() -> Dict[str, List[tuple]]:
    url = "http://brg-files.ai.sri.com/public/dist/meta.tar.gz"
    fd, filepath = mkstemp()
    os.close(fd)

    with open(filepath, "wb") as fh, urlopen(url) as res:
        if res.getcode() != 200:
            raise HTTPError(res.geturl(), res.getcode(), '', res.info(), None)

        for chunk in res:
            fh.write(chunk)

    path = mkdtemp()
    pathways = None
    xrefs = None
    with tarfile.open(filepath, "r") as tar:
        for name in tar.getnames():
            if name.endswith("pathways.dat"):
                tar.extract(name, path=path, set_attrs=False)
                pathways = load_pathways(os.path.join(path, name))
            elif name.endswith("reactions.dat"):
                tar.extract(name, path=path, set_attrs=False)
                xrefs = load_enzyme_xrefs(os.path.join(path, name))

    os.remove(filepath)
    shutil.rmtree(path)

    if pathways is None:
        raise RuntimeError("'pathways.dat' not in the archive")
    elif xrefs is None:
        raise RuntimeError("'reactions.dat' not in the archive")

    ec2pathways = {}
    for pathway_id in xrefs:
        name = pathways[pathway_id]

        for ecno in xrefs[pathway_id]:
            try:
                ec2pathways[ecno].append((pathway_id, name, "metacyc"))
            except KeyError:
                ec2pathways[ecno] = [(pathway_id, name, "metacyc")]

    return ec2pathways
