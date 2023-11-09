import os
import re
import shutil
import tarfile
from tempfile import mkdtemp, mkstemp
from urllib.request import HTTPBasicAuthHandler, build_opener, urlopen
from urllib.error import HTTPError

from interpro7dw.utils import logger


def load_enzyme_xrefs(filepath) -> dict[str, set[str]]:
    pathways = {}
    with open(filepath, "rt", encoding="latin-1") as fh:
        reaction_pathways = []
        ecno = None

        for line in fh:
            if line[0] == '#':
                continue

            line = line.rstrip()
            m = re.match(r"EC-NUMBER - EC-(\d+\.\d+\.\d+(\.\d+)?)", line)
            if m:
                g1, g2 = m.groups()
                if g2 is None:
                    # 3.4.19 -> 3.4.19.-
                    ecno = g1 + ".-"
                else:
                    ecno = g1

            m = re.match(r"IN-PATHWAY - (PWYG?-\d+)", line)
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


def load_pathways(filepath) -> dict[str, str]:
    pathways = {}
    with open(filepath, "rt", encoding="latin-1") as fh:
        pathway_id = None
        pathway_name = None
        for line in fh:
            if line[0] == '#':
                continue

            line = line.rstrip()
            m = re.match(r"UNIQUE-ID - (PWYG?-\d+)", line)
            if m:
                pathway_id = m.group(1)

            m = re.match(r"COMMON-NAME - (.+)", line)
            if m:
                pathway_name = m.group(1)

                # Replacing HTML symbols for greek letters, e.g. &alpha; -> alpha
                pathway_name = re.sub(r"&([a-z]+);", r"\1", pathway_name, flags=re.I)

                # Remove HTML tags (e.g. <em>c</em> -> c)
                pathway_name = re.sub(r"</?.+?>", r'', pathway_name)

            if line == "//":
                if pathway_id:
                    pathways[pathway_id] = pathway_name

                pathway_id = None
                pathway_name = None

    return pathways


def download_archive(username: str, password: str) -> str:
    url = "http://brg-files.ai.sri.com/public/dist/meta.tar.gz"

    try:
        res = urlopen(url)
    except HTTPError as err:
        res = err

    if res.getcode() == 401:
        """
        Retry with basic authentication
        Find realm in headers:
            e.g. Www-authenticate: Basic realm="SRI BRG Restricted access"
        """
        header = res.info()["WWW-Authenticate"]
        realm = re.match("Basic realm=[\"'](.+?)[\"']", header).group(1)

        auth_handler = HTTPBasicAuthHandler()
        auth_handler.add_password(realm=realm,
                                  uri=url,
                                  user=username,
                                  passwd=password)

        opener = build_opener(auth_handler)

        try:
            res = opener.open(url)
        except HTTPError as err:
            res = err

    if res.getcode() != 200:
        raise HTTPError(res.geturl(), res.getcode(), '', res.info(), None)

    fd, filepath = mkstemp()
    os.close(fd)

    with open(filepath, "wb") as fh:
        for chunk in res:
            fh.write(chunk)

    return filepath


def get_ec2pathways(filepath: str) -> dict[str, list[tuple]]:
    path = mkdtemp()
    pathways = None
    xrefs = None
    with tarfile.open(filepath, "r") as tar:
        for name in tar.getnames():
            if "pathways.dat" in name:
                tar.extract(name, path=path, set_attrs=False)
                pathways = load_pathways(os.path.join(path, name))
            elif "reactions.dat" in name:
                tar.extract(name, path=path, set_attrs=False)
                xrefs = load_enzyme_xrefs(os.path.join(path, name))

    shutil.rmtree(path)

    if pathways is None:
        raise RuntimeError("'pathways.dat' not in the archive")
    elif xrefs is None:
        raise RuntimeError("'reactions.dat' not in the archive")

    ec2pathways = {}
    for pathway_id, enzymes in xrefs.items():
        try:
            name = pathways[pathway_id]
        except KeyError:
            logger.error(f"missing pathway {pathway_id}")
            continue

        for ecno in enzymes:
            try:
                ec2pathways[ecno].append((pathway_id, name))
            except KeyError:
                ec2pathways[ecno] = [(pathway_id, name)]

    return ec2pathways
