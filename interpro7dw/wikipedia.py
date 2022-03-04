import json
import time
from base64 import b64encode
from typing import Callable, Optional
from http.client import IncompleteRead
from urllib.parse import quote, unquote, urlencode
from urllib.error import HTTPError
from urllib.request import Request, urlopen
import xml.etree.ElementTree as ET


API = "https://en.wikipedia.org/w/api.php"
REST_API = "https://en.wikipedia.org/api/rest_v1/page/summary/"

# See https://meta.wikimedia.org/wiki/User-Agent_policy
USER_AGENT = "InterPro/1.0 (https://www.ebi.ac.uk/interpro/) Python-urllib/3"


def get_ext_links(query: str, validate: Optional[Callable] = None) -> set[str]:
    # See https://www.mediawiki.org/wiki/API:Exturlusage
    params = {
        "action": "query",
        "format": "json",
        "list": "exturlusage",
        "euquery": query,
        "eulimit": 100,
        # "eunamespace": 0  # checked
    }

    pages = set()
    while True:
        url = f"{API}?{urlencode(params)}"
        with urlopen(url) as f:
            data = json.loads(f.read().decode("utf-8"))

            for ext_obj in data["query"]["exturlusage"]:
                if ext_obj["ns"] == 0:
                    # Page is an article
                    page_title = ext_obj["title"]
                    ext_url = ext_obj['url']

                    if validate is None or validate(ext_url):
                        pages.add(page_title)

        try:
            params.update(data["continue"])
        except KeyError:
            break

    return pages


def get_summary(title: str, max_retries: int = 4):
    # Some records contains HTML %xx escapes: we need to replace them
    title = unquote(title)

    # default `safe` is '/' but we *want* to replace it
    url = REST_API + quote(title, safe='')

    obj = None
    num_retries = 0
    req = Request(url)
    req.add_header("User-Agent", USER_AGENT)
    while True:
        try:
            res = urlopen(req)
            data = res.read()
        except HTTPError as e:
            # Content can still be retrieved with e.fp.read()
            break
        except IncompleteRead:
            if num_retries == max_retries:
                break
            else:
                num_retries += 1
                time.sleep(3)
        else:
            obj = json.loads(data.decode("utf-8"))
            break

    return obj


def get_thumbnail(summary: dict, max_retries: int = 4):
    try:
        thumbnail = summary["thumbnail"]
    except KeyError:
        return None

    num_retries = 0
    b64string = None
    req = Request(thumbnail["source"])
    req.add_header("User-Agent", USER_AGENT)
    while True:
        try:
            res = urlopen(req)
            data = res.read()
        except HTTPError as e:
            break
        except IncompleteRead:
            if num_retries == max_retries:
                break
            else:
                num_retries += 1
                time.sleep(3)
        else:
            b64string = b64encode(data).decode("utf-8")
            break

    return b64string


def parse_infobox(title: str,
                  validate: Optional[Callable] = None) -> dict[str, set[str]]:
    params = {
        "action": "parse",
        "page": title,
        "format": "json",
        "prop": "parsetree",
    }

    url = f"{API}?{urlencode(params)}"
    with urlopen(url) as f:
        data = json.loads(f.read().decode("utf-8"))

    root = ET.fromstring(data["parse"]["parsetree"]["*"])
    props = {}
    for template in root.findall("template"):
        for part in template.findall("part"):
            try:
                name = part.find("name").text.strip().lower()
                value = part.find("value").text.strip()
            except AttributeError:
                continue

            if name and value:
                if validate is None or validate(name, value):
                    if name in props:
                        props[name].add(value)
                    else:
                        props[name] = {value}

    return props
