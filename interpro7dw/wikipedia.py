import json
import time
from base64 import b64encode
from http.client import IncompleteRead
from urllib.parse import quote, unquote
from urllib.error import HTTPError
from urllib.request import Request, urlopen


URL = "https://en.wikipedia.org/api/rest_v1/page/summary/"

# See https://meta.wikimedia.org/wiki/User-Agent_policy
USER_AGENT = "InterPro/1.0 (https://www.ebi.ac.uk/interpro/) Python-urllib/3"


def get_summary(title: str, max_retries: int = 4):
    # Some records contains HTML %xx escapes: we need to replace them
    title = unquote(title)

    # default `safe` is '/' but we *want* to replace it
    url = URL + quote(title, safe='')

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
