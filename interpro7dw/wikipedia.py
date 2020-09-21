import json
import time
from base64 import b64encode
from http.client import IncompleteRead
from urllib.parse import quote, unquote
from urllib.error import HTTPError
from urllib.request import urlopen

from interpro7dw import logger


URL = "https://en.wikipedia.org/api/rest_v1/page/summary/"


def get_summary(title: str, max_retries: int = 4):
    # Some records contains HTML %xx escapes: we need to replace them
    title = unquote(title)

    # default `safe` is '/' but we *want* to replace it
    url = URL + quote(title, safe='')

    obj = None
    num_retries = 0
    while True:
        try:
            res = urlopen(url)
            data = res.read()
        except HTTPError as e:
            # Content can still be retrieved with e.fp.read()
            logger.error(f"{title}: {e.code} ({e.reason})")
            break
        except IncompleteRead:
            if num_retries == max_retries:
                logger.error(f"{title}: incomplete")
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
    while True:
        try:
            res = urlopen(thumbnail["source"])
            data = res.read()
        except HTTPError as e:
            logger.error(f"{summary['title']} (thumbnail): "
                         f"{e.code} ({e.reason})")
            break
        except IncompleteRead:
            if num_retries == max_retries:
                logger.error(f"{summary['title']} (thumbnail): incomplete")
                break
            else:
                num_retries += 1
                time.sleep(3)
        else:
            b64string = b64encode(data).decode("utf-8")
            break

    return b64string
