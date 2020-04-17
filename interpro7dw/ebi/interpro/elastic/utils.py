# -*- coding: utf-8 -*-

import glob
import logging
import os
import time
import shutil
from typing import Callable, Optional, Sequence

from elasticsearch import Elasticsearch, exceptions
from elasticsearch.helpers import parallel_bulk as pbulk

from interpro7dw import logger
from interpro7dw.utils import DirectoryTree, loadobj, dumpobj


DEFAULT_SHARDS = 5
EXTENSION = ".dat"
LOAD_SUFFIX = ".load"
DONE_SUFFIX = ".done"


def delete_index(es: Elasticsearch, name: str):
    # Make sure the index is deleted
    while True:
        try:
            es.indices.delete(index=name)
        except exceptions.ConnectionTimeout:
            time.sleep(30)
        except exceptions.NotFoundError:
            break
        except Exception as exc:
            logger.error(f"{name}: {exc}")
            time.sleep(30)
        else:
            break


def add_alias(es: Elasticsearch, indices: Sequence[str], name: str):
    if es.indices.exists_alias(name=name):
        # Alias already exists: update it

        # Current indices pointed by alias
        current_indices = set(es.indices.get_alias(name=name))

        actions = []
        for index in set(indices):
            if index in current_indices:
                # index is already pointed by alias
                current_indices.remove(index)
                continue

            # Point alias to index
            actions.append({"add": {"index": index, "alias": name}})

        # Remove alias from old indices
        for index in current_indices:
            actions.append({"remove": {"index": index, "alias": name}})

        if actions:
            """
            Atomic operation:
            Alias removed from the old indices
            at the same time it's added to the new ones
            """
            es.indices.update_aliases(body={"actions": actions})
    else:
        # Creat new alias, then point to indices
        es.indices.put_alias(index=','.join(indices), name=name)

    # Update index settings
    for index in indices:
        es.indices.put_settings(
            body={
                # "number_of_replicas": 1,
                "refresh_interval": None  # default (1s)
            },
            index=index
        )


def connect(hosts: Sequence[str], verbose: bool=True) -> Elasticsearch:
    es = Elasticsearch(hosts=hosts)

    if not verbose:
        # Disable Elastic logger
        tracer = logging.getLogger("elasticsearch")
        tracer.setLevel(logging.CRITICAL + 1)

    return es


def create_index(es: Elasticsearch, name: str, body: dict):
    delete_index(es, name)

    # Then make sure it's created
    while True:
        try:
            es.indices.create(index=name, body=body)
        except exceptions.ConnectionTimeout:
            time.sleep(30)
        except exceptions.RequestError as exc:
            raise exc
        except Exception as exc:
            logger.warning(f"{name}: {exc}")
            time.sleep(30)
        else:
            break


def find_files(root: str, version: str):
    if version:
        load_sentinel = os.path.join(root, f"{version}{LOAD_SUFFIX}")
        done_sentinel = os.path.join(root, f"{version}{DONE_SUFFIX}")

        if not os.path.isfile(done_sentinel):
            while not os.path.isfile(load_sentinel):
                time.sleep(60)

    pathname = os.path.join(root, "**", f"*{EXTENSION}")

    files = set()
    active = True
    while True:
        for path in glob.iglob(pathname, recursive=True):
            if path in files:
                continue

            files.add(path)
            yield path

        if not version or not active:
            break
        elif os.path.isfile(done_sentinel):
            # All files ready: they will all be found in the next iteration
            active = False
        else:
            # Files are still being written
            time.sleep(60)


def index_documents(es: Elasticsearch, indir: str, version: str,
                    callback: Optional[Callable]=None,
                    outdir: Optional[str]=None, threads: int=4,
                    writeback: bool=False):
    logger.info("starting")

    num_documents = 0
    num_indexed = 0
    first_pass = True

    if outdir:
        try:
            # We don't want to use old documents: ensure the directory is deleted
            shutil.rmtree(outdir)
        except FileNotFoundError:
            pass

        organizer = DirectoryTree(outdir)
    else:
        organizer = None

    kwargs = {
        "thread_count": threads,
        "queue_size": threads,
        "raise_on_exception": False,
        "raise_on_error": False
    }

    ts_then = time.time()
    while True:
        """
        On the first pass:
            * we read from `indir` (where documents files are created)
            * if `writeback` is False (default), we use `version` to know
              when it's OK to start indexing, and when all files were written.
            * if `writeback` is False and `outdir` is passed,
              failed documents are written in an output directory.
            * if `writeback` is True, `version` is ignored and files in `indir`
              are overwritten with failed documents. 
              If all documents were indexed, the file is deleted.
              
        On the subsequent passes:
            * we read from `outdir` if passed, from `indir` otherwise
            * if `outdir` is passed, we set `writeback` to True
              so files in `outdir` are updated/deleted,
              and `version` is ignored.

        /!\ 
        If `outdir` is not passed, and `writeback` is False, we will loop
        until all documents in all files are successfully indexed in bulk:
            * `writeback` is False: we do not update files
            * `outdir` is not set: we do not save failed documents       
        """
        if first_pass:
            root = indir
            if writeback:
                version = None
        else:
            if outdir:
                root = outdir
                writeback = True
                version = None
            else:
                root = indir

        for filepath in find_files(root, version):
            docs = loadobj(filepath)

            if first_pass:
                # Count the number of documents to index only once
                num_documents += len(docs)

            if callback:
                actions = map(callback, docs)
            else:
                actions = docs

            failed = []
            pause = False
            for i, (ok, info) in enumerate(pbulk(es, actions, **kwargs)):
                if ok:
                    num_indexed += 1
                    continue

                failed.append(docs[i])

                try:
                    is_429 = info["index"]["status"] == 429
                except (KeyError, IndexError):
                    is_429 = False

                try:
                    exc = info["index"]["exception"]
                except (KeyError, TypeError):
                    exc = None

                if is_429 or isinstance(exc, exceptions.ConnectionTimeout):
                    pause = True
                else:
                    logger.debug(info)

            if failed:
                if writeback:
                    # Overwrite file with failed documents
                    dumpobj(filepath, failed)
                elif organizer:
                    # Write failed documents to new file
                    filepath = organizer.mktemp(suffix=EXTENSION)
                    dumpobj(filepath, failed)
            elif writeback:
                # Remove file as all documents have been successfully indexed
                os.remove(filepath)

            ts_now = time.time()
            if ts_now - ts_then >= 3600:
                logger.info(f"{num_indexed:>14,} / {num_documents:,}")
                ts_then = ts_now

            if pause:
                time.sleep(30)

        logger.info(f"{num_indexed:>14,} / {num_documents:,}")
        first_pass = False

        if num_indexed == num_documents:
            break


def publish(hosts: Sequence[str], staging: str, live: str, previous: str):
    es = connect(hosts, verbose=False)

    if es.indices.exists_alias(name=live):
        # Make 'live' indices also 'previous'
        indices = es.indices.get_alias(name=live)
        add_alias(es, indices, previous)

    # Make 'staging' indices public ('live')
    indices = es.indices.get_alias(name=staging)
    add_alias(es, indices, live)
