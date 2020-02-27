# -*- coding: utf-8 -*-

import glob
import logging
import os
import time
import shutil
from typing import Callable, Optional, Sequence

from elasticsearch import Elasticsearch, exceptions
from elasticsearch.helpers import parallel_bulk

from interpro7dw import logger
from interpro7dw.utils import DirectoryTree, dataload, datadump


DEFAULT_SHARDS = 5
EXTENSION = ".dat"
LOADING = "loading"


def _delete_index(es: Elasticsearch, name: str):
    # Make sure the index is deleted
    while True:
        try:
            es.indices.delete(index=name)
        except exceptions.NotFoundError:
            break
        except Exception as exc:
            time.sleep(10)
        else:
            break


def add_alias(es: Elasticsearch, indices: Sequence[str], name: str,
              delete_indices: bool):
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

        if delete_indices:
            # Delete indices previously pointed by alias
            for index in current_indices:
                _delete_index(es, index)
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
    _delete_index(es, name)

    # Then make sure it's created
    while True:
        try:
            es.indices.create(index=name, body=body)
        except exceptions.RequestError as exc:
            raise exc
        except Exception as exc:
            time.sleep(10)
        else:
            break


def find_files(root: str, wait_until_init: bool=False):
    sentinel = os.path.join(root, LOADING)

    while wait_until_init and not os.path.isfile(sentinel):
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

        if os.path.isfile(sentinel):
            # Sentinel exists (files are still written): keep going
            time.sleep(60)
        elif active:
            # All files ready: they will all be found in the next iteration
            active = False
        else:
            break


def index_documents(es: Elasticsearch, indir: str,
                    callback: Optional[Callable]=None,
                    outdir: Optional[str]=None, threads: int=4,
                    writeback: bool=False):
    logger.info("starting")

    num_documents = 0
    num_indexed = 0
    num_iter = 0

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

    while True:
        if num_iter and outdir:
            root = outdir
            wait = False
            _writeback = True
        else:
            root = indir
            wait = True  # we expect a sentinel to be in indir
            _writeback = writeback

        for filepath in find_files(root, wait_until_init=wait):
            docs = dataload(filepath)

            if not num_iter:
                """
                Count the number of documents to index only during
                the first iteration as subsequent iterations will use
                documents that already have been tried to be indexed
                """
                num_documents += len(docs)

            if callback:
                docs = list(map(callback, docs))

            failed = []
            for i, (ok, info) in enumerate(parallel_bulk(es, docs, **kwargs)):
                if ok:
                    num_indexed += 1
                else:
                    failed.append(docs[i])

            if failed:
                if _writeback:
                    # Overwrite file with failed documents
                    datadump(filepath, failed)
                elif organizer:
                    # Write failed documents to new file
                    filepath = organizer.mktemp(suffix=EXTENSION)
                    datadump(filepath, failed)
            elif _writeback:
                # Remove file as all documents have been successfully indexed
                os.remove(filepath)

            logger.info(f"{num_indexed:>12,} / {num_documents:,}")

        num_iter += 1

        if num_indexed == num_documents:
            break

