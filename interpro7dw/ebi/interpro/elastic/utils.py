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


def find_files(root: str, wait_for_sentinel: bool=False):
    sentinel = os.path.join(root, LOADING)

    while wait_for_sentinel and not os.path.isfile(sentinel):
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

    while True:
        """
        On the first pass:
            * we read from `indir` (where documents files are created)
            * if `outdir` is passed, we set `wait_for_sentinel` to True, 
              i.e. we do not index documents until the sentinel file exists. 
            * if `writeback` is True, documents files are updated with
              the documents that failed to be indexed. If all documents were
              indexed, the file is deleted. 
            * if `writeback` is False and `outdir` is passed, 
              documents that failed to be indexed are written in an output
              directory
              
        On the subsequent passes:
            * we read from `outdir` if passed, from `indir` otherwise
            * we set `wait_for_sentinel` to False 
              (if we were expecting a sentinel file, 
              it must have been created during the first pass)
            * if `outdir` is passed, we set `writeback` to True, 
              so files in `outdir` are updated/deleted

        /!\ 
        If `outdir` is not passed, and `writeback` is False, we will loop
        until all documents in all files are successfully indexed in bulk:
            * `writeback` is False: we do not update files
            * `outdir` is not set: we do not save failed documents       
        """
        if first_pass:
            root = indir
            wait_for_sentinel = outdir is not None
        else:
            if outdir:
                root = outdir
                writeback = True
            else:
                root = indir

            wait_for_sentinel = False

        for filepath in find_files(root, wait_for_sentinel):
            docs = dataload(filepath)

            if first_pass:
                # Count the number of documents to index only once
                num_documents += len(docs)

            if callback:
                actions = map(callback, docs)
            else:
                actions = docs

            failed = []
            for i, (ok, info) in enumerate(pbulk(es, actions, **kwargs)):
                if ok:
                    num_indexed += 1
                else:
                    failed.append(docs[i])

            if failed:
                if writeback:
                    # Overwrite file with failed documents
                    datadump(filepath, failed)
                elif organizer:
                    # Write failed documents to new file
                    filepath = organizer.mktemp(suffix=EXTENSION)
                    datadump(filepath, failed)
            elif writeback:
                # Remove file as all documents have been successfully indexed
                os.remove(filepath)

            logger.info(f"{num_indexed:>12,} / {num_documents:,}")

        first_pass = False

        if num_indexed == num_documents:
            break
