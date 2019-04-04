import glob
import json
import logging
import os
import shutil
import time
from multiprocessing import Process, Queue
from typing import List, Optional

from elasticsearch import Elasticsearch, helpers, exceptions

from .organize import JsonFileOrganizer, is_ready
from .. import mysql
from ... import logger

NODB_INDEX = "others"
SRCH_INDEX = "iprsearch"


class DocumentLoader(Process):
    def __init__(self, hosts: List, doc_type: str, task_queue: Queue,
                 done_queue: Queue, **kwargs):
        super().__init__()
        self.hosts = hosts
        self.type = doc_type
        self.task_queue = task_queue
        self.done_queue = done_queue
        self.chunk_size = kwargs.get("chunk_size", 500)
        self.max_bytes = kwargs.get("max_bytes", 100 * 1024 * 1024)
        self.suffix = kwargs.get("suffix", "")
        self.threads = kwargs.get("threads", 4)
        self.err_log = kwargs.get("err_log")

    def run(self):
        es = Elasticsearch(self.hosts, timeout=30)

        # Disable Elastic logger
        tracer = logging.getLogger("elasticsearch")
        tracer.setLevel(logging.CRITICAL + 1)

        if self.err_log:
            err_dst = self.err_log
        else:
            err_dst = os.devnull

        with open(err_dst, "wt") as err:
            for filepath in iter(self.task_queue.get, None):
                with open(filepath, "rt") as fh:
                    documents = json.load(fh)

                actions = {}
                for doc in documents:
                    # Define which index to use
                    try:
                        index = doc["entry_db"] + self.suffix
                    except KeyError:
                        index = SRCH_INDEX + self.suffix
                    except TypeError:
                        """
                        `entry_db` is None
                            -> use the index for documents without entry
                        """
                        index = NODB_INDEX + self.suffix

                    _id = doc["id"]
                    actions[_id] = {
                        "_op_type": "index",
                        "_index": index,
                        "_type": self.type,
                        "_id": _id,
                        "_source": doc
                    }

                gen = helpers.parallel_bulk(
                    es, actions.values(),
                    thread_count=self.threads,
                    queue_size=self.threads,
                    chunk_size=self.chunk_size,
                    max_chunk_bytes=self.max_bytes,
                    raise_on_exception=False,
                    raise_on_error=False
                )

                num_successful = 0
                failed = []
                for status, item in gen:
                    if status:
                        num_successful += 1
                    else:
                        err.write("{}\n".format(item))

                        """
                        Some items have a `data` key under `index`...
                        but some don't (so a KeyError would be raised)
                        """
                        _id = item["index"]["_id"]
                        failed.append(actions[_id]["_source"])

                self.done_queue.put((filepath, num_successful, failed))


def create_indices(my_ippro: str, hosts: List[str], body_f: str,
                   shards: Optional[str], num_shards: int=5, suffix: str=""):

    logger.info("creating indices in: {}".format(", ".join(hosts)))

    # Load default settings and property mapping
    with open(body_f, "rt") as fh:
        body = json.load(fh)

    # Pop mappings (only one mapping type per index)
    mappings = body.pop("mappings")

    if shards:
        # Custom number of shards
        with open(shards, "rt") as fh:
            shards = json.load(fh)

        # Force keys to be in lower case
        shards = {k.lower(): v for k, v in shards.items()}
    else:
        shards = {}

    # Load databases (base of indices)
    databases = list(mysql.database.get_databases(my_ippro).keys())

    # Establish connection
    es = Elasticsearch(hosts, timeout=30)

    # Disable Elastic logger
    tracer = logging.getLogger("elasticsearch")
    tracer.setLevel(logging.CRITICAL + 1)

    # Create indices
    for index in databases + [NODB_INDEX, SRCH_INDEX]:
        try:
            n_shards = shards[index]
        except KeyError:
            # No custom number of shards for this index
            n_shards = num_shards

        # Update body with mappings to use for current index
        _type = "search" if index == SRCH_INDEX else "relationship"
        body["mappings"] = {_type: mappings[_type]}

        """
        Change settings for large bulk imports:

        https://www.elastic.co/guide/en/elasticsearch/reference/current/indices-update-settings.html
        https://www.elastic.co/guide/en/elasticsearch/guide/current/indexing-performance.html
        """
        try:
            settings = body["settings"]
        except KeyError:
            settings = body["settings"] = {}
        finally:
            settings.update({
                "number_of_shards": n_shards,
                "number_of_replicas": 0,    # default: 1
                "refresh_interval": -1,     # default: 1s
                "codec": "best_compression"
            })

        index += suffix

        # Make sure the index is deleted
        while True:
            try:
                res = es.indices.delete(index)
            except exceptions.NotFoundError:
                break
            except Exception as exc:
                logger.error("{}: {}".format(type(exc), exc))
                time.sleep(10)
            else:
                break

        # And make sure it's created
        while True:
            try:
                res = es.indices.create(index, body=body)
            except exceptions.RequestError as exc:
                break  # raised if index exists
            except Exception as exc:
                logger.error("{}: {}".format(type(exc), exc))
                time.sleep(10)
            else:
                break

    logger.info("complete")


def find_json_files(src: str, seconds: int=60):
    pathname = os.path.join(src, "**", "*.json")
    files = set()
    stop = False
    while True:
        for path in glob.iglob(pathname, recursive=True):
            if path not in files:
                files.add(path)
                yield path

        if stop:
            break
        elif is_ready(src):
            # All files ready, but loop one last time
            stop = True
        else:
            time.sleep(seconds)


def index_documents(my_ippro: str, hosts: List[str], doc_type: str,
                    src: str, **kwargs):
    alias = kwargs.get("alias")
    err_prefix = kwargs.get("err_prefix")
    max_retries = kwargs.get("max_retries", 0)
    processes = kwargs.get("processes", 1)
    raise_on_error = kwargs.get("raise_on_error", True)
    suffix = kwargs.get("suffix", "").lower()
    failed_docs_dir = kwargs.get("failed_docs_dir")
    write_back = kwargs.get("write_back", False)

    if failed_docs_dir:
        try:
            """
            Ensure the directory does not exist
            as we don't want files from a previous run to be considered
            """
            shutil.rmtree(failed_docs_dir)
        except FileNotFoundError:
            pass
        finally:
            os.makedirs(failed_docs_dir)
            organizer = JsonFileOrganizer(failed_docs_dir)
    else:
        organizer = None

    processes = max(1, processes - 1)  # consider parent process
    num_retries = 0
    while True:
        logger.info("indexing documents (try #{})".format(num_retries + 1))
        task_queue = Queue()
        done_queue = Queue(maxsize=processes)
        workers = []

        for i in range(processes):
            if err_prefix:
                err_log = err_prefix + str(i + 1) + ".err"
            else:
                err_log = None

            w = DocumentLoader(hosts, doc_type, task_queue, done_queue,
                               suffix=suffix, err_log=err_log)
            w.start()
            workers.append(w)

        num_files = 0
        for path in find_json_files(src):
            task_queue.put(path)
            num_files += 1

        for _ in workers:
            task_queue.put(None)

        tot_successful = 0
        tot_failed = 0
        for i in range(num_files):
            path, num_successful, failed = done_queue.get()
            tot_successful += num_successful
            tot_failed += len(failed)

            if failed:
                if organizer:
                    for doc in failed:
                        organizer.add(doc)
                    organizer.flush()
                elif write_back:
                    with open(path, "wt") as fh:
                        json.dump(failed, fh)
            elif write_back:
                os.remove(path)

            if not (i + 1 ) % 1000:
                logger.info("documents indexed: {:>15,} "
                            "({:,} failed)".format(tot_successful,
                                                   tot_failed))

        logger.info("documents indexed: {:>15,} "
                    "({:,} failed)".format(tot_successful, tot_failed))

        # Wait for workers to complete
        for w in workers:
            w.join()

        if not tot_failed or num_retries == max_retries:
            break
        else:
            num_retries += 1

    if tot_failed and raise_on_error:
        raise RuntimeError("{:,} documents not indexed".format(tot_failed))
    elif not tot_failed and alias:
        """
        We DO NOT want to delete old indices as they are used in production
        They will be deleted when we switch transparently between them
            and the new indices
        """
        update_alias(my_ippro, hosts, alias, suffix=suffix, delete=False)

    logger.info("complete")


def update_alias(my_ippro: str, hosts: List, alias: str, **kwargs):
    suffix = kwargs.get("suffix", "").lower()
    delete = kwargs.get("delete", False)

    databases = list(mysql.database.get_databases(my_ippro).keys())
    new_indices = set()
    for index in databases + [NODB_INDEX, SRCH_INDEX]:
        new_indices.add(index + suffix)

    es = Elasticsearch(hosts, timeout=30)

    # Disable Elastic logger
    tracer = logging.getLogger("elasticsearch")
    tracer.setLevel(logging.CRITICAL + 1)

    exists = es.indices.exists_alias(name=alias)
    if exists:
        # Alias already exists: update it

        # Indices currently using the alias
        indices = set(es.indices.get_alias(name=alias))

        actions = []
        for index in new_indices:
            try:
                # If passes: new index is already using the alias
                indices.remove(index)
            except KeyError:
                # Otherwise, add the alias to the new index
                actions.append({
                    "add": {
                        "index": index,
                        "alias": alias
                    }
                })

        # Remove the alias from the old indices
        for index in indices:
            actions.append({
                "remove": {
                    "index": index,
                    "alias": alias
                }
            })

        if actions:
            """
            Atomic operation:
            Alias removed from the old indices
            at the same time it's added to the new ones
            """
            es.indices.update_aliases(body={"actions": actions})

        if delete:
            # Delete old indices that used the alias
            for index in indices:
                while True:
                    try:
                        res = es.indices.delete(index)
                    except exceptions.NotFoundError:
                        break
                    except Exception as exc:
                        logger.error("{}: {}".format(type(exc), exc))
                    else:
                        break
    else:
        # Create alias, and add it to new indices
        es.indices.put_alias(index=','.join(new_indices), name=alias)

    # Update index settings
    for index in new_indices:
        es.indices.put_settings({
            # "number_of_replicas": 1,
            "refresh_interval": None  # default (1s)
        }, index)
