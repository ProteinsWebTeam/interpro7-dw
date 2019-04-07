import glob
import json
import logging
import os
import shutil
import time
from abc import ABC, abstractmethod
from multiprocessing import Process, Queue
from typing import List, Optional

from elasticsearch import Elasticsearch, helpers, exceptions

from .organize import JsonFileOrganizer, is_ready
from ... import logger


class DocumentController(ABC):
    def __init__(self, **kwargs):
        self.suffix = kwargs.get("suffix", "")
        super().__init__()

    @abstractmethod
    def dump(self, doc: dict) -> dict:
        pass

    @abstractmethod
    def parse(self, item: dict):
        pass


class DocumentLoader(Process):
    def __init__(self, hosts: List, controller: DocumentController,
                task_queue: Queue, done_queue: Queue, **kwargs):
        super().__init__()
        self.hosts = hosts
        self.controller = controller
        self.task_queue = task_queue
        self.done_queue = done_queue
        self.err_dst = kwargs.get("err", os.path.devnull)

        # elasticsearch-py defaults
        self.chunk_size = kwargs.get("chunk_size", 500)
        self.max_bytes = kwargs.get("max_bytes", 100 * 1024 * 1024)
        self.threads = kwargs.get("threads", 4)
        self.timeout = kwargs.get("timeout", 10)

    def run(self):
        es = Elasticsearch(self.hosts, timeout=self.timeout)

        # Disable Elastic logger
        tracer = logging.getLogger("elasticsearch")
        tracer.setLevel(logging.CRITICAL + 1)

        with open(self.err_dst, "wt") as err:
            for filepath in iter(self.task_queue.get, None):
                with open(filepath, "rt") as fh:
                    documents = json.load(fh)

                bulk = helpers.parallel_bulk(
                    es, map(controller.dump, documents),
                    thread_count=self.threads,
                    queue_size=self.threads,
                    chunk_size=self.chunk_size,
                    max_chunk_bytes=self.max_bytes,
                    raise_on_exception=False,
                    raise_on_error=False
                )

                num_successful = 0
                failed = []
                for i, (status, item) in enumerate(bulk):
                    if status:
                        num_successful += 1
                    else:
                        controller.parse(item)
                        print(item)
                        print(i)
                        printe(documents[i])
                        raise RuntimeError()


                        try:
                            doc = items[_id]
                        except KeyError as exc:
                            print(item)
                            raise exc

                        failed.append(doc["_source"])

                        try:
                            _id = item["index"]["_id"]
                        except KeyError:
                            key = "update"
                            _id = item[key]["_id"]
                        else:
                            key = "index"
                        finally:
                            failed.append(items[_id]["_source"])

                        try:
                            """
                            Some items have a `data` property,
                            we do not need to store it in the err file
                            """
                            del item[key]["data"]
                        except KeyError:
                            pass
                        finally:
                            err.write("{}\n".format(item))

                self.done_queue.put((filepath, num_successful, failed))


def create_indices(hosts: List[str], indices: List[str], body_path: str,
                   doc_type: str, **kwargs):
    delete_all = kwargs.get("delete_all", False)
    num_shards = kwargs.get("num_shards", 5)
    shards_path = kwargs.get("shards_path")
    suffix = kwargs.get("suffix", "")

    # Load default settings and property mapping
    with open(body_path, "rt") as fh:
        body = json.load(fh)

    # Select desired mapping type
    body["mappings"] = {doc_type: body["mappings"][doc_type]}

    # Custom number of shards
    if shards_path:
        with open(shards_path, "rt") as fh:
            shards = json.load(fh)

        # Force keys to be in lower case
        shards = {k.lower(): v for k, v in shards.items()}
    else:
        shards = {}

    # Establish connection
    es = Elasticsearch(hosts, timeout=30)

    # Disable Elastic logger
    tracer = logging.getLogger("elasticsearch")
    tracer.setLevel(logging.CRITICAL + 1)

    if delete_all:
        # Delete all existing indices
        # to_delete = es.indices.get('*')
        # for index in to_delete:
        #     while True:
        #         try:
        #             res = es.indices.delete(index)
        #         except Exception as exc:
        #             logger.error("{}: {}".format(type(exc), exc))
        #             time.sleep(10)
        #         else:
        #             break
        es.indices.delete('*', allow_no_indices=True)

    # Create indices
    for index in indices:
        try:
            n_shards = shards[index]
        except KeyError:
            # No custom number of shards for this index
            n_shards = num_shards

        """
        Change settings for large bulk imports:

        https://www.elastic.co/guide/en/elasticsearch/reference/current/indices-update-settings.html
        https://www.elastic.co/guide/en/elasticsearch/guide/current/indexing-performance.html
        """
        body["settings"].update({
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
                raise exc
            except Exception as exc:
                logger.error("{}: {}".format(type(exc), exc))
                time.sleep(10)
            else:
                break


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


def organize_failed_docs(task_queue: Queue, done_queue: Queue,
                         dst: Optional[str]=None, write_back: bool=False):
    if dst:
        try:
            """
            Ensure the directory does not exist
            as we don't want files from a previous run to be considered
            """
            shutil.rmtree(dst)
        except FileNotFoundError:
            pass
        finally:
            os.makedirs(dst)
            organizer = JsonFileOrganizer(dst)
    else:
        organizer = None

    while True:
        total_successful = 0
        total_failed = 0
        num_files = 0
        for path, num_successful, failed in iter(task_queue.get, None):
            total_successful += num_successful
            total_failed += len(failed)

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

            num_files += 1
            if not num_files % 1000:
                logger.info("documents indexed: {:>15,} "
                            "({:,} failed)".format(total_successful,
                                                   total_failed))

        logger.info("documents indexed: {:>15,} "
                    "({:,} failed)".format(total_successful, total_failed))

        done_queue.put(total_failed)

        """
        Wait for instruction from parent:
            - True: we continue
            - False: max retries reached: quit
        """
        if not task_queue.get():
            break


def index_documents(hosts: List[str], controller: DocumentController,
                    src: str, **kwargs) -> bool:
    dst = kwargs.get("dst")
    max_retries = kwargs.get("max_retries", 0)
    processes = kwargs.get("processes", 1)
    raise_on_error = kwargs.get("raise_on_error", True)
    write_back = kwargs.get("write_back", False)

    file_queue = Queue()
    fail_queue = Queue(maxsize=processes)
    count_queue = Queue()
    organizer = Process(target=organize_failed_docs,
                        args=(fail_queue, count_queue, dst, write_back))
    organizer.start()

    processes = max(1, processes - 2)  # parent process + organizer
    num_retries = 0
    while True:
        logger.info("indexing documents (try #{})".format(num_retries + 1))
        workers = []

        for i in range(processes):
            w = DocumentLoader(hosts, controller, file_queue, fail_queue,
                               **kwargs)
            w.start()
            workers.append(w)

        for path in find_json_files(src):
            file_queue.put(path)

        for _ in workers:
            file_queue.put(None)

        # Wait for workers to complete
        for w in workers:
            w.join()

        # Inform organizer that we want the count of failed documents
        fail_queue.put(None)
        num_failed = count_queue.get()

        if num_failed and num_retries < max_retries:
            num_retries += 1
            fail_queue.put(True)   # continue
        else:
            fail_queue.put(False)  # stop
            organizer.join()
            break

    if num_failed:
        if raise_on_error:
            raise RuntimeError("{:,} documents not indexed".format(num_failed))
        else:
            logger.error("{:,} documents not indexed".format(num_failed))
            return False
    else:
        logger.info("complete")
        return True


def update_alias(hosts: List[str], indices: List[str], alias: str, **kwargs):
    delete_removed = kwargs.get("delete_removed", False)
    suffix = kwargs.get("suffix", "")

    new_indices = {index + suffix for index in indices}
    es = Elasticsearch(hosts, timeout=30)

    # Disable Elastic logger
    tracer = logging.getLogger("elasticsearch")
    tracer.setLevel(logging.CRITICAL + 1)

    exists = es.indices.exists_alias(name=alias)
    if exists:
        # Alias already exists: update it

        # Indices currently using the alias
        cur_indices = set(es.indices.get_alias(name=alias))

        actions = []
        for index in new_indices:
            try:
                # If passes: new index is already using the alias
                cur_indices.remove(index)
            except KeyError:
                # Otherwise, add the alias to the new index
                actions.append({
                    "add": {
                        "index": index,
                        "alias": alias
                    }
                })

        # Remove the alias from the current indices
        for index in cur_indices:
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

        if delete_removed:
            # Delete old indices that used the alias
            for index in cur_indices:
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
