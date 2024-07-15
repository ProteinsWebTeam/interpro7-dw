import glob
import logging
import multiprocessing as mp
import os
import pickle
import time

from elasticsearch import Elasticsearch, exceptions
from elasticsearch.helpers import streaming_bulk as sbulk

from interpro7dw.utils import logger
from . import config


def connect(hosts: list[str], user: str, password: str, fingerprint: str,
            timeout: int = 10, verbose: bool = True) -> Elasticsearch:
    _hosts = []
    for host in hosts:
        if host.startswith("https"):
            _hosts.append(host)
        else:
            _hosts.append(f"https://{host}")

    es = Elasticsearch(hosts=_hosts,
                       ssl_assert_fingerprint=fingerprint,
                       basic_auth=(user, password),
                       request_timeout=timeout)

    if not verbose:
        # Disable Elastic logger
        tracer = logging.getLogger("elasticsearch")
        tracer.setLevel(logging.CRITICAL + 1)

    return es


def delete_index(es: Elasticsearch, name: str):
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


def add_alias(es: Elasticsearch, indices: set[str], alias: str):
    if es.indices.exists_alias(name=alias):
        # Alias already exists: update it

        # Indices currently pointed by the alias
        old_indices = set(es.indices.get_alias(name=alias))

        actions = []
        for index in indices:
            if index in old_indices:
                # Index is already pointed by alias
                old_indices.remove(index)
            else:
                # Add alias to index
                actions.append({"add": {"index": index, "alias": alias}})

        # Remove alias from old indices
        for index in old_indices:
            actions.append({"remove": {"index": index, "alias": alias}})

        if actions:
            """
            Atomic operation:
            Alias removed from the old indices
            at the same time it's added to the new ones
            """
            es.indices.update_aliases(body={"actions": actions})
    else:
        # Creat new alias, then point to indices
        es.indices.put_alias(index=','.join(indices), name=alias)


def create_indices(databases_file: str, hosts: list[str], user: str,
                   password: str, fingerprint: str, version: str,
                   suffix: str = ""):
    es = connect(hosts, user, password, fingerprint, verbose=False)

    """
    Assuming we are creating indices for version 100.0:
      * version 99.0 should have the 'live' alias
      * version 98.0 should have the 'previous' alias: delete the indices
    """
    for alias in (config.IDA_ALIAS, config.REL_ALIAS):
        alias += config.PREVIOUS_ALIAS_SUFFIX
        if es.indices.exists_alias(name=alias):
            for idx in es.indices.get_alias(name=alias):
                delete_index(es, idx)

    # Create a list of all new indices to create
    indices = [
        (config.IDA_INDEX, config.IDA_BODY, config.IDA_ALIAS),
        (config.REL_DEFAULT_INDEX, config.REL_BODY, config.REL_ALIAS)
    ]
    with open(databases_file, "rb") as fh:
        for key, info in pickle.load(fh).items():
            if info["type"] == "entry":
                indices.append((key.lower(),
                                config.REL_BODY,
                                config.REL_ALIAS))

    # Create new indices
    alias2indices = {}
    for index, body, alias in indices:
        try:
            settings = body["settings"]
        except KeyError:
            settings = body["settings"] = {}

        num_shards = config.NUM_SHARDS.get(index, config.DEFAULT_SHARDS)
        settings.update({
            "index": {
                # Static settings
                "number_of_shards": num_shards,
                "codec": "best_compression",

                # Dynamic settings
                "number_of_replicas": 0,  # defaults to 1
                "refresh_interval": -1  # defaults to 1s
            }
        })

        index += version  # Use InterPro version as suffix
        index += suffix   # Additional debugging suffix

        delete_index(es, index)
        while True:
            try:
                es.indices.create(index=index, body=body)
            except exceptions.ConnectionTimeout:
                time.sleep(30)
            except exceptions.RequestError as exc:
                raise exc
            except Exception as exc:
                logger.warning(f"{index}: {exc}")
                time.sleep(30)
            else:
                break

        try:
            alias2indices[alias].add(index)
        except KeyError:
            alias2indices[alias] = {index}

    # Add an 'staging' alias to all newly created indices
    for alias, new_indices in alias2indices.items():
        add_alias(es, new_indices, alias + config.STAGING_ALIAS_SUFFIX)


def update_indices(hosts: list[str], user: str, password: str,
                   fingerprint: str):
    es = connect(hosts, user, password, fingerprint, verbose=False)

    # Update index settings
    for alias in (config.IDA_ALIAS, config.REL_ALIAS):
        alias += config.STAGING_ALIAS_SUFFIX

        # This assumes there are indices with the 'staging' alias
        for index in es.indices.get_alias(name=alias):
            es.indices.put_settings(
                body={
                    "number_of_replicas": 1,
                    "refresh_interval": None  # default (1s)
                },
                index=index
            )


def iter_files(root: str, version: str):
    load_sentinel = os.path.join(root, f"{version}{config.LOAD_SUFFIX}")
    done_sentinel = os.path.join(root, f"{version}{config.DONE_SUFFIX}")

    if not os.path.isfile(done_sentinel):
        logger.info("pending")
        while not os.path.isfile(load_sentinel):
            # Wait until files start being generated
            time.sleep(15)

    pathname = os.path.join(root, "**", f"*{config.EXTENSION}")
    files = set()
    done = os.path.isfile(done_sentinel)
    while True:
        for path in glob.iglob(pathname, recursive=True):
            if path in files:
                continue

            files.add(path)
            yield path

        if done:
            break
        elif os.path.isfile(done_sentinel):
            # All files ready: last iteration
            done = True
        else:
            # Files are still being written
            time.sleep(60)


def index_documents(hosts: list[str], user: str, password: str,
                    fingerprint: str, indir: str, version: str,
                    processes: int = 8, suffix: str = ""):
    logger.info("starting")
    progress = 0
    milestone = step = 1e8
    while True:
        inqueue = mp.Queue()
        outqueue = mp.Queue()
        indexers = []
        for _ in range(max(1, processes - 2)):
            p = mp.Process(
                target=run_consumer,
                args=(hosts, user, password, fingerprint, inqueue, outqueue,
                      suffix),
            )
            p.start()
            indexers.append(p)

        producer = mp.Process(target=run_producer,
                              args=(indir, version, inqueue, len(indexers)))
        producer.start()

        running = len(indexers)
        files_processed = 0
        while running:
            ongoing, count = outqueue.get()
            if ongoing:
                progress += count
                if progress >= milestone:
                    logger.info(f"{progress:>15,}")
                    milestone += step
            else:
                files_processed += count
                running -= 1

        if not files_processed:
            break

    logger.info(f"{progress:>15,}")
    update_indices(hosts, user, password, fingerprint)
    logger.info("done")


def run_producer(indir: str, version: str, queue: mp.Queue, num_workers: int):
    for filepath in iter_files(indir, version):
        queue.put(filepath)

    for _ in range(num_workers):
        queue.put(None)


def run_consumer(hosts: list[str], user: str, password: str, fingerprint: str,
                 inqueue: mp.Queue, outqueue: mp.Queue, suffix: str = ""):
    es = connect(hosts, user, password, fingerprint,
                 timeout=120, verbose=False)
    files = 0
    for filepath in iter(inqueue.get, None):
        files += 1

        # Load documents from file
        with open(filepath, "rb") as fh:
            documents = pickle.load(fh)

        # Prepare body for bulk request
        actions = []
        for idx, doc_id, doc in documents:
            actions.append({
                "_op_type": "index",
                "_index": idx + suffix,
                "_id": doc_id,
                "_source": doc
            })

        # Index documents
        failed = []
        n_docs = n_indexed = 0
        results = sbulk(es, actions, max_retries=5,
                        raise_on_error=False, raise_on_exception=False)
        for i, (ok, info) in enumerate(results):
            n_docs += 1
            if ok:
                n_indexed += 1
            else:
                failed.append(documents[i])

        if n_indexed == n_docs:
            # Remove file as all documents have been successfully indexed
            assert len(failed) == 0
            os.unlink(filepath)
        elif failed:
            # Overwrite file with failed documents
            with open(filepath, "wb") as fh:
                pickle.dump(failed, fh)
        else:
            raise RuntimeError(f"{filepath}: {n_docs} documents, "
                               f"{n_indexed} indexed, {len(failed)} failed")

        outqueue.put((True, n_indexed))

    outqueue.put((False, files))


def publish(hosts: list[str], user: str, password: str, fingerprint: str):
    es = connect(hosts, user, password, fingerprint, timeout=60, verbose=False)

    for alias in (config.IDA_ALIAS, config.REL_ALIAS):
        staging_alias = alias + config.STAGING_ALIAS_SUFFIX
        live_alias = alias + config.LIVE_ALIAS_SUFFIX

        # ObjectApiResponse, which is like a dict (index as str -> dict)
        response = es.indices.get_alias(name=staging_alias)
        staging_indices = set(response)

        add_previous = set()
        if es.indices.exists_alias(name=live_alias):
            response = es.indices.get_alias(name=live_alias)
            live_indices = set(response)

            # If an index is in staging and live, do not add it to previous
            add_previous |= live_indices - staging_indices

        if add_previous:
            prev_alias = alias + config.PREVIOUS_ALIAS_SUFFIX
            add_alias(es, add_previous, prev_alias)

        add_alias(es, staging_indices, live_alias)
