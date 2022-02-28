import glob
import logging
import os
import pickle
import time

from elasticsearch import Elasticsearch, exceptions
from elasticsearch.helpers import parallel_bulk as pbulk

from interpro7dw.utils import logger
from . import config


def connect(hosts: list[str], timeout: int = 10,
            verbose: bool = True) -> Elasticsearch:
    es = Elasticsearch(hosts=hosts, timeout=timeout)

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


def add_alias(es: Elasticsearch, indices: list[str], alias: str):
    if es.indices.exists_alias(name=alias):
        # Alias already exists: update it

        # Indices currently pointed by the alias
        old_indices = set(es.indices.get_alias(name=alias))

        actions = []
        for index in set(indices):
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


def create_indices(databases_file: str, hosts: list[str], version: str):
    es = connect(hosts, verbose=False)

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

        num_shards = config.REL_SHARDS.get(index, config.DEFAULT_SHARDS)
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
            alias2indices[alias].append(index)
        except KeyError:
            alias2indices[alias] = [index]

    # Add an 'staging' alias to all newly created indices
    for alias, new_indices in alias2indices.items():
        add_alias(es, new_indices, alias + config.STAGING_ALIAS_SUFFIX)


def iter_files(root: str, version: str):
    load_sentinel = os.path.join(root, f"{version}{config.LOAD_SUFFIX}")
    done_sentinel = os.path.join(root, f"{version}{config.DONE_SUFFIX}")

    if not os.path.isfile(done_sentinel):
        while not os.path.isfile(load_sentinel):
            # Wait until files start being generated
            time.sleep(60)

    pathname = os.path.join(root, "**", f"*{config.EXTENSION}")
    files = set()
    active = True
    while True:
        for path in glob.iglob(pathname, recursive=True):
            if path in files:
                continue

            files.add(path)
            yield path

        if not active:
            break
        elif os.path.isfile(done_sentinel):
            # All files ready: they will all be found in the next iteration
            active = False
        else:
            # Files are still being written
            time.sleep(60)


def index_documents(hosts: list[str], indir: str, version: str,
                    threads: int = 4):
    logger.info("starting")
    kwargs = {
        "thread_count": threads,
        "queue_size": threads,
        "raise_on_exception": False,
        "raise_on_error": False
    }

    es = connect(hosts, timeout=30, verbose=False)
    num_documents = 0
    num_indexed = 0
    first_pass = True
    while True:
        for filepath in iter_files(indir, version):
            with open(filepath, "rb") as fh:
                documents = pickle.load(fh)

            if first_pass:
                # Count only once the number of documents to index
                num_documents += len(documents)

            actions = []
            for idx, doc_id, doc in documents:
                actions.append({
                    "_op_type": "index",
                    "_index": idx,
                    "_id": doc_id,
                    "_source": doc
                })

            failed = []
            for i, (ok, info) in enumerate(pbulk(es, actions, **kwargs)):
                if ok:
                    num_indexed += 1
                    if not num_indexed % 100e6:
                        logger.info(f"{num_indexed:>15,}")
                else:
                    failed.append(documents[i])

                    # try:
                    #     is_429 = info["index"]["status"] == 429
                    # except (KeyError, IndexError):
                    #     is_429 = False
                    #
                    # try:
                    #     exc = info["index"]["exception"]
                    # except (KeyError, TypeError):
                    #     exc = None
                    #
                    # if is_429 or isinstance(exc, exceptions.ConnectionTimeout):
                    #     pause = True
                    # else:
                    #     logger.debug(info)

            if failed:
                # Overwrite file with failed documents
                with open(filepath, "wb") as fh:
                    pickle.dump(failed, fh)
            else:
                # Remove file as all documents have been successfully indexed
                os.unlink(filepath)

        logger.info(f"{num_indexed:>15,}")
        first_pass = False

        if num_indexed == num_documents:
            break

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

    logger.info("done")


def publish(hosts: list[str]):
    es = connect(hosts, verbose=False)

    for alias in (config.IDA_ALIAS, config.REL_ALIAS):
        live_alias = alias + config.LIVE_ALIAS_SUFFIX

        # Add the 'previous' alias to current 'live' indices
        if es.indices.exists_alias(name=live_alias):
            indices = es.indices.get_alias(name=live_alias)

            prev_alias = alias + config.PREVIOUS_ALIAS_SUFFIX
            add_alias(es, indices, prev_alias)

        # Add the 'live' alias to current 'staging' indices
        staging_alias = alias + config.STAGING_ALIAS_SUFFIX
        indices = es.indices.get_alias(name=staging_alias)
        add_alias(es, indices, live_alias)
