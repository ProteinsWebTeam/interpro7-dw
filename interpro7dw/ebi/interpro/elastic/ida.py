# -*- coding: utf-8 -*-

import os
import shutil
from typing import Optional, Sequence

from interpro7dw import logger
from interpro7dw.utils import DirectoryTree, Store, datadump
from . import utils


BODY = {
    "mappings": {
        "properties": {
            "ida_id": {"type": "keyword"},
            "ida": {"type": "keyword"},
            "counts": {"type": "integer"}
        }
    }
}
INDEX = "ida"

# Aliases
STAGING = "ida_staging"
LIVE = "ida_current"
PREVIOUS = "ida_previous"


def dump_documents(src_uniprot2ida: str, outdir: str, cache_size: int=1000000):
    logger.info("preparing data")
    try:
        shutil.rmtree(outdir)
    except FileNotFoundError:
        pass
    finally:
        os.makedirs(outdir)
        organizer = DirectoryTree(outdir)
        open(os.path.join(outdir, utils.LOADING), "w").close()

    uniprot2ida = Store(src_uniprot2ida)

    ida_count = {}
    for dom_arch, dom_arch_id in uniprot2ida.values():
        try:
            ida_count[dom_arch_id] += 1
        except KeyError:
            ida_count[dom_arch_id] = 1

    logger.info("starting")
    i = 0
    num_documents = 0
    cached_documents = []
    for dom_arch, dom_arch_id in uniprot2ida.values():
        cached_documents.append({
            "ida_id": dom_arch_id,
            "ida": dom_arch,
            "counts": ida_count[dom_arch_id]
        })

        if len(cached_documents) == cache_size:
            filepath = organizer.mktemp()
            datadump(filepath, cached_documents[:cache_size])
            os.rename(filepath, f"{filepath}{utils.EXTENSION}")
            cached_documents = []
            num_documents += cache_size

        i += 1
        if not i % 10000000:
            logger.info(f"{i:>12,}")

    logger.info(f"{i:>12,}")

    num_documents += len(cached_documents)
    if cached_documents:
        filepath = organizer.mktemp()
        datadump(filepath, cached_documents[:cache_size])
        os.rename(filepath, f"{filepath}{utils.EXTENSION}")

    uniprot2ida.close()

    # Delete flag file to notify loaders that all files are ready
    os.remove(os.path.join(outdir, utils.LOADING))

    logger.info(f"complete ({num_documents:,} documents)")


def index_documents(hosts: Sequence[str], indir: str, version: str,
                    outdir: Optional[str]=None, writeback: bool=False):
    index = f"{INDEX}{version}"

    def wrap(doc: dict) -> dict:
        return {
            "_op_type": "index",
            "_index": index,
            "_id": doc["ida_id"],
            "_source": doc
        }

    es = utils.connect(hosts, verbose=False)
    if outdir:
        logger.info("creating indices")
        body = BODY.copy()
        body["settings"] = {
            "index": {
                # Static settings
                "number_of_shards": utils.DEFAULT_SHARDS,

                # Dynamic settings
                "number_of_replicas": 0,  # defaults to 1
                "refresh_interval": -1  # defaults to 1s
            }
        }

        utils.create_index(es, index, body)


    utils.index_documents(es, indir, callback=wrap, outdir=outdir,
                          writeback=writeback)

    utils.add_alias(es, [index], STAGING, delete_indices=False)


def publish(hosts: Sequence[str]):
    es = utils.connect(hosts, verbose=False)

    # Make LIVE indices pointed by PREVIOUS
    live = es.indices.get_alias(name=LIVE)
    utils.add_alias(es, live, PREVIOUS, delete_indices=True)

    # Make STAGING indices pointed by LIVE
    staging = es.indices.get_alias(name=STAGING)
    utils.add_alias(es, staging, LIVE, delete_indices=False)
