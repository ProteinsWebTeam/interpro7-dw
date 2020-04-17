# -*- coding: utf-8 -*-

import os
import shutil
from typing import Optional, Sequence

from interpro7dw import logger
from interpro7dw.utils import DirectoryTree, Store, dumpobj
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


def dump_documents(src_uniprot2ida: str, outdir: str, version: str,
                   cache_size: int=100000):
    logger.info("preparing data")
    try:
        shutil.rmtree(outdir)
    except FileNotFoundError:
        pass
    finally:
        os.makedirs(outdir)
        organizer = DirectoryTree(outdir)
        open(os.path.join(outdir, f"{version}{utils.LOAD_SUFFIX}"), "w").close()

    uniprot2ida = Store(src_uniprot2ida)

    logger.info("loading domain architectures")
    domains = {}
    for dom_members, dom_arch, dom_arch_id in uniprot2ida.values():
        try:
            dom = domains[dom_arch_id]
        except KeyError:
            domains[dom_arch_id] = {
                "ida_id": dom_arch_id,
                "ida": dom_arch,
                "counts": 1
            }
        else:
            dom["counts"] += 1

    logger.info("writing documents")
    domains = list(domains.values())
    for i in range(0, len(domains), cache_size):
        filepath = organizer.mktemp()
        dumpobj(filepath, domains[i:i+cache_size])
        os.rename(filepath, f"{filepath}{utils.EXTENSION}")

    uniprot2ida.close()

    open(os.path.join(outdir, f"{version}{utils.DONE_SUFFIX}"), "w").close()
    os.remove(os.path.join(outdir, f"{version}{utils.LOAD_SUFFIX}"))

    logger.info(f"complete ({len(domains):,} documents)")


def index_documents(hosts: Sequence[str], indir: str, version: str,
                    outdir: Optional[str]=None, writeback: bool=False,
                    create_indices: bool=True):
    index = f"{INDEX}{version}"

    def wrap(doc: dict) -> dict:
        return {
            "_op_type": "index",
            "_index": index,
            "_id": doc["ida_id"],
            "_source": doc
        }

    es = utils.connect(hosts, verbose=False)
    if create_indices:
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

    if es.indices.exists_alias(name=PREVIOUS):
        for prev_index in es.indices.get_alias(name=PREVIOUS):
            utils.delete_index(es, prev_index)

    utils.index_documents(es, indir, version, callback=wrap, outdir=outdir,
                          writeback=writeback)

    utils.add_alias(es, [index], STAGING)


def publish(hosts: Sequence[str]):
    utils.publish(hosts, STAGING, LIVE, PREVIOUS)
