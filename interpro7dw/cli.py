#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import configparser
import os

from mundone import Task, Workflow

from interpro7dw import __version__
from interpro7dw.ebi.interpro import production as ippro
from interpro7dw.ebi import uniprot


def build():
    parser = argparse.ArgumentParser(
        description="Build InterPro7 data warehouse"
    )
    parser.add_argument("config",
                        metavar="config.ini",
                        help="configuration file")
    parser.add_argument("-t", "--tasks",
                        nargs="*",
                        default=None,
                        metavar="TASK",
                        help="tasks to run")
    parser.add_argument("--dry-run",
                        action="store_true",
                        default=False,
                        help="list tasks to run and exit")
    parser.add_argument("--resume",
                        action="store_true",
                        help="skip completed tasks")
    parser.add_argument("--detach",
                        action="store_true",
                        help="enqueue tasks to run and exit")
    parser.add_argument("-v", "--version", action="version",
                        version=f"%(prog)s {__version__}",
                        help="show the version and exit")
    args = parser.parse_args()

    if not os.path.isfile(args.config):
        parser.error(f"cannot open '{args.config}': no such file or directory")

    config = configparser.ConfigParser()
    config.read(args.config)

    version = config["release"]["version"]

    ipr_pro = config["databases"]["production"]

    stores_dir = config["stores"]["path"]

    lsf_queue = config["workflow"]["lsf_queue"]
    workflow_dir = config["workflow"]["path"]

    os.makedirs(stores_dir, exist_ok=True)
    tasks = [
        Task(fn=ippro.chunk_proteins,
            args=(ipr_pro, os.path.join(stores_dir, "keys")),
            name="init-export",
            scheduler=dict(mem=16000, queue=lsf_queue)
        ),

        # Export data from InterPro Oracle database
        Task(
            fn=ippro.export_features,
            args=(ipr_pro,
                  os.path.join(stores_dir, "keys"),
                  os.path.join(stores_dir, "features")),
            kwargs=dict(dir=stores_dir, processes=8),
            name="export-features",
            requires=["init-export"],
            scheduler=dict(mem=4000, cpu=8)
        ),
        Task(
            fn=ippro.export_matches,
            args=(ipr_pro,
                  os.path.join(stores_dir, "keys"),
                  os.path.join(stores_dir, "matches")),
            kwargs=dict(dir=stores_dir, processes=8),
            name="export-matches",
            requires=["init-export"],
            scheduler=dict(mem=8000, cpu=8)
        ),
        Task(
            fn=ippro.export_proteins,
            args=(ipr_pro,
                  os.path.join(stores_dir, "keys"),
                  os.path.join(stores_dir, "proteins")),
            kwargs=dict(dir=stores_dir, processes=8),
            name="export-proteins",
            requires=["init-export"],
            scheduler=dict(mem=4000, cpu=8)
        ),
        Task(
            fn=ippro.export_residues,
            args=(ipr_pro,
                  os.path.join(stores_dir, "keys"),
                  os.path.join(stores_dir, "residues")),
            kwargs=dict(dir=stores_dir, processes=8),
            name="export-residues",
            requires=["init-export"],
            scheduler=dict(mem=8000, cpu=8)
        ),
        Task(
            fn=ippro.export_sequences,
            args=(ipr_pro,
                  os.path.join(stores_dir, "keys"),
                  os.path.join(stores_dir, "sequences")),
            kwargs=dict(dir=stores_dir, processes=8),
            name="export-sequences",
            requires=["init-export"],
            scheduler=dict(mem=4000, cpu=8)
        ),
        # Export data from UniProt Oracle database
        Task(
            fn=uniprot.export_comments,
            args=(ipr_pro,
                  os.path.join(stores_dir, "keys"),
                  os.path.join(stores_dir, "comments")),
            kwargs=dict(dir=stores_dir, processes=8),
            name="export-comments",
            requires=["init-export"],
            scheduler=dict(mem=32000, cpu=8)
        ),
        Task(
            fn=uniprot.export_descriptions,
            args=(ipr_pro,
                  os.path.join(stores_dir, "keys"),
                  os.path.join(stores_dir, "descriptions")),
            kwargs=dict(dir=stores_dir, processes=8),
            name="export-descriptions",
            requires=["init-export"],
            scheduler=dict(mem=32000, cpu=8)
        ),
        Task(
            fn=uniprot.export_evidences,
            args=(ipr_pro,
                  os.path.join(stores_dir, "keys"),
                  os.path.join(stores_dir, "evidences")),
            kwargs=dict(dir=stores_dir, processes=8),
            name="export-evidences",
            requires=["init-export"],
            scheduler=dict(mem=32000, cpu=8)
        ),
        Task(
            fn=uniprot.export_proteomes,
            args=(ipr_pro,
                  os.path.join(stores_dir, "keys"),
                  os.path.join(stores_dir, "proteomes")),
            kwargs=dict(dir=stores_dir, processes=8),
            name="export-proteomes",
            requires=["init-export"],
            scheduler=dict(mem=32000, cpu=8)
        ),
    ]

    database = os.path.join(workflow_dir, f"{version}.sqlite")
    with Workflow(tasks, dir=workflow_dir, database=database) as workflow:
        workflow.run(args.tasks, dry_run=args.dry_run, monitor=not args.detach)
