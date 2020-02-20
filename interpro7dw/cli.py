#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import configparser
import os

from mundone import Task, Workflow

from interpro7dw import __version__
from interpro7dw.ebi.interpro import production as ippro
from interpro7dw.ebi.interpro import staging
from interpro7dw.ebi import pdbe, uniprot
from interpro7dw.utils import DataFiles


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

    ipr_pro_url = config["databases"]["production"]
    ipr_stg_url = config["databases"]["staging"]
    pfam_url = config["databases"]["pfam"]

    data_dir = config["stores"]["path"]

    lsf_queue = config["workflow"]["lsf_queue"]
    workflow_dir = config["workflow"]["path"]

    df = DataFiles(data_dir)
    tasks = [
        Task(fn=ippro.chunk_proteins,
            args=(ipr_pro_url, df.keys),
            name="init-export",
            scheduler=dict(mem=16000, queue=lsf_queue)
        ),

        # Export data from InterPro Oracle database
        Task(
            fn=ippro.export_features,
            args=(ipr_pro_url, df.keys, df.features),
            kwargs=dict(dir=data_dir, processes=8),
            name="export-features",
            requires=["init-export"],
            scheduler=dict(mem=4000, cpu=8, queue=lsf_queue)
        ),
        Task(
            fn=ippro.export_matches,
            args=(ipr_pro_url, df.keys, df.matches),
            kwargs=dict(dir=data_dir, processes=8),
            name="export-matches",
            requires=["init-export"],
            scheduler=dict(mem=8000, cpu=8, queue=lsf_queue)
        ),
        Task(
            fn=ippro.export_proteins,
            args=(ipr_pro_url, df.keys, df.proteins),
            kwargs=dict(dir=data_dir, processes=8),
            name="export-proteins",
            requires=["init-export"],
            scheduler=dict(mem=4000, cpu=8, queue=lsf_queue)
        ),
        Task(
            fn=ippro.export_residues,
            args=(ipr_pro_url, df.keys, df.residues),
            kwargs=dict(dir=data_dir, processes=8),
            name="export-residues",
            requires=["init-export"],
            scheduler=dict(mem=8000, cpu=8, queue=lsf_queue)
        ),
        Task(
            fn=ippro.export_sequences,
            args=(ipr_pro_url, df.keys, df.sequences),
            kwargs=dict(dir=data_dir, processes=8),
            name="export-sequences",
            requires=["init-export"],
            scheduler=dict(mem=4000, cpu=8, queue=lsf_queue)
        ),
        # Export data from UniProt Oracle database
        Task(
            fn=uniprot.export_comments,
            args=(ipr_pro_url, df.keys, df.comments),
            kwargs=dict(dir=data_dir, processes=8),
            name="export-comments",
            requires=["init-export"],
            scheduler=dict(mem=32000, cpu=8, queue=lsf_queue)
        ),
        Task(
            fn=uniprot.export_descriptions,
            args=(ipr_pro_url, df.keys, df.descriptions),
            kwargs=dict(dir=data_dir, processes=8),
            name="export-descriptions",
            requires=["init-export"],
            scheduler=dict(mem=32000, cpu=8, queue=lsf_queue)
        ),
        Task(
            fn=uniprot.export_evidences,
            args=(ipr_pro_url, df.keys, df.evidences),
            kwargs=dict(dir=data_dir, processes=8),
            name="export-evidences",
            requires=["init-export"],
            scheduler=dict(mem=32000, cpu=8, queue=lsf_queue)
        ),
        Task(
            fn=uniprot.export_proteomes,
            args=(ipr_pro_url, df.keys, df.proteomes),
            kwargs=dict(dir=data_dir, processes=8),
            name="export-proteomes",
            requires=["init-export"],
            scheduler=dict(mem=32000, cpu=8, queue=lsf_queue)
        ),


        Task(
            fn=ippro.export_taxonomy,
            args=(ipr_pro_url, df.taxonomy),
            name="export-taxonomy",
            scheduler=dict(mem=4000, queue=lsf_queue)
        ),
        Task(
            fn=uniprot.export_ref_proteomes,
            args=(ipr_pro_url, df.ref_proteomes),
            name="export-ref-proteomes",
            scheduler=dict(mem=100, queue=lsf_queue)
        ),
        Task(
            fn=staging.misc.insert_databases,
            args=(ipr_pro_url, ipr_stg_url),
            name="insert-databases",
            scheduler=dict(mem=100, queue=lsf_queue)
        ),
        Task(
            fn=ippro.export_entries,
            args=(ipr_pro_url, df.entries),
            name="export-entries",
            scheduler=dict(mem=1000, queue=lsf_queue)
        ),
        Task(
            fn=staging.entry.insert_annotations,
            args=(pfam_url, ipr_stg_url),
            name="insert-annotations",
            scheduler=dict(mem=4000, queue=lsf_queue)
        ),
        Task(
            fn=staging.entry.init_sets,
            args=(ipr_pro_url, ipr_stg_url, df.sets),
            name="init-sets",
            scheduler=dict(mem=16000, queue=lsf_queue)
        ),
        Task(
            fn=pdbe.export_structures,
            args=(ipr_pro_url, df.structures),
            name="export-structures",
            scheduler=dict(mem=8000, queue=lsf_queue)
        ),
        Task(
            fn=staging.protein.insert_isoforms,
            args=(df.entries, ipr_pro_url, ipr_stg_url),
            name="insert-isoforms",
            scheduler=dict(mem=4000, queue=lsf_queue),
            requires=["export-entries"]
        ),

        Task(
            fn=staging.protein.export_ida,
            args=(df.entries, df.matches, df.ida),
            kwargs=dict(dir=data_dir, processes=8),
            name="export-ida",
            scheduler=dict(mem=16000, queue=lsf_queue),
            requires=["export-entries", "export-matches"]
        ),

        # Task(
        #     fn=staging.proteome.init,
        #     args=(ipr_pro_url, ipr_stg_url),
        #     name="init-proteomes",
        #     scheduler=dict(mem=4000)
        # ),
        # Task(
        #     fn=staging.taxonomy.init,
        #     args=(ipr_pro_url, ipr_stg_url),
        #     name="init-taxonomy",
        #     scheduler=dict(mem=4000)
        # ),
    ]

    database = os.path.join(workflow_dir, f"{version}.sqlite")
    with Workflow(tasks, dir=workflow_dir, database=database) as workflow:
        workflow.run(args.tasks, dry_run=args.dry_run, monitor=not args.detach)
