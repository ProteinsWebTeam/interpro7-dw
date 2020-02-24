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


class DataFiles(object):
    def __init__(self, path: str):
        os.makedirs(path, exist_ok=True)

        self.clans = os.path.join(path, "sets")
        self.comments = os.path.join(path, "comments")
        self.descriptions = os.path.join(path, "descriptions")
        self.entries = os.path.join(path, "entries")
        self.evidences = os.path.join(path, "evidences")
        self.features = os.path.join(path, "features")
        self.keys = os.path.join(path, "keys")
        self.ida = os.path.join(path, "ida")
        self.matches = os.path.join(path, "matches")
        self.overlapping = os.path.join(path, "overlapping")
        self.proteins = os.path.join(path, "proteins")
        self.proteomes = os.path.join(path, "proteomes")
        self.ref_proteomes = os.path.join(path, "ref-proteomes")
        self.residues = os.path.join(path, "residues")
        self.sequences = os.path.join(path, "sequences")
        self.structures = os.path.join(path, "structures")
        self.taxonomy = os.path.join(path, "taxonomy")


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
            fn=staging.insert_databases,
            args=(ipr_pro_url, ipr_stg_url),
            name="insert-databases",
            scheduler=dict(mem=100, queue=lsf_queue)
        ),
        Task(
            fn=staging.insert_annotations,
            args=(pfam_url, ipr_stg_url),
            name="insert-annotations",
            scheduler=dict(mem=4000, queue=lsf_queue)
        ),
        Task(
            fn=staging.init_clans,
            args=(ipr_pro_url, ipr_stg_url, df.clans),
            name="init-sets",
            scheduler=dict(mem=1000, queue=lsf_queue)
        ),
        Task(
            fn=ippro.export_entries,
            args=(ipr_pro_url, df.clans, df.entries),
            name="export-entries",
            scheduler=dict(mem=1000, queue=lsf_queue),
            requires=["init-sets"]
        ),
        Task(
            fn=pdbe.export_structures,
            args=(ipr_pro_url, df.structures),
            name="export-structures",
            scheduler=dict(mem=8000, queue=lsf_queue)
        ),
        Task(
            fn=staging.insert_isoforms,
            args=(df.entries, ipr_pro_url, ipr_stg_url),
            name="insert-isoforms",
            scheduler=dict(mem=4000, queue=lsf_queue),
            requires=["export-entries"]
        ),

        Task(
            fn=staging.export_ida,
            args=(df.entries, df.matches, df.ida),
            kwargs=dict(dir=data_dir, processes=8),
            name="export-ida",
            scheduler=dict(mem=16000, cpu=8, queue=lsf_queue),
            requires=["export-entries", "export-matches"]
        ),

        Task(
            fn=staging.export_overlapping_entries,
            args=(df.entries, df.matches, df.overlapping),
            # kwargs=dict(url=ipr_pro_url),
            name="overlapping-entries",
            scheduler=dict(mem=4000, queue=lsf_queue),
            requires=["export-entries", "export-matches"]
        ),

        Task(
            fn=staging.insert_proteins,
            args=(df.proteins, df.comments, df.descriptions, df.entries,
                  df.evidences, df.features, df.ida, df.matches, df.proteomes,
                  df.residues, df.sequences, df.structures, df.taxonomy,
                  ipr_pro_url, ipr_stg_url),
            name="insert-proteins",
            scheduler=dict(mem=16000, queue=lsf_queue),
            requires=["export-proteins", "export-comments",
                      "export-descriptions", "export-entries",
                      "export-evidences", "export-features", "export-ida",
                      "export-matches", "export-proteomes", "export-residues",
                      "export-sequences", "export-structures",
                      "export-taxonomy"]
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
