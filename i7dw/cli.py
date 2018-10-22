#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import configparser
import os

from mundone import Task, Workflow

from i7dw import __version__, elastic, mysql
from i7dw.ebi import ebisearch, goa, interpro, uniprot


def cli():
    parser = argparse.ArgumentParser(
        description="Build InterPro7 data warehouse"
    )
    parser.add_argument("config",
                        metavar="config.ini",
                        help="configuration file")
    parser.add_argument("-t", "--tasks",
                        nargs="*",
                        default=None,
                        help="tasks to run")
    parser.add_argument("--dry-run",
                        action="store_true",
                        default=False,
                        help="list tasks to run and quit")
    parser.add_argument("--resume",
                        action="store_true",
                        default=False,
                        help="skip completed tasks")
    parser.add_argument("--detach",
                        action="store_true",
                        default=False,
                        help="enqueue tasks to run and quit")
    parser.add_argument("--retry",
                        action="store_true",
                        default=False,
                        help="rerun failed tasks (once)")
    parser.add_argument("--daemon",
                        action="store_true",
                        default=False,
                        help="do not kill running tasks "
                             "when exiting the program")
    parser.add_argument("-v", "--version", action="version",
                        version="%(prog)s {}".format(__version__),
                        help="show the version and quit")
    args = parser.parse_args()

    if not os.path.isfile(args.config):
        raise RuntimeError(
            "cannot open '{}': "
            "no such file or directory".format(args.config)
        )

    config = configparser.ConfigParser()
    config.read(args.config)

    export_dir = config["export"]["dir"]
    if not os.path.isdir(export_dir):
        os.makedirs(export_dir)

    ora_ipro = config["databases"]["interpro_oracle"]
    my_ipro_stg = config["databases"]["interpro_mysql_stg"]
    my_ipro_rel = config["databases"]["interpro_mysql_rel"]
    my_pfam = config["databases"]["pfam_mysql"]
    queue = config["workflow"]["queue"]

    elastic_hosts = config["elastic"]["hosts"].split(',')
    elastic_dir = config["elastic"]["dir"]

    threshold = config.getfloat("jaccard", "threshold")

    tasks = [
        # Export data to stores
        Task(
            name="export-comments",
            fn=uniprot.export_protein_comments,
            args=(ora_ipro, os.path.join(export_dir, "comments.bs")),
            kwargs=dict(chunk_size=100000),
            scheduler=dict(queue=queue, mem=1000)
        ),
        Task(
            name="export-descriptions",
            fn=uniprot.export_protein_descriptions,
            args=(ora_ipro, os.path.join(export_dir, "descriptions.bs")),
            kwargs=dict(chunk_size=100000),
            scheduler=dict(queue=queue, mem=1000)
        ),
        Task(
            name="export-evidences",
            fn=uniprot.export_protein_evidence,
            args=(ora_ipro, os.path.join(export_dir, "evidences.bs")),
            kwargs=dict(chunk_size=100000),
            scheduler=dict(queue=queue, mem=500)
        ),
        Task(
            name="export-genes",
            fn=uniprot.export_protein_gene,
            args=(ora_ipro, os.path.join(export_dir, "genes.bs")),
            kwargs=dict(chunk_size=100000),
            scheduler=dict(queue=queue, mem=500)
        ),
        Task(
            name="export-proteomes",
            fn=uniprot.export_protein_proteomes,
            args=(ora_ipro, os.path.join(export_dir, "proteomes.bs")),
            kwargs=dict(chunk_size=100000),
            scheduler=dict(queue=queue, mem=500)
        ),
        Task(
            name="export-annotations",
            fn=goa.export_annotations,
            args=(ora_ipro, os.path.join(export_dir, "annotations.bs")),
            kwargs=dict(chunk_size=100000),
            scheduler=dict(queue=queue, mem=2000)
        ),
        Task(
            name="export-structures",
            fn=interpro.export_struct_matches,
            args=(ora_ipro, os.path.join(export_dir, "struct_matches.bs")),
            kwargs=dict(chunk_size=100000),
            scheduler=dict(queue=queue, mem=4000)
        ),
        Task(
            name="export-matches",
            fn=interpro.export_prot_matches,
            args=(ora_ipro, os.path.join(export_dir, "prot_matches.bs")),
            kwargs=dict(chunk_size=100000),
            scheduler=dict(queue=queue, mem=3000)
        ),
        Task(
            name="export-features",
            fn=interpro.export_prot_matches_extra,
            args=(ora_ipro, os.path.join(export_dir, "prot_matches_extra.bs")),
            kwargs=dict(chunk_size=100000),
            scheduler=dict(queue=queue, mem=2000)
        ),
        Task(
            name="export-residues",
            fn=interpro.export_residues,
            args=(ora_ipro, os.path.join(export_dir, "residues.bs")),
            kwargs=dict(chunk_size=100000),
            scheduler=dict(queue=queue, mem=3000)
        ),
        Task(
            name="export-proteins",
            fn=interpro.export_proteins,
            args=(
                ora_ipro, 
                os.path.join(export_dir, "proteins.bs"), 
                os.path.join(export_dir, "sequences.bs")
            ),
            kwargs=dict(chunk_size=100000),
            scheduler=dict(queue=queue, mem=2000)
        ),

        # Load data to MySQL
        Task(
            name="init-tables",
            fn=mysql.init_tables,
            args=(my_ipro_stg,),
            scheduler=dict(queue=queue)
        ),
        Task(
            name="insert-taxa",
            fn=mysql.insert_taxa,
            args=(ora_ipro, my_ipro_stg),
            scheduler=dict(queue=queue, mem=4000),
            requires=["init-tables"]
        ),
        Task(
            name="insert-proteomes",
            fn=mysql.insert_proteomes,
            args=(ora_ipro, my_ipro_stg),
            scheduler=dict(queue=queue, mem=1000),
            requires=["insert-taxa"]
        ),
        Task(
            name="insert-databases",
            fn=mysql.insert_databases,
            args=(ora_ipro, my_ipro_stg),
            scheduler=dict(queue=queue),
            requires=["init-tables"]
        ),
        Task(
            name="insert-entries",
            fn=mysql.insert_entries,
            args=(ora_ipro, my_pfam, my_ipro_stg),
            scheduler=dict(queue=queue, mem=2000),
            requires=["insert-databases"]
        ),
        Task(
            name="insert-annotations",
            fn=mysql.insert_annotations,
            args=(my_pfam, my_ipro_stg),
            scheduler=dict(queue=queue, mem=4000),
            requires=["insert-entries"]
        ),
        Task(
            name="insert-structures",
            fn=mysql.insert_structures,
            args=(ora_ipro, my_ipro_stg),
            scheduler=dict(queue=queue, mem=1000),
            requires=["insert-databases"]
        ),
        Task(
            name="insert-sets",
            fn=mysql.insert_sets,
            args=(my_pfam, my_ipro_stg),
            scheduler=dict(queue=queue),
            requires=["insert-databases"]
        ),
        Task(
            name="insert-proteins",
            fn=mysql.insert_proteins,
            args=(
                my_ipro_stg,
                os.path.join(export_dir, "proteins.bs"),
                os.path.join(export_dir, "sequences.bs"),
                os.path.join(export_dir, "evidences.bs"),
                os.path.join(export_dir, "descriptions.bs"),
                os.path.join(export_dir, "comments.bs"),
                os.path.join(export_dir, "proteomes.bs"),
                os.path.join(export_dir, "genes.bs"),
                os.path.join(export_dir, "annotations.bs"),
                os.path.join(export_dir, "residues.bs"),
                os.path.join(export_dir, "struct_matches.bs"),
                os.path.join(export_dir, "prot_matches_extra.bs")
            ),
            scheduler=dict(queue=queue, mem=32000),
            requires=[
                "insert-taxa", "export-proteins", "export-evidences",
                "export-descriptions", "export-comments",
                "export-proteomes", "export-genes",
                "export-annotations", "export-residues",
                "export-structures", "export-features"
            ]
        ),

        # Release notes
        Task(
            name="release-notes",
            fn=mysql.make_release_notes,
            args=(
                my_ipro_stg,
                my_ipro_rel,
                os.path.join(export_dir, "proteins.bs"),
                os.path.join(export_dir, "prot_matches.bs"),
                os.path.join(export_dir, "struct_matches.bs"),
                os.path.join(export_dir, "proteomes.bs"),
                config["meta"]["release"],
                config["meta"]["release_date"],
            ),
            scheduler=dict(queue=queue, mem=4000),
            requires=[
                "insert-entries", "export-proteins", "export-matches",
                "export-structures", "export-proteomes",
                "insert-proteomes", "insert-structures"
            ]
        ),

        # Create EBI Search index
        Task(
            name="ebisearch",
            fn=ebisearch.create_index,
            args=(
                my_ipro_stg,
                os.path.join(export_dir, "proteins.bs"),
                os.path.join(export_dir, "prot_matches.bs"),
                os.path.join(export_dir, "struct_matches.bs"),
                os.path.join(export_dir, "proteomes.bs"),
                config["meta"]["name"],
                config["meta"]["release"],
                config["meta"]["release_date"],
                config["ebisearch"]["dir"]
            ),
            kwargs=dict(writers=3),
            scheduler=dict(queue=queue, mem=32000, tmp=15000, cpu=4),
            requires=[
                "insert-entries", "export-proteins", "export-matches",
                "export-structures", "export-proteomes",
            ],
        ),

        # Indexing Elastic documents
        Task(
            name="init-es-dir",
            fn=elastic.init_dir,
            args=(elastic_dir,),
            scheduler=dict(queue=queue),
        ),
        Task(
            name="create-documents",
            fn=elastic.create_documents,
            args=(
                ora_ipro,
                my_ipro_stg,
                os.path.join(export_dir, "proteins.bs"),
                os.path.join(export_dir, "descriptions.bs"),
                os.path.join(export_dir, "comments.bs"),
                os.path.join(export_dir, "proteomes.bs"),
                os.path.join(export_dir, "prot_matches.bs"),
                elastic_dir
            ),
            kwargs=dict(producers=3, threshold=threshold),
            scheduler=dict(queue=queue, cpu=5, mem=48000),
            requires=(
                "insert-entries", "insert-sets", "insert-taxa",
                "insert-proteomes",
                "export-proteins", "export-descriptions",
                "export-comments", "export-proteomes",
                "export-matches", "init-es-dir"
            )
        )
    ]

    index_tasks = []
    for i, host in enumerate(elastic_hosts):
        task_name = "index-{}".format(i+1)
        index_tasks.append(task_name)
        tasks.append(
            Task(
                name=task_name,
                fn=elastic.index_documents,
                args=(
                    my_ipro_stg,
                    host,
                    config["elastic"]["type"],
                    config["elastic"]["properties"],
                    elastic_dir
                ),
                kwargs=dict(
                    indices=config["elastic"]["indices"],
                    suffix=config["meta"]["release"],
                    shards=config.getint("elastic", "shards"),
                    loaders=4
                ),
                scheduler=dict(queue=queue, cpu=5, mem=4000),
                requires=("init-es-dir", "insert-databases")
            )
        )

    tasks.append(
        Task(
            name="update-alias",
            fn=elastic.update_alias,
            args=(my_ipro_stg, elastic_hosts, config["elastic"]["alias"]),
            kwargs=dict(
                suffix=config["meta"]["release"],
                delete=True
            ),
            scheduler=dict(queue=queue),
            requires=index_tasks
        )
    )

    task_names = []
    for t in tasks:
        if t.name not in task_names:
            task_names.append(t.name)

    if args.tasks is None:
        # Run all tasks
        args.tasks = task_names
    elif args.tasks:
        for arg in args.tasks:
            if arg not in task_names:
                parser.error(
                    "argument -t/--tasks: "
                    "invalid choice: '{}' (choose from {})\n".format(
                        arg,
                        ", ".join(map("'{}'".format, task_names))
                    )
                )
    else:
        args.tasks = []

    wdir = config["workflow"]["dir"]
    wname = "InterPro7 DW"
    with Workflow(tasks, name=wname, dir=wdir, daemon=args.daemon) as w:
        w.run(
            args.tasks,
            secs=0 if args.detach else 10,
            resume=args.resume,
            dry=args.dry_run,
            resubmit=1 if args.retry else 0
        )
