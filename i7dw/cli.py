#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import configparser
import os

from mundone import Task, Workflow

from i7dw import __version__, elastic, mysql, xref
from i7dw.ebi import ebisearch, interpro, uniprot


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
    os.makedirs(export_dir, exist_ok=True)

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
            name="chunk-proteins",
            fn=interpro.chunk_proteins,
            args=(ora_ipro, os.path.join(export_dir, "chunks.json")),
            kwargs=dict(order_by=False),
            scheduler=dict(queue=queue, mem=12000),
        ),
        Task(
            name="export-comments",
            fn=uniprot.export_protein2comments,
            args=(
                ora_ipro,
                os.path.join(export_dir, "chunks.json"),
                os.path.join(export_dir, "comments.dat")
            ),
            kwargs=dict(processes=4),
            scheduler=dict(queue=queue, mem=1000, cpu=4),
            requires=["chunk-proteins"]
        ),
        Task(
            name="export-names",
            fn=uniprot.export_protein2names,
            args=(
                ora_ipro,
                os.path.join(export_dir, "chunks.json"),
                os.path.join(export_dir, "names.dat")
            ),
            kwargs=dict(processes=4),
            scheduler=dict(queue=queue, mem=2000, tmp=3000, cpu=4),
            requires=["chunk-proteins"]
        ),
        Task(
            name="export-misc",
            fn=uniprot.export_protein2supplementary,
            args=(
                ora_ipro,
                os.path.join(export_dir, "chunks.json"),
                os.path.join(export_dir, "misc.dat")
            ),
            kwargs=dict(processes=4),
            scheduler=dict(queue=queue, mem=1000, tmp=1000, cpu=4),
            requires=["chunk-proteins"]
        ),
        Task(
            name="export-proteomes",
            fn=uniprot.export_protein2proteome,
            args=(
                ora_ipro,
                os.path.join(export_dir, "chunks.json"),
                os.path.join(export_dir, "proteomes.dat")
            ),
            kwargs=dict(processes=4),
            scheduler=dict(queue=queue, mem=500, tmp=500, cpu=4),
            requires=["chunk-proteins"]
        ),
        Task(
            name="export-structures",
            fn=interpro.export_protein2structures,
            args=(
                ora_ipro,
                os.path.join(export_dir, "chunks.json"),
                os.path.join(export_dir, "structures.dat")
            ),
            kwargs=dict(processes=4),
            scheduler=dict(queue=queue, mem=4000, tmp=200, cpu=4),
            requires=["chunk-proteins"]
        ),
        Task(
            name="export-matches",
            fn=interpro.export_protein2matches,
            args=(
                ora_ipro,
                os.path.join(export_dir, "chunks.json"),
                os.path.join(export_dir, "matches.dat")
            ),
            kwargs=dict(processes=4),
            scheduler=dict(queue=queue, mem=8000, tmp=20000, cpu=4),
            requires=["chunk-proteins"]
        ),
        Task(
            name="export-features",
            fn=interpro.export_protein2features,
            args=(
                ora_ipro,
                os.path.join(export_dir, "chunks.json"),
                os.path.join(export_dir, "features.dat")
            ),
            kwargs=dict(processes=4),
            scheduler=dict(queue=queue, mem=2000, tmp=8000, cpu=4),
            requires=["chunk-proteins"]
        ),
        Task(
            name="export-residues",
            fn=interpro.export_protein2residues,
            args=(
                ora_ipro,
                os.path.join(export_dir, "chunks.json"),
                os.path.join(export_dir, "residues.dat")
            ),
            kwargs=dict(processes=4),
            scheduler=dict(queue=queue, mem=3000, tmp=8000, cpu=4),
            requires=["chunk-proteins"]
        ),
        Task(
            name="export-proteins",
            fn=interpro.export_proteins,
            args=(
                ora_ipro,
                os.path.join(export_dir, "chunks.json"),
                os.path.join(export_dir, "proteins.dat"),
                os.path.join(export_dir, "sequences.dat")
            ),
            kwargs=dict(processes=4),
            scheduler=dict(queue=queue, mem=2000, tmp=30000, cpu=4),
            requires=["chunk-proteins"]
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
            args=(ora_ipro, my_pfam, my_ipro_stg),
            scheduler=dict(queue=queue, mem=1000),
            requires=["insert-databases"]
        ),
        Task(
            name="insert-proteins",
            fn=mysql.insert_proteins,
            args=(
                my_ipro_stg,
                os.path.join(export_dir, "proteins.dat"),
                os.path.join(export_dir, "sequences.dat"),
                os.path.join(export_dir, "misc.dat"),
                os.path.join(export_dir, "names.dat"),
                os.path.join(export_dir, "comments.dat"),
                os.path.join(export_dir, "proteomes.dat"),
                os.path.join(export_dir, "residues.dat"),
                os.path.join(export_dir, "structures.dat"),
                os.path.join(export_dir, "features.dat"),
                os.path.join(export_dir, "matches.dat")
            ),
            scheduler=dict(queue=queue, mem=24000),
            requires=[
                "insert-entries", "insert-structures", "insert-taxa",
                "insert-proteomes", "export-proteins", "export-misc",
                "export-names", "export-comments", "export-proteomes",
                "export-residues", "export-structures", "export-features",
                "export-matches"
            ]
        ),

        # Cross-references (counts in MySQL)
        Task(
            name="export-xrefs",
            fn=xref.export,
            args=(
                my_ipro_stg,
                os.path.join(export_dir, "proteins.dat"),
                os.path.join(export_dir, "matches.dat"),
                os.path.join(export_dir, "proteomes.dat"),
                os.path.join(export_dir, "entries_xref.dat"),
                os.path.join(export_dir, "taxa_xref.dat"),
                os.path.join(export_dir, "proteomes_xref.dat"),
                os.path.join(export_dir, "sets_xref.dat"),
                os.path.join(export_dir, "structures_xref.dat")
            ),
            scheduler=dict(queue=queue, mem=24000, tmp=16000, cpu=6),
            requires=[
                "insert-entries", "insert-taxa", "insert-proteomes",
                "insert-sets", "insert-structures", "export-proteins",
                "export-matches", "export-proteomes"
            ]
        ),
        Task(
            name="update-counts",
            fn=mysql.update_counts,
            args=(
                my_ipro_stg,
                os.path.join(export_dir, "entries_xref.dat"),
                os.path.join(export_dir, "taxa_xref.dat"),
                os.path.join(export_dir, "proteomes_xref.dat"),
                os.path.join(export_dir, "sets_xref.dat"),
                os.path.join(export_dir, "structures_xref.dat")
            ),
            scheduler=dict(queue=queue, mem=24000),
            requires=["insert-proteins", "export-xrefs"]
        ),

        # Release notes
        Task(
            name="release-notes",
            fn=mysql.make_release_notes,
            args=(
                my_ipro_stg,
                my_ipro_rel,
                os.path.join(export_dir, "proteins.dat"),
                os.path.join(export_dir, "matches.dat"),
                os.path.join(export_dir, "structures.dat"),
                os.path.join(export_dir, "proteomes.dat"),
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
                os.path.join(export_dir, "proteins.dat"),
                os.path.join(export_dir, "matches.dat"),
                os.path.join(export_dir, "structures.dat.bs"),
                os.path.join(export_dir, "proteomes.dat"),
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
                os.path.join(export_dir, "proteins.dat"),
                os.path.join(export_dir, "descriptions.bs"),
                os.path.join(export_dir, "comments.dat"),
                os.path.join(export_dir, "proteomes.dat"),
                os.path.join(export_dir, "matches.dat"),
                elastic_dir
            ),
            kwargs=dict(producers=3, threshold=threshold),
            scheduler=dict(queue=queue, cpu=5, mem=48000),
            requires=(
                "insert-entries", "insert-sets", "insert-taxa",
                "insert-proteomes",
                "export-proteins", "export-names",
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
