#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import configparser
import os
import sys

from mundone import Task, Workflow

from i7dw import elastic, mysql
from i7dw.ebi import ebisearch, goa, interpro, uniprot


def parse_elastic_hosts(str_hosts):
    hosts = []
    for host in str_hosts.split(','):
        host = host.strip()

        if not host:
            continue

        pair = host.split(':')
        if len(pair) == 2:
            host = pair[0]
            port = int(pair[1])
        else:
            host = host
            port = 9200

        hosts.append({'host': host, 'port': port})

    return hosts


def cli():
    parser = argparse.ArgumentParser(
        description="Build InterPro7 data warehouse"
    )
    parser.add_argument("config",
                        metavar="config.ini",
                        help="configuration file")
    parser.add_argument("-t", "--tasks",
                        nargs="*",
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
                        help="collect completed tasks, "
                             "submit tasks ready, and quit")
    args = parser.parse_args()

    if not os.path.isfile(args.config):
        raise RuntimeError(
            "cannot open '{}': "
            "no such file or directory".format(args.config)
        )

    config = configparser.ConfigParser()
    config.read(args.config)

    export_dir = config['export']['dir']
    if not os.path.isdir(export_dir):
        os.makedirs(export_dir)

    ora_ipro = config['databases']['interpro_oracle']
    my_ipro = config['databases']['interpro_mysql']
    my_pfam = config['databases']['pfam_mysql']
    queue = config['workflow']['queue']

    elastic_hosts = config['elastic']['hosts'].split(',')
    elastic_dir = config['elastic']['dir']

    threshold = config.getfloat("jaccard", "threshold")

    tasks = [
        # Export data to stores
        Task(
            name="export_comments",
            fn=uniprot.export_protein_comments,
            args=(ora_ipro, os.path.join(export_dir, "comments.bs")),
            kwargs=dict(chunk_size=100000),
            scheduler=dict(queue=queue, mem=1000)
        ),
        Task(
            name="export_descriptions",
            fn=uniprot.export_protein_descriptions,
            args=(ora_ipro, os.path.join(export_dir, 'descriptions.bs')),
            kwargs=dict(chunk_size=100000),
            scheduler=dict(queue=queue, mem=1000)
        ),
        Task(
            name="export_evidences",
            fn=uniprot.export_protein_evidence,
            args=(ora_ipro, os.path.join(export_dir, 'evidences.bs')),
            kwargs=dict(chunk_size=100000),
            scheduler=dict(queue=queue, mem=500)
        ),
        Task(
            name="export_gene",
            fn=uniprot.export_protein_gene,
            args=(ora_ipro, os.path.join(export_dir, 'genes.bs')),
            kwargs=dict(chunk_size=100000),
            scheduler=dict(queue=queue, mem=500)
        ),
        Task(
            name="export_proteomes",
            fn=uniprot.export_protein_proteomes,
            args=(ora_ipro, os.path.join(export_dir, 'proteomes.bs')),
            kwargs=dict(chunk_size=100000),
            scheduler=dict(queue=queue, mem=500)
        ),
        Task(
            name="export_annotations",
            fn=goa.export_annotations,
            args=(ora_ipro, os.path.join(export_dir, 'annotations.bs')),
            kwargs=dict(chunk_size=100000),
            scheduler=dict(queue=queue, mem=2000)
        ),
        Task(
            name="export_structures",
            fn=interpro.export_struct_matches,
            args=(ora_ipro, os.path.join(export_dir, 'struct_matches.bs')),
            kwargs=dict(chunk_size=100000),
            scheduler=dict(queue=queue, mem=2000)
        ),
        Task(
            name="export_matches",
            fn=interpro.export_prot_matches,
            args=(ora_ipro, os.path.join(export_dir, 'prot_matches.bs')),
            kwargs=dict(chunk_size=100000),
            scheduler=dict(queue=queue, mem=3000)
        ),
        Task(
            name="export_features",
            fn=interpro.export_prot_matches_extra,
            args=(ora_ipro, os.path.join(export_dir, 'prot_matches_extra.bs')),
            kwargs=dict(chunk_size=100000),
            scheduler=dict(queue=queue, mem=2000)
        ),
        Task(
            fn=interpro.export_residues,
            args=(ora_ipro, os.path.join(export_dir, 'residues.bs')),
            kwargs=dict(chunk_size=100000),
            scheduler=dict(queue=queue, mem=3000)
        ),
        Task(
            fn=interpro.export_proteins,
            args=(ora_ipro, os.path.join(export_dir, 'proteins.bs')),
            kwargs=dict(chunk_size=100000),
            scheduler=dict(queue=queue, mem=2000)
        ),

        # Load data to MySQL
        Task(
            fn=mysql.init_tables,
            args=(my_ipro,),
            scheduler=dict(queue=queue)
        ),
        Task(
            fn=mysql.insert_taxa,
            args=(ora_ipro, my_ipro),
            scheduler=dict(queue=queue, mem=4000),
            requires=["init_tables"]
        ),
        Task(
            fn=mysql.insert_proteomes,
            args=(ora_ipro, my_ipro),
            scheduler=dict(queue=queue, mem=1000),
            requires=["insert_taxa"]
        ),
        Task(
            fn=mysql.insert_databases,
            args=(ora_ipro, my_ipro),
            scheduler=dict(queue=queue),
            requires=["init_tables"]
        ),
        Task(
            fn=mysql.insert_entries,
            args=(ora_ipro, my_pfam, my_ipro),
            scheduler=dict(queue=queue, mem=2000),
            requires=["insert_databases"]
        ),
        Task(
            fn=mysql.insert_annotations,
            args=(my_pfam, my_ipro),
            scheduler=dict(queue=queue, mem=4000),
            requires=["insert_entries"]
        ),
        Task(
            fn=mysql.insert_structures,
            args=(ora_ipro, my_ipro),
            scheduler=dict(queue=queue, mem=1000),
            requires=["insert_databases"]
        ),
        Task(
            fn=mysql.insert_sets,
            args=(my_pfam, my_ipro),
            scheduler=dict(queue=queue),
            requires=["insert_databases"]
        ),
        Task(
            fn=mysql.insert_proteins,
            args=(
                my_ipro,
                os.path.join(export_dir, "proteins.bs"),
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
                "insert_taxa", "export_proteins", "export_evidences",
                "export_descriptions", "export_comments",
                "export_proteomes", "export_genes",
                "export_annotations", "export_residues",
                "export_structures", "export_features"
            ]
        ),

        # Create EBI Search index
        Task(
            fn=ebisearch.export,
            args=(
                my_ipro,
                os.path.join(export_dir, "proteins.bs"),
                os.path.join(export_dir, 'prot_matches.bs'),
                os.path.join(export_dir, "struct_matches.bs"),
                os.path.join(export_dir, "proteomes.bs"),
                config["meta"]["name"],
                config["meta"]["release"],
                config["meta"]["release_date"],
                config["ebisearch"]["dir"]
            ),
            scheduler=dict(
                queue=config['workflow']['queue'], mem=24000, tmp=10000
            ),
            requires=[
                "insert_entries", "export_proteins", "export_matches",
                "export_structures", "export_proteomes",
            ],
        ),

        # Indexing Elastic documents
        Task(
            fn=elastic.init_dir,
            args=(elastic_dir,),
            scheduler=dict(queue=queue),
        ),
        Task(
            fn=elastic.create_documents,
            args=(
                ora_ipro,
                my_ipro,
                os.path.join(export_dir, "proteins.bs"),
                os.path.join(export_dir, "descriptions.bs"),
                os.path.join(export_dir, "comments.bs"),
                os.path.join(export_dir, "proteomes.bs"),
                os.path.join(export_dir, 'prot_matches.bs'),
                elastic_dir
            ),
            kwargs=dict(producers=3, threshold=threshold),
            scheduler=dict(queue=queue, cpu=4, mem=64000),
            requires=(
                "insert_entries", "insert_sets", "insert_taxa",
                "export_proteins", "export_descriptions",
                "export_comments", "export_proteomes",
                "export_matches", "init_dir"
            )
        )
    ]

    for i, host in enumerate(elastic_hosts):
        tasks.append(
            Task(
                name="index-{}".format(i+1),
                fn=elastic.index_documents,
                args=(
                    my_ipro,
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
                scheduler=dict(queue=queue, cpu=5, mem=16000),
                requires=("init_dir", "insert_databases")
            )
        )

    task_names = []
    for t in tasks:
        if t.name not in task_names:
            task_names.append(t.name)

    if args.tasks:
        for arg in args.tasks:
            if arg not in task_names:
                sys.stderr.write(
                    "{}: error: "
                    "argument -t/--tasks: "
                    "invalid choice: '{}' (choose from {})\n".format(
                        os.path.basename(sys.argv[0]),
                        arg,
                        ', '.join(map("'{}'".format, task_names))
                    )
                )
                exit(1)
    else:
        args.tasks = []

    w = Workflow(tasks, name="InterPro7 DW", dir=config["workflow"]["dir"])
    w.run(
        args.tasks,
        secs=0 if args.detach else 60,
        resume=args.resume,
        dry=args.dry_run
    )
