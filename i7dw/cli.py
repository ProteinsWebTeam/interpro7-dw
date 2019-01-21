#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import configparser
import os
import re

from mundone import Task, Workflow

from . import (
    __version__,
    ebisearch,
    interpro,
    uniprot
)


def parse_emails(emails: list, server: dict):
    p = re.compile(r"[a-z0-9]+[a-z0-9\-_.]*@[a-z0-9\-_.]+\.[a-z]{2,}$", re.I)
    for i, email in enumerate(emails):
        if p.match(email) is None:
            raise ValueError("'{}': invalid email address".format(email))
        elif not i:
            server.update({
                "user": email,
                "to": [email]
            })
        else:
            server["to"].append(email)


def parse_server(s: str) -> dict:
    m = re.match(r"([a-z0-9]+[a-z0-9\-_.]*[a-z]{2,})(:\d+)?$", s, re.I)
    if m is None:
        raise ValueError("'{}': invalid server address".format(s))

    host, port = m.groups()
    if port is None:
        return {"host": host}
    else:
        return {"host": host, "port": int(port[1:])}


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
                        metavar="TASK",
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
                        action="store_const",
                        const=0,
                        default=10,
                        help="enqueue tasks to run and quit")
    parser.add_argument("--retry",
                        action="store_const",
                        const=1,
                        default=0,
                        help="rerun failed tasks (once)")
    parser.add_argument("--daemon",
                        action="store_true",
                        default=False,
                        help="do not kill running tasks "
                             "when exiting the program")
    parser.add_argument("--smtp-server",
                        metavar="SERVER:PORT",
                        help="SMTP server for mail notification "
                             "(format: host[:port])")
    parser.add_argument("--send-mail",
                        nargs="+",
                        metavar="ADDRESS",
                        help="recipients' addresses to send an email to "
                             " on completion")
    parser.add_argument("-v", "--version", action="version",
                        version="%(prog)s {}".format(__version__),
                        help="show the version and quit")
    args = parser.parse_args()

    if not os.path.isfile(args.config):
        raise RuntimeError(
            "cannot open '{}': "
            "no such file or directory".format(args.config)
        )

    if args.send_mail:
        if not args.smtp_server:
            parser.error("argument --smtp-server required with --send-mail")

        notif = parse_server(args.smtp_server)
        parse_emails(args.send_mail, notif)
    else:
        notif = None

    config = configparser.ConfigParser()
    config.read(args.config)

    export_dir = config["export"]["dir"]
    os.makedirs(export_dir, exist_ok=True)

    ora_ipro = config["databases"]["interpro_oracle"]
    my_ipro_stg = config["databases"]["interpro_mysql_stg"]
    my_ipro_rel = config["databases"]["interpro_mysql_rel"]
    ora_pdbe = config["databases"]["pdbe_oracle"]
    my_pfam = config["databases"]["pfam_mysql"]
    queue = config["workflow"]["queue"]

    es_clusters = config["elastic"]["clusters"].split(';')
    es_dir = config["elastic"]["dir"]

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
                os.path.join(export_dir, "proteins.dat")
            ),
            kwargs=dict(processes=4),
            scheduler=dict(queue=queue, mem=2000, tmp=3000, cpu=4),
            requires=["chunk-proteins"]
        ),
        Task(
            name="export-sequences",
            fn=interpro.export_sequences,
            args=(
                ora_ipro,
                os.path.join(export_dir, "chunks.json"),
                os.path.join(export_dir, "sequences.dat")
            ),
            kwargs=dict(processes=4),
            scheduler=dict(queue=queue, mem=2000, tmp=30000, cpu=4),
            requires=["chunk-proteins"]
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

        # Load data to MySQL
        Task(
            name="init-tables",
            fn=interpro.init_tables,
            args=(my_ipro_stg,),
            scheduler=dict(queue=queue)
        ),
        Task(
            name="insert-taxa",
            fn=interpro.insert_taxa,
            args=(ora_ipro, my_ipro_stg),
            scheduler=dict(queue=queue, mem=4000),
            requires=["init-tables"]
        ),
        Task(
            name="insert-proteomes",
            fn=interpro.insert_proteomes,
            args=(ora_ipro, my_ipro_stg),
            scheduler=dict(queue=queue, mem=2000),
            requires=["insert-taxa"]
        ),
        Task(
            name="insert-databases",
            fn=interpro.insert_databases,
            args=(ora_ipro, my_ipro_stg, config["meta"]["release"],
                  config["meta"]["release_date"],),
            scheduler=dict(queue=queue),
            requires=["init-tables"]
        ),
        Task(
            name="insert-entries",
            fn=interpro.insert_entries,
            args=(ora_ipro, my_pfam, my_ipro_stg),
            scheduler=dict(queue=queue, mem=2000),
            requires=["insert-databases"]
        ),
        Task(
            name="insert-annotations",
            fn=interpro.insert_annotations,
            args=(my_pfam, my_ipro_stg),
            scheduler=dict(queue=queue, mem=4000),
            requires=["insert-entries"]
        ),
        Task(
            name="insert-structures",
            fn=interpro.insert_structures,
            args=(ora_ipro, my_ipro_stg),
            scheduler=dict(queue=queue, mem=1000),
            requires=["insert-databases"]
        ),
        Task(
            name="insert-sets",
            fn=interpro.insert_sets,
            args=(ora_ipro, my_pfam, my_ipro_stg),
            scheduler=dict(queue=queue, mem=3000, tmp=3000),
            requires=["insert-entries"]
        ),
        Task(
            name="insert-proteins",
            fn=interpro.insert_proteins,
            args=(
                ora_ipro,
                ora_pdbe,
                my_ipro_stg,
                os.path.join(export_dir, "proteins.dat"),
                os.path.join(export_dir, "sequences.dat"),
                os.path.join(export_dir, "misc.dat"),
                os.path.join(export_dir, "names.dat"),
                os.path.join(export_dir, "comments.dat"),
                os.path.join(export_dir, "proteomes.dat"),
                os.path.join(export_dir, "residues.dat"),
                os.path.join(export_dir, "features.dat"),
                os.path.join(export_dir, "matches.dat")
            ),
            scheduler=dict(queue=queue, mem=24000),
            requires=[
                "insert-entries", "insert-structures", "insert-taxa",
                "insert-sets", "export-proteins", "export-sequences",
                "export-misc", "export-names", "export-comments",
                "export-proteomes", "export-residues", "export-features",
                "export-matches"
            ]
        ),
        Task(
            name="release-notes",
            fn=interpro.make_release_notes,
            args=(
                my_ipro_stg,
                my_ipro_rel,
                os.path.join(export_dir, "proteins.dat"),
                os.path.join(export_dir, "matches.dat"),
                os.path.join(export_dir, "proteomes.dat"),
                config["meta"]["release"],
                config["meta"]["release_date"],
            ),
            scheduler=dict(queue=queue, mem=4000),
            requires=[
                "insert-entries", "insert-proteomes", "insert-structures",
                # insert-annotations only so it's not forgiven
                "insert-annotations",
                "export-proteins", "export-matches", "export-proteomes"
            ]
        ),

        # Overlapping homologous superfamilies
        Task(
            name="overlapping-families",
            fn=interpro.calculate_relationships,
            args=(
                my_ipro_stg,
                os.path.join(export_dir, "proteins.dat"),
                os.path.join(export_dir, "matches.dat"),
                threshold
            ),
            # kwargs=dict(ora_uri=ora_ipro),
            scheduler=dict(queue=queue, mem=6000),
            requires=["export-proteins", "export-matches", "insert-entries"]
        ),

        # Cross-references
        Task(
            name="export-xrefs",
            fn=interpro.export_xrefs,
            args=(
                my_ipro_stg,
                os.path.join(export_dir, "proteins.dat"),
                os.path.join(export_dir, "matches.dat"),
                os.path.join(export_dir, "proteomes.dat"),
                os.path.join(export_dir, "entries_xref.dat"),
                os.path.join(export_dir, "proteomes_xref.dat"),
                os.path.join(export_dir, "structures_xref.dat"),
                os.path.join(export_dir, "taxa_xref.dat")
            ),
            scheduler=dict(queue=queue, mem=24000, tmp=20000, cpu=5),
            requires=[
                "export-proteins", "export-matches", "export-proteomes",
                "insert-entries", "insert-proteomes", "insert-sets",
                "insert-structures"
            ]
        ),

        Task(
            name="update-counts",
            fn=interpro.update_counts,
            args=(
                my_ipro_stg,
                os.path.join(export_dir, "entries_xref.dat"),
                os.path.join(export_dir, "proteomes_xref.dat"),
                os.path.join(export_dir, "structures_xref.dat"),
                os.path.join(export_dir, "taxa_xref.dat")
            ),
            scheduler=dict(queue=queue, mem=16000, tmp=15000),
            requires=[
                "export-xrefs", "insert-proteins", "overlapping-families"
            ]
        ),

        # Create EBI Search index
        Task(
            name="ebi-search",
            fn=ebisearch.dump,
            args=(
                my_ipro_stg,
                os.path.join(export_dir, "entries_xref.dat"),
                config["meta"]["name"],
                config["meta"]["release"],
                config["meta"]["release_date"],
                config["ebisearch"]["dir"]
            ),
            scheduler=dict(queue=queue, mem=16000, cpu=4),
            requires=["export-xrefs"],
        ),

        # Indexing Elastic documents
        Task(
            name="init-elastic",
            fn=interpro.init_elastic,
            args=(es_dir,),
            scheduler=dict(queue=queue),
            requires=[
                "insert-entries", "insert-sets", "insert-proteomes",
                "export-proteins", "export-names", "export-comments",
                "export-proteomes", "export-matches"
            ]
        ),
        Task(
            name="create-documents",
            fn=interpro.create_documents,
            args=(
                ora_ipro,
                my_ipro_stg,
                os.path.join(export_dir, "proteins.dat"),
                os.path.join(export_dir, "names.dat"),
                os.path.join(export_dir, "comments.dat"),
                os.path.join(export_dir, "proteomes.dat"),
                os.path.join(export_dir, "matches.dat"),
                es_dir
            ),
            kwargs=dict(processes=8),
            scheduler=dict(queue=queue, cpu=8, mem=32000),
            requires=["init-elastic"]
        )
    ]

    for i, hosts in enumerate(es_clusters):
        hosts = list(set(hosts.split(',')))
        dst = es_dir.rstrip('/') + "-" + str(i+1)

        t = Task(
            name="index-" + str(i+1),
            fn=interpro.index_documents,
            args=(
                my_ipro_stg,
                hosts,
                config["elastic"]["type"],
                es_dir
            ),
            kwargs=dict(
                body=config["elastic"]["body"],
                create_indices=True,
                custom_shards=config["elastic"]["indices"],
                default_shards=config.getint("elastic", "shards"),
                processes=4,
                raise_on_error=False,
                suffix=config["meta"]["release"],
                outdir=dst
            ),
            scheduler=dict(queue=queue, cpu=4, mem=8000),
            requires=["init-elastic"]
        )

        tasks.append(t)

        tasks.append(
            Task(
                name="complete-index-{}".format(i + 1),
                fn=interpro.index_documents,
                args=(
                    my_ipro_stg,
                    hosts,
                    config["elastic"]["type"],
                    dst
                ),
                kwargs=dict(
                    alias="staging",
                    files=t.output,
                    max_retries=5,
                    processes=4,
                    suffix=config["meta"]["release"],
                    writeback=True
                ),
                scheduler=dict(queue=queue, cpu=4, mem=8000),
                requires=["index-{}".format(i + 1), "create-documents"]
            )
        )

        tasks.append(
            Task(
                name="update-alias-{}".format(i+1),
                fn=interpro.update_alias,
                args=(my_ipro_stg, hosts, config["elastic"]["alias"]),
                kwargs=dict(
                    suffix=config["meta"]["release"],
                    delete=True
                ),
                scheduler=dict(queue=queue),
                requires=["complete-index-{}".format(i + 1)]
            )
        )

    task_names = []
    for t in tasks:
        if t.name not in task_names:
            task_names.append(t.name)

    if args.tasks:
        for arg in args.tasks:
            if arg not in task_names:
                parser.error(
                    "argument -t/--tasks: "
                    "invalid choice: '{}' (choose from {})\n".format(
                        arg,
                        ", ".join(map("'{}'".format, task_names))
                    )
                )

    wdir = config["workflow"]["dir"]
    wname = "InterPro7 DW"
    with Workflow(tasks, name=wname, dir=wdir, daemon=args.daemon,
                  mail=notif) as w:
        w.run(
            args.tasks,
            secs=args.detach,
            resume=args.resume,
            dry=args.dry_run,
            resubmit=args.retry
        )
