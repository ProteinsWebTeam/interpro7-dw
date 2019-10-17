#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import configparser
import os
import re
import sys

from mundone import Task, Workflow

from i7dw import __version__, ebisearch, goa, uniprot
from i7dw.interpro import elastic, mysql, oracle


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


def build_dw():
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
                        default=False,
                        help="skip completed tasks")
    parser.add_argument("--detach",
                        action="store_const",
                        const=0,
                        default=10,
                        help="enqueue tasks to run and exit")
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
                        help="show the version and exit")
    args = parser.parse_args()

    if not os.path.isfile(args.config):
        parser.error(f"cannot open '{args.config}': no such file or directory")

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

    ora_ipro = config["databases"]["interpro_production"]
    my_ipro_stg = config["databases"]["interpro_staging"]
    my_ipro_rel = config["databases"]["interpro_offsite"]
    ora_pdbe = config["databases"]["pdbe"]
    my_pfam = config["databases"]["pfam"]
    queue = config["workflow"]["queue"]

    es_clusters = config["elastic"]["clusters"].split(';')
    es_dir = config["elastic"]["dir"]

    tasks = [
        Task(
            name="chunk-proteins",
            fn=oracle.export.chunk_proteins,
            args=(ora_ipro, os.path.join(export_dir, "chunks.json")),
            scheduler=dict(queue=queue, mem=16000),
        ),
        Task(
            name="export-matches",
            fn=oracle.export.export_matches,
            args=(
                ora_ipro,
                os.path.join(export_dir, "chunks.json"),
                os.path.join(export_dir, "matches.dat")
            ),
            kwargs=dict(processes=4, tmpdir="/scratch"),
            scheduler=dict(queue=queue, mem=8000, scratch=25000, cpu=4),
            requires=["chunk-proteins"]
        ),
        Task(
            name="export-features",
            fn=oracle.export.export_features,
            args=(
                ora_ipro,
                os.path.join(export_dir, "chunks.json"),
                os.path.join(export_dir, "features.dat")
            ),
            kwargs=dict(processes=4, tmpdir="/scratch"),
            scheduler=dict(queue=queue, mem=4000, scratch=10000, cpu=4),
            requires=["chunk-proteins"]
        ),
        Task(
            name="export-residues",
            fn=oracle.export.export_residues,
            args=(
                ora_ipro,
                os.path.join(export_dir, "chunks.json"),
                os.path.join(export_dir, "residues.dat")
            ),
            kwargs=dict(processes=4, tmpdir="/scratch"),
            scheduler=dict(queue=queue, mem=3000, scratch=10000, cpu=4),
            requires=["chunk-proteins"]
        ),
        Task(
            name="export-proteins",
            fn=oracle.export.export_proteins,
            args=(
                ora_ipro,
                os.path.join(export_dir, "chunks.json"),
                os.path.join(export_dir, "proteins.dat")
            ),
            kwargs=dict(processes=4, tmpdir="/scratch"),
            scheduler=dict(queue=queue, mem=2000, scratch=3000, cpu=4),
            requires=["chunk-proteins"]
        ),
        Task(
            name="export-sequences",
            fn=oracle.export.export_sequences,
            args=(
                ora_ipro,
                os.path.join(export_dir, "chunks.json"),
                os.path.join(export_dir, "sequences.dat")
            ),
            kwargs=dict(processes=4, tmpdir="/scratch"),
            scheduler=dict(queue=queue, mem=4000, scratch=35000, cpu=4),
            requires=["chunk-proteins"]
        ),
        Task(
            name="export-comments",
            fn=uniprot.export_comments,
            args=(
                ora_ipro,
                os.path.join(export_dir, "chunks.json"),
                os.path.join(export_dir, "comments.dat")
            ),
            kwargs=dict(processes=4),
            scheduler=dict(queue=queue, mem=1000, scratch=1000, cpu=4),
            requires=["chunk-proteins"]
        ),
        Task(
            name="export-names",
            fn=uniprot.export_descriptions,
            args=(
                ora_ipro,
                os.path.join(export_dir, "chunks.json"),
                os.path.join(export_dir, "names.dat")
            ),
            kwargs=dict(processes=4, tmpdir="/scratch"),
            scheduler=dict(queue=queue, mem=2000, scratch=4000, cpu=4),
            requires=["chunk-proteins"]
        ),
        Task(
            name="export-misc",
            fn=uniprot.export_misc,
            args=(
                ora_ipro,
                os.path.join(export_dir, "chunks.json"),
                os.path.join(export_dir, "misc.dat")
            ),
            kwargs=dict(processes=4, tmpdir="/scratch"),
            scheduler=dict(queue=queue, mem=1000, scratch=2000, cpu=4),
            requires=["chunk-proteins"]
        ),
        Task(
            name="export-proteomes",
            fn=uniprot.export_proteomes,
            args=(
                ora_ipro,
                os.path.join(export_dir, "chunks.json"),
                os.path.join(export_dir, "proteomes.dat")
            ),
            kwargs=dict(processes=4, tmpdir="/scratch"),
            scheduler=dict(queue=queue, mem=1000, scratch=1000, cpu=4),
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
            name="insert-databases",
            fn=mysql.databases.insert_databases,
            args=(my_ipro_stg, ora_ipro, config["meta"]["release"],
                  config["meta"]["release_date"],),
            scheduler=dict(queue=queue),
            requires=["init-tables"]
        ),
        Task(
            name="insert-taxa",
            fn=mysql.taxonomy.insert_taxa,
            args=(my_ipro_stg, ora_ipro),
            scheduler=dict(queue=queue, mem=4000),
            requires=["init-tables"]
        ),
        Task(
            name="insert-proteomes",
            fn=mysql.proteomes.insert_proteomes,
            args=(my_ipro_stg, ora_ipro),
            scheduler=dict(queue=queue, mem=2000),
            requires=["insert-taxa"]
        ),
        Task(
            name="insert-entries",
            fn=mysql.entries.insert_entries,
            args=(my_ipro_stg, ora_ipro, my_pfam),
            scheduler=dict(queue=queue, mem=2000),
            requires=["insert-databases"]
        ),
        Task(
            name="insert-annotations",
            fn=mysql.entries.insert_annotations,
            args=(my_ipro_stg, my_pfam),
            scheduler=dict(queue=queue, mem=4000),
            requires=["insert-entries"]
        ),
        Task(
            name="insert-sets",
            fn=mysql.entries.insert_sets,
            args=(my_ipro_stg, ora_ipro, my_pfam),
            scheduler=dict(queue=queue, mem=1000),
            requires=["insert-entries"]
        ),
        Task(
            name="insert-structures",
            fn=mysql.structures.insert_structures,
            args=(my_ipro_stg, ora_ipro),
            scheduler=dict(queue=queue, mem=4000),
            requires=["insert-databases"]
        ),
        Task(
            name="insert-isoforms",
            fn=mysql.proteins.insert_isoforms,
            args=(my_ipro_stg, ora_ipro),
            scheduler=dict(queue=queue, mem=4000),
            requires=["insert-entries"]
        ),
        Task(
            name="insert-proteins",
            fn=mysql.proteins.insert_proteins,
            args=(
                my_ipro_stg,
                ora_ipro,
                ora_pdbe,
                os.path.join(export_dir, "comments.dat"),
                os.path.join(export_dir, "features.dat"),
                os.path.join(export_dir, "matches.dat"),
                os.path.join(export_dir, "misc.dat"),
                os.path.join(export_dir, "names.dat"),
                os.path.join(export_dir, "proteins.dat"),
                os.path.join(export_dir, "proteomes.dat"),
                os.path.join(export_dir, "residues.dat"),
                os.path.join(export_dir, "sequences.dat")
            ),
            kwargs=dict(processes=4),
            scheduler=dict(queue=queue, mem=24000, cpu=4),
            requires=["export-comments", "export-features", "export-matches",
                      "export-misc", "export-names", "export-proteins",
                      "export-proteomes", "export-residues",
                      "export-sequences", "insert-entries", "insert-isoforms",
                      "insert-sets", "insert-structures", "insert-taxa"
                      ]
        ),
        Task(
            name="release-notes",
            fn=mysql.relnotes.make_release_notes,
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
            requires=["export-matches", "export-proteins", "export-proteomes",
                      "insert-databases", "insert-entries", "insert-proteomes",
                      "insert-sets", "insert-structures", "insert-taxa"]
        ),

        # Overlapping entries
        Task(
            name="overlapping-families",
            fn=mysql.entries.find_overlapping_entries,
            args=(
                my_ipro_stg,
                os.path.join(export_dir, "matches.dat")
            ),
            kwargs=dict(ora_url=ora_ipro),
            scheduler=dict(queue=queue, mem=6000),
            requires=["export-matches", "insert-entries"]
        ),

        # Mappings for GOA team
        Task(
            name="export-goa",
            fn=goa.export_mappings,
            args=(my_ipro_stg, ora_ipro, config["export"]["goa"]),
            scheduler=dict(queue=queue, mem=2000),
            requires=["insert-databases"]
        ),

        # Export entries
        Task(
            name="export-entries",
            fn=mysql.entries.export,
            args=(
                my_ipro_stg,
                os.path.join(export_dir, "proteins.dat"),
                os.path.join(export_dir, "proteomes.dat"),
                os.path.join(export_dir, "matches.dat"),
                os.path.join(export_dir, "entries.dat")
            ),
            kwargs=dict(processes=4, tmpdir="/scratch"),
            scheduler=dict(queue=queue, mem=12000, scratch=30000, cpu=4),
            requires=["export-matches", "export-proteins", "export-proteomes",
                      "insert-entries", "insert-sets", "insert-structures"]
        ),
        Task(
            name="update-entries",
            fn=mysql.entries.update_counts,
            args=(
                my_ipro_stg,
                os.path.join(export_dir, "entries.dat")
            ),
            kwargs=dict(tmpdir="/scratch"),
            scheduler=dict(queue=queue, mem=8000, scratch=40000),
            requires=["export-entries", "overlapping-families"]
        ),
        Task(
            name="update-proteomes",
            fn=mysql.proteomes.update_counts,
            args=(
                my_ipro_stg,
                os.path.join(export_dir, "proteins.dat"),
                os.path.join(export_dir, "proteomes.dat"),
                os.path.join(export_dir, "matches.dat"),
            ),
            kwargs=dict(processes=4, tmpdir="/scratch"),
            scheduler=dict(queue=queue, mem=8000, scratch=2000, cpu=4),
            requires=["export-matches", "export-proteins", "export-proteomes",
                      "insert-entries", "insert-sets", "insert-proteomes",
                      "insert-structures"]
        ),
        Task(
            name="update-structures",
            fn=mysql.structures.update_counts,
            args=(
                my_ipro_stg,
                os.path.join(export_dir, "proteins.dat"),
                os.path.join(export_dir, "proteomes.dat"),
                os.path.join(export_dir, "matches.dat")
            ),
            kwargs=dict(processes=4, tmpdir="/scratch"),
            scheduler=dict(queue=queue, mem=8000, scratch=100, cpu=4),
            requires=["export-matches", "export-proteins", "export-proteomes",
                      "insert-entries", "insert-sets", "insert-structures"]
        ),
        Task(
            name="update-taxa",
            fn=mysql.taxonomy.update_counts,
            args=(
                my_ipro_stg,
                os.path.join(export_dir, "proteins.dat"),
                os.path.join(export_dir, "proteomes.dat"),
                os.path.join(export_dir, "matches.dat")
            ),
            kwargs=dict(processes=4, tmpdir="/scratch"),
            scheduler=dict(queue=queue, mem=32000, scratch=30000, cpu=4),
            requires=["export-matches", "export-proteins", "export-proteomes",
                      "insert-entries", "insert-sets", "insert-structures",
                      "insert-taxa"]
        ),

        # Create EBI Search index
        Task(
            name="ebi-search",
            fn=ebisearch.dump,
            args=(
                my_ipro_stg,
                os.path.join(export_dir, "entries.dat"),
                config["meta"]["name"],
                config["meta"]["release"],
                config["meta"]["release_date"],
                config["ebisearch"]["stg"]
            ),
            kwargs=dict(processes=4, by_type=True),
            scheduler=dict(queue=queue, mem=24000, cpu=4),
            requires=["update-entries"],
        ),

        # Replace previous dump by new one
        Task(
            name="publish-ebi-search",
            fn=ebisearch.exchange,
            args=(
                config["ebisearch"]["stg"],
                config["ebisearch"]["rel"]
            ),
            scheduler=dict(queue=queue),
            requires=["ebi-search"],
        ),

        # Indexing Elastic documents
        Task(
            name="init-elastic",
            fn=elastic.init,
            args=(os.path.join(es_dir, "documents"),),
            scheduler=dict(queue=queue),
            requires=["export-comments", "export-matches", "export-names",
                      "export-proteins", "export-proteomes", "insert-entries",
                      "insert-proteomes", "insert-sets", "insert-structures",
                      "insert-taxa"]
        ),
        Task(
            name="create-documents",
            fn=elastic.write_documents,
            args=(
                my_ipro_stg,
                os.path.join(export_dir, "comments.dat"),
                os.path.join(export_dir, "matches.dat"),
                os.path.join(export_dir, "names.dat"),
                os.path.join(export_dir, "proteomes.dat"),
                os.path.join(export_dir, "proteins.dat"),
                os.path.join(es_dir, "documents")
            ),
            kwargs=dict(processes=8),
            scheduler=dict(queue=queue, cpu=8, mem=32000),
            requires=["init-elastic"]
        )
    ]

    for i, hosts in enumerate(es_clusters):
        hosts = list(set(hosts.split(',')))
        outdir = os.path.join(es_dir, "cluster-" + str(i+1))

        tasks += [
            Task(
                name="index-" + str(i+1),
                fn=elastic.index_documents,
                args=(
                    my_ipro_stg,
                    hosts,
                    os.path.join(es_dir, "documents"),
                ),
                kwargs=dict(
                    suffix=config["meta"]["release"],
                    create_indices=True,
                    outdir=outdir,
                    max_retries=0,
                    processes=6,
                    raise_on_error=False,
                    write_back=False,
                    alias="staging"
                ),
                scheduler=dict(queue=queue, cpu=6, mem=16000),
                requires=["init-elastic"]
            ),
            Task(
                name="complete-index-{}".format(i+1),
                fn=elastic.index_documents,
                args=(
                    my_ipro_stg,
                    hosts,
                    outdir,
                ),
                kwargs=dict(
                    suffix=config["meta"]["release"],
                    create_indices=False,
                    outdir=None,
                    max_retries=5,
                    processes=6,
                    raise_on_error=True,
                    write_back=True,
                    alias="staging"
                ),
                scheduler=dict(queue=queue, cpu=6, mem=16000),
                requires=[f"index-{i+1}", "create-documents"]
            ),
            Task(
                name="update-alias-{}".format(i+1),
                fn=elastic.update_alias,
                args=(my_ipro_stg, hosts),
                kwargs=dict(
                    suffix=config["meta"]["release"],
                    keep_prev_indices=False
                ),
                scheduler=dict(queue=queue),
                requires=[f"complete-index-{i+1}"]
            )
        ]

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
    wdb = os.path.join(wdir, "interpro7dw.db")
    with Workflow(tasks, db=wdb, name=wname, dir=wdir, daemon=args.daemon,
                  mail=notif) as w:
        success = w.run(args.tasks, secs=args.detach, resume=args.resume,
                        dry=args.dry_run, resubmit=args.retry)

    sys.exit(0 if success else 1)


def test_database_links():
    parser = argparse.ArgumentParser(
        description="Test Oracle public database links")
    parser.add_argument("config", metavar="config.ini",
                        help="configuration file")
    args = parser.parse_args()

    if not os.path.isfile(args.config):
        parser.error(f"cannot open '{args.config}': no such file or directory")

    config = configparser.ConfigParser()
    config.read(args.config)
    url = config["databases"]["interpro_production"]
    sys.exit(0 if oracle.utils.test_database_links(url) else 1)


def drop_database():
    parser = argparse.ArgumentParser(
        description="Drop offsite/fallback MySQL database")
    parser.add_argument("config", metavar="config.ini",
                        help="configuration file")
    parser.add_argument("-d", "--database", choices=("offsite", "fallback"))
    args = parser.parse_args()

    if not os.path.isfile(args.config):
        parser.error(f"cannot open '{args.config}': no such file or directory")

    config = configparser.ConfigParser()
    config.read(args.config)

    s = input(f"Do you want to drop the {args.database} database [y/N]? ")
    if s not in ('y', 'Y'):
        print("Aborted")
        return

    print("Dropping database")
    mysql.drop_database(config["database"]["interpro_" + args.database])
