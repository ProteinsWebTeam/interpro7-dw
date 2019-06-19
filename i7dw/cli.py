#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import configparser
import os
import re

from mundone import Task, Workflow

from . import ebisearch, uniprot, __version__
from .interpro import elastic, export, mysql, supermatch, xref


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

    tasks = [
        # Export data to stores
        Task(
            name="chunk-proteins",
            fn=export.chunk_proteins,
            args=(ora_ipro, os.path.join(export_dir, "chunks.json")),
            kwargs=dict(order_by=False),
            scheduler=dict(queue=queue, mem=12000),
        ),
        Task(
            name="export-matches",
            fn=export.export_protein2matches,
            args=(
                ora_ipro,
                os.path.join(export_dir, "chunks.json"),
                os.path.join(export_dir, "matches.dat")
            ),
            kwargs=dict(processes=4, tmpdir="/scratch"),
            scheduler=dict(queue=queue, mem=8000, scratch=20000, cpu=4),
            requires=["chunk-proteins"]
        ),
        Task(
            name="export-features",
            fn=export.export_protein2features,
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
            fn=export.export_protein2residues,
            args=(
                ora_ipro,
                os.path.join(export_dir, "chunks.json"),
                os.path.join(export_dir, "residues.dat")
            ),
            kwargs=dict(processes=4, tmpdir="/scratch"),
            scheduler=dict(queue=queue, mem=3000, scratch=8000, cpu=4),
            requires=["chunk-proteins"]
        ),
        Task(
            name="export-proteins",
            fn=export.export_proteins,
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
            fn=export.export_sequences,
            args=(
                ora_ipro,
                os.path.join(export_dir, "chunks.json"),
                os.path.join(export_dir, "sequences.dat")
            ),
            kwargs=dict(processes=4, tmpdir="/scratch"),
            scheduler=dict(queue=queue, mem=2000, scratch=30000, cpu=4),
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
            kwargs=dict(processes=4, tmpdir="/scratch"),
            scheduler=dict(queue=queue, mem=2000, scratch=3000, cpu=4),
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
            kwargs=dict(processes=4, tmpdir="/scratch"),
            scheduler=dict(queue=queue, mem=1000, scratch=1000, cpu=4),
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
            kwargs=dict(processes=4, tmpdir="/scratch"),
            scheduler=dict(queue=queue, mem=500, scratch=500, cpu=4),
            requires=["chunk-proteins"]
        ),

        # Load data to MySQL
        Task(
            name="init-tables",
            fn=mysql.init,
            args=(my_ipro_stg,),
            scheduler=dict(queue=queue)
        ),
        Task(
            name="insert-taxa",
            fn=mysql.taxonomy.insert_taxa,
            args=(ora_ipro, my_ipro_stg),
            scheduler=dict(queue=queue, mem=4000),
            requires=["init-tables"]
        ),
        Task(
            name="insert-proteomes",
            fn=mysql.proteome.insert_proteomes,
            args=(ora_ipro, my_ipro_stg),
            scheduler=dict(queue=queue, mem=2000),
            requires=["insert-taxa"]
        ),
        Task(
            name="insert-databases",
            fn=mysql.database.insert_databases,
            args=(ora_ipro, my_ipro_stg, config["meta"]["release"],
                  config["meta"]["release_date"],),
            scheduler=dict(queue=queue),
            requires=["init-tables"]
        ),
        Task(
            name="insert-entries",
            fn=mysql.entry.insert_entries,
            args=(ora_ipro, my_pfam, my_ipro_stg),
            scheduler=dict(queue=queue, mem=2000),
            requires=["insert-databases"]
        ),
        Task(
            name="insert-annotations",
            fn=mysql.entry.insert_annotations,
            args=(my_pfam, my_ipro_stg),
            scheduler=dict(queue=queue, mem=4000),
            requires=["insert-entries"]
        ),
        Task(
            name="insert-structures",
            fn=mysql.structure.insert_structures,
            args=(ora_ipro, my_ipro_stg),
            scheduler=dict(queue=queue, mem=4000),
            requires=["insert-databases"]
        ),
        Task(
            name="insert-sets",
            fn=mysql.entry.insert_sets,
            args=(ora_ipro, my_pfam, my_ipro_stg),
            kwargs=dict(tmpdir="/scratch"),
            scheduler=dict(queue=queue, mem=3000, scratch=3000),
            requires=["insert-entries"]
        ),

        Task(
            name="export-ida",
            fn=export.export_ida,
            args=(
                my_ipro_stg,
                os.path.join(export_dir, "matches.dat"),
                os.path.join(export_dir, "ida.dat")
            ),
            kwargs=dict(processes=4, tmpdir="/scratch"),
            scheduler=dict(queue=queue, mem=8000, scratch=4000, cpu=4),
            requires=["export-matches", "insert-entries"]
        ),
        Task(
            name="insert-isoforms",
            fn=mysql.protein.insert_isoforms,
            args=(ora_ipro, my_ipro_stg),
            scheduler=dict(queue=queue, mem=4000),
            requires=["insert-entries"]
        ),
        Task(
            name="insert-proteins",
            fn=mysql.protein.insert_proteins,
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
                os.path.join(export_dir, "matches.dat"),
                os.path.join(export_dir, "ida.dat")
            ),
            scheduler=dict(queue=queue, mem=24000),
            requires=[
                "insert-structures", "insert-taxa", "insert-sets",
                "insert-isoforms", "export-proteins", "export-sequences",
                "export-misc", "export-names", "export-comments",
                "export-proteomes", "export-residues", "export-features",
                "export-ida"
            ]
        ),
        Task(
            name="release-notes",
            fn=mysql.relnote.make_release_notes,
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
                "export-proteins", "export-matches", "export-proteomes"
            ]
        ),

        # Overlapping homologous superfamilies
        Task(
            name="overlapping-families",
            fn=supermatch.calculate_relationships,
            args=(
                my_ipro_stg,
                os.path.join(export_dir, "proteins.dat"),
                os.path.join(export_dir, "matches.dat"),
                config.getfloat("jaccard", "threshold")
            ),
            # kwargs=dict(ora_uri=ora_ipro),
            scheduler=dict(queue=queue, mem=6000),
            requires=["export-proteins", "export-matches", "insert-entries"]
        ),

        # Mappings for GOA team
        Task(
            name="export-goa",
            fn=xref.export_goa_mappings,
            args=(my_ipro_stg, ora_ipro, config["export"]["goa"]),
            scheduler=dict(queue=queue, mem=2000),
            requires=["insert-databases"]
        ),

        # Cross-references
        Task(
            name="export-xrefs",
            fn=xref.export,
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
            kwargs=dict(tmpdir="/scratch"),
            scheduler=dict(queue=queue, mem=24000, scratch=20000, cpu=5),
            requires=[
                "export-proteins", "export-matches", "export-proteomes",
                "insert-entries", "insert-proteomes", "insert-sets",
                "insert-structures"
            ]
        ),

        Task(
            name="update-counts",
            fn=mysql.update_counts,
            args=(
                my_ipro_stg,
                os.path.join(export_dir, "entries_xref.dat"),
                os.path.join(export_dir, "proteomes_xref.dat"),
                os.path.join(export_dir, "structures_xref.dat"),
                os.path.join(export_dir, "taxa_xref.dat")
            ),
            kwargs=dict(tmpdir="/scratch"),
            scheduler=dict(queue=queue, mem=16000, scratch=15000),
            requires=["export-xrefs", "overlapping-families"]
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
                config["ebisearch"]["stg"]
            ),
            kwargs=dict(processes=4, by_type=True),
            scheduler=dict(queue=queue, mem=16000, cpu=4),
            requires=["export-xrefs"],
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
            fn=elastic.init_dir,
            args=(os.path.join(es_dir, "documents"),),
            scheduler=dict(queue=queue),
            requires=[
                "insert-entries", "insert-sets", "insert-proteomes",
                "export-proteins", "export-names", "export-comments",
                "export-proteomes", "export-matches"
            ]
        ),
        Task(
            name="create-documents",
            fn=elastic.relationship.create_documents,
            args=(
                ora_ipro,
                my_ipro_stg,
                os.path.join(export_dir, "proteins.dat"),
                os.path.join(export_dir, "names.dat"),
                os.path.join(export_dir, "comments.dat"),
                os.path.join(export_dir, "proteomes.dat"),
                os.path.join(export_dir, "matches.dat"),
                os.path.join(es_dir, "documents")
            ),
            kwargs=dict(processes=8),
            scheduler=dict(queue=queue, cpu=8, mem=32000),
            requires=["init-elastic"]
        )
    ]

    for i, hosts in enumerate(es_clusters):
        hosts = list(set(hosts.split(',')))
        dst = os.path.join(es_dir, "cluster-" + str(i+1))

        tasks += [
            Task(
                name="index-" + str(i+1),
                fn=elastic.relationship.index_documents,
                args=(
                    my_ipro_stg,
                    hosts,
                    os.path.join(es_dir, "documents")
                ),
                kwargs=dict(
                    body_path=config["elastic"]["body"],
                    num_shards=config.getint("elastic", "shards"),
                    shards_path=config["elastic"]["indices"],
                    suffix=config["meta"]["release"],
                    dst=dst,
                    processes=6,
                    raise_on_error=False
                ),
                scheduler=dict(queue=queue, cpu=6, mem=16000),
                requires=["init-elastic"]
            ),
            Task(
                name="complete-index-{}".format(i+1),
                fn=elastic.relationship.index_documents,
                args=(
                    my_ipro_stg,
                    hosts,
                    dst
                ),
                kwargs=dict(
                    suffix=config["meta"]["release"],
                    processes=6,
                    write_back=True,
                    max_retries=5,
                    alias="staging"
                ),
                scheduler=dict(queue=queue, cpu=6, mem=16000),
                requires=["index-{}".format(i+1), "create-documents"]
            ),
            Task(
                name="update-alias-{}".format(i+1),
                fn=elastic.relationship.update_alias,
                args=(my_ipro_stg, hosts),
                kwargs=dict(
                    suffix=config["meta"]["release"],
                    delete_removed=True
                ),
                scheduler=dict(queue=queue),
                requires=["complete-index-{}".format(i+1)]
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
    with Workflow(tasks, name=wname, dir=wdir, daemon=args.daemon,
                  mail=notif) as w:
        w.run(
            args.tasks,
            secs=args.detach,
            resume=args.resume,
            dry=args.dry_run,
            resubmit=args.retry
        )


def test_db_links():
    # Lazy loading
    from cx_Oracle import DatabaseError
    from .dbms import connect

    parser = argparse.ArgumentParser(
        description="Build InterPro7 data warehouse"
    )
    parser.add_argument("config",
                        metavar="config.ini",
                        help="configuration file")

    args = parser.parse_args()

    if not os.path.isfile(args.config):
        raise RuntimeError(
            "cannot open '{}': "
            "no such file or directory".format(args.config)
        )

    config = configparser.ConfigParser()
    config.read(args.config)

    has_errors = False
    con, cur = connect(config["databases"]["interpro_oracle"])
    for link in ("GOAPRO", "PDBE_LIVE", "SWPREAD"):
        try:
            cur.execute("SELECT * FROM DUAL@{}".format(link))
        except DatabaseError:
            print("{:<15} error".format(link))
            has_errors = True
        else:
            print("{:<15} ok".format(link))

    cur.close()
    con.close()

    exit(1 if has_errors else 0)
