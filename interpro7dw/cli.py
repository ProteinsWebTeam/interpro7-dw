#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import configparser
import os
from typing import List, Mapping

from mundone import Task, Workflow

from interpro7dw import __version__
from interpro7dw.ebi.interpro import elastic
from interpro7dw.ebi.interpro import email
from interpro7dw.ebi.interpro import ftp
from interpro7dw.ebi.interpro import production as ippro
from interpro7dw.ebi.interpro import staging
from interpro7dw.ebi import ebisearch, goa, pdbe, uniprot


class DataFiles:
    def __init__(self, path: str):
        os.makedirs(path, exist_ok=True)

        self.alignments = os.path.join(path, "alignments")
        self.clans = os.path.join(path, "sets")
        self.entries = os.path.join(path, "entries")
        self.entry2xrefs = os.path.join(path, "entry2xrefs")
        self.keys = os.path.join(path, "keys")
        self.interpro2taxonomy = os.path.join(path, "interpro2taxonomy")
        self.proteins = os.path.join(path, "proteins")
        self.proteomes = os.path.join(path, "proteomes")
        self.structures = os.path.join(path, "structures")
        self.taxonomy = os.path.join(path, "taxonomy")
        self.uniprot2comments = os.path.join(path, "uniprot2comments")
        self.uniprot2evidence = os.path.join(path, "uniprot2evidence")
        self.uniprot2features = os.path.join(path, "uniprot2features")
        self.uniprot2ida = os.path.join(path, "uniprot2ida")
        self.uniprot2matches = os.path.join(path, "uniprot2matches")
        self.uniprot2name = os.path.join(path, "uniprot2name")
        self.uniprot2proteome = os.path.join(path, "uniprot2proteome")
        self.uniprot2sequence = os.path.join(path, "uniprot2sequence")

        self.elastic = os.path.join(path, "elastic")
        self.ebisearch = os.path.join(path, "ebisearch")
        self.goa = os.path.join(path, "goa")
        self.pdbe = os.path.join(path, "pdbe")

        self.announcements = os.path.join(path, "announcements.txt")


def gen_tasks(config: configparser.ConfigParser) -> List[Task]:
    version = config["release"]["version"]
    release_date = config["release"]["date"]
    update_release = config.getboolean("release", "update")
    data_dir = config["data"]["path"]
    tmp_dir = config["data"]["tmp"]
    ipr_pro_url = config["databases"]["production"]
    ipr_stg_url = config["databases"]["staging"]
    ipr_rel_url = config["databases"]["release"]
    pfam_url = config["databases"]["pfam"]
    lsf_queue = config["workflow"]["lsf_queue"]
    pub_dir = os.path.join(config["exchange"]["interpro"], version)
    os.makedirs(pub_dir, mode=0o775, exist_ok=True)
    df = DataFiles(data_dir)

    es_dirs = [os.path.join(df.elastic, "default")]
    for cluster in config["elasticsearch"]:
        es_dirs.append(os.path.join(df.elastic, cluster))

    tasks = [
        # Export PDBe data
        Task(
            fn=pdbe.export_structures,
            args=(ipr_pro_url, df.structures),
            name="export-structures",
            scheduler=dict(mem=8000, queue=lsf_queue)
        ),

        # Export data from InterPro Oracle database
        Task(
            fn=ippro.export_clans,
            args=(ipr_pro_url, pfam_url, df.clans, df.alignments),
            name="export-clans",
            scheduler=dict(mem=1000, queue=lsf_queue)
        ),
        Task(
            fn=ippro.chunk_proteins,
            args=(ipr_pro_url, df.keys),
            name="init-export",
            scheduler=dict(mem=24000, queue=lsf_queue)
        ),
        Task(
            fn=ippro.export_proteins,
            args=(ipr_pro_url, df.keys, df.proteins),
            kwargs=dict(processes=8, tmpdir=tmp_dir),
            name="export-proteins",
            requires=["init-export"],
            scheduler=dict(cpu=8, mem=4000, scratch=4000, queue=lsf_queue)
        ),
        Task(
            fn=ippro.export_features,
            args=(ipr_pro_url, df.keys, df.uniprot2features),
            kwargs=dict(processes=8, tmpdir=tmp_dir),
            name="uniprot2features",
            requires=["init-export"],
            scheduler=dict(cpu=8, mem=4000, scratch=8000, queue=lsf_queue)
        ),
        Task(
            fn=ippro.export_matches,
            args=(ipr_pro_url, df.keys, df.uniprot2matches),
            kwargs=dict(processes=8, tmpdir=tmp_dir),
            name="uniprot2matches",
            requires=["init-export"],
            scheduler=dict(cpu=8, mem=10000, scratch=35000, queue=lsf_queue)
        ),
        Task(
            fn=ippro.export_sequences,
            args=(ipr_pro_url, df.keys, df.uniprot2sequence),
            kwargs=dict(processes=8, tmpdir=tmp_dir),
            name="uniprot2sequence",
            requires=["init-export"],
            scheduler=dict(cpu=8, mem=4000, scratch=40000, queue=lsf_queue)
        ),

        # Export data from UniProt Oracle database
        Task(
            fn=uniprot.export_proteomes,
            args=(ipr_pro_url, df.proteomes),
            name="export-proteomes",
            scheduler=dict(mem=100, queue=lsf_queue)
        ),
        Task(
            fn=ippro.export_taxonomy,
            args=(ipr_pro_url, df.taxonomy),
            name="export-taxonomy",
            scheduler=dict(mem=8000, queue=lsf_queue)
        ),
        Task(
            fn=uniprot.export_comments,
            args=(ipr_pro_url, df.keys, df.uniprot2comments),
            kwargs=dict(processes=8, tmpdir=tmp_dir),
            name="uniprot2comments",
            requires=["init-export"],
            scheduler=dict(cpu=8, mem=4000, scratch=2000, queue=lsf_queue)
        ),
        Task(
            fn=uniprot.export_name,
            args=(ipr_pro_url, df.keys, df.uniprot2name),
            kwargs=dict(processes=8, tmpdir=tmp_dir),
            name="uniprot2name",
            requires=["init-export"],
            scheduler=dict(cpu=8, mem=4000, scratch=2000, queue=lsf_queue)
        ),
        Task(
            fn=uniprot.export_evidence,
            args=(ipr_pro_url, df.keys, df.uniprot2evidence),
            kwargs=dict(processes=8, tmpdir=tmp_dir),
            name="uniprot2evidence",
            requires=["init-export"],
            scheduler=dict(cpu=8, mem=4000, scratch=2000, queue=lsf_queue)
        ),
        Task(
            fn=uniprot.export_proteome,
            args=(ipr_pro_url, df.keys, df.uniprot2proteome),
            kwargs=dict(processes=8, tmpdir=tmp_dir),
            name="uniprot2proteome",
            requires=["init-export"],
            scheduler=dict(cpu=8, mem=4000, scratch=1000, queue=lsf_queue)
        ),

        # Export signatures/entries with cross-references
        Task(
            fn=ippro.export_entries,
            args=(ipr_pro_url, config["metacyc"]["path"], df.clans,
                  df.proteins, df.structures, df.uniprot2matches,
                  df.uniprot2proteome, df.uniprot2ida, df.entry2xrefs,
                  df.entries),
            kwargs=dict(processes=8, tmpdir=tmp_dir),
            name="export-entries",
            scheduler=dict(cpu=8, mem=24000, scratch=60000, queue=lsf_queue),
            requires=["export-clans", "export-proteins", "export-structures",
                      "uniprot2matches", "uniprot2proteome"]
        ),

        # MySQL tables
        Task(
            fn=staging.insert_annotations,
            args=(ipr_pro_url, df.uniprot2matches, pfam_url, ipr_stg_url),
            name="insert-annotations",
            kwargs=dict(tmpdir=tmp_dir),
            scheduler=dict(cpu=2, mem=4000, scratch=40000, queue=lsf_queue),
            requires=["uniprot2matches"]
        ),
        Task(
            fn=staging.insert_clans,
            args=(ipr_stg_url, df.alignments, df.clans, df.entries,
                  df.entry2xrefs),
            kwargs=dict(tmpdir=tmp_dir),
            name="insert-clans",
            scheduler=dict(mem=16000, scratch=15000, queue=lsf_queue),
            requires=["export-entries"]
        ),
        Task(
            fn=staging.insert_databases,
            args=(ipr_pro_url, ipr_stg_url, version, release_date),
            kwargs=dict(update_prod=update_release),
            name="insert-databases",
            scheduler=dict(mem=100, queue=lsf_queue)
        ),
        Task(
            fn=staging.insert_structural_models,
            args=(ipr_pro_url, ipr_stg_url, df.entry2xrefs),
            name="insert-struct-models",
            scheduler=dict(mem=4000, queue=lsf_queue),
            requires=["export-entries"]
        ),
        Task(
            fn=staging.insert_entries,
            args=(pfam_url, ipr_stg_url, df.entries, df.entry2xrefs),
            name="insert-entries",
            scheduler=dict(mem=10000, queue=lsf_queue),
            requires=["insert-struct-models"]
        ),
        Task(
            fn=staging.insert_isoforms,
            args=(df.entries, ipr_pro_url, ipr_stg_url),
            name="insert-isoforms",
            scheduler=dict(mem=4000, queue=lsf_queue),
            requires=["export-entries"]
        ),
        Task(
            fn=staging.insert_proteins,
            args=(df.entries, df.proteins, df.structures, df.taxonomy,
                  df.uniprot2comments, df.uniprot2name, df.uniprot2evidence,
                  df.uniprot2ida, df.uniprot2matches, df.uniprot2proteome,
                  df.uniprot2sequence, ipr_pro_url, ipr_stg_url),
            name="insert-proteins",
            scheduler=dict(mem=8000, queue=lsf_queue),
            requires=["export-entries",  "export-taxonomy",
                      "uniprot2comments", "uniprot2name", "uniprot2evidence",
                      "uniprot2sequence", "insert-isoforms"]
        ),
        Task(
            fn=staging.insert_extra_features,
            args=(ipr_stg_url, df.uniprot2features),
            name="insert-features",
            scheduler=dict(mem=1000, queue=lsf_queue),
            requires=["uniprot2features"]
        ),
        Task(
            fn=staging.insert_residues,
            args=(ipr_pro_url, ipr_stg_url),
            kwargs=dict(tmpdir=tmp_dir),
            name="insert-residues",
            scheduler=dict(mem=2000, scratch=10000, queue=lsf_queue),
        ),
        Task(
            fn=staging.insert_proteomes,
            args=(df.entries, df.proteins, df.proteomes, df.structures,
                  df.uniprot2ida, df.uniprot2matches, df.uniprot2proteome,
                  ipr_stg_url),
            name="insert-proteomes",
            scheduler=dict(mem=28000, queue=lsf_queue),
            requires=["export-entries", "export-proteomes",
                      "export-structures"]
        ),
        Task(
            fn=staging.insert_structures,
            args=(df.entries, df.proteins, df.structures, df.uniprot2ida,
                  df.uniprot2matches, df.uniprot2proteome, ipr_stg_url),
            name="insert-structures",
            scheduler=dict(mem=8000, queue=lsf_queue),
            requires=["export-entries", "export-structures"]
        ),
        Task(
            fn=staging.insert_taxonomy,
            args=(df.entries, df.proteins, df.structures, df.taxonomy,
                  df.uniprot2matches, df.uniprot2proteome, ipr_stg_url,
                  df.interpro2taxonomy),
            kwargs=dict(tmpdir=tmp_dir),
            name="insert-taxonomy",
            scheduler=dict(mem=16000, scratch=20000, queue=lsf_queue),
            requires=["export-entries", "export-structures", "export-taxonomy"]
        ),

        Task(
            fn=staging.insert_release_notes,
            args=(df.entries, df.proteins, df.proteomes, df.structures,
                  df.taxonomy, df.uniprot2matches, df.uniprot2proteome,
                  ipr_rel_url, ipr_stg_url, df.announcements),
            name="insert-release-notes",
            scheduler=dict(mem=12000, queue=lsf_queue),
            requires=["export-entries", "export-proteomes",
                      "export-structures", "export-taxonomy",
                      "insert-databases"]
        ),

        # EBI Search
        Task(
            fn=ebisearch.export,
            args=(ipr_stg_url, df.entries, df.entry2xrefs, df.taxonomy,
                  df.ebisearch),
            name="export-ebisearch",
            scheduler=dict(mem=12000, queue=lsf_queue),
            requires=["insert-databases", "export-entries", "export-taxonomy"]
        ),
        Task(
            fn=ebisearch.publish,
            args=(df.ebisearch, config["exchange"]["ebisearch"]),
            name="publish-ebisearch",
            scheduler=dict(queue=lsf_queue),
            requires=["export-ebisearch"]
        ),

        # Export data for GOA
        Task(
            fn=goa.export,
            args=(ipr_pro_url, ipr_stg_url, df.goa),
            name="export-goa",
            scheduler=dict(mem=2000, queue=lsf_queue),
            requires=["insert-databases"]
        ),
        Task(
            fn=goa.publish,
            args=(df.goa, config["exchange"]["goa"]),
            name="publish-goa",
            scheduler=dict(queue=lsf_queue),
            requires=["export-goa"]
        ),

        # Export data from PDBe
        Task(
            fn=pdbe.export_pdb_matches,
            args=(ipr_pro_url, ipr_stg_url, df.pdbe),
            name="export-pdbe",
            # TODO: update resource requirements
            scheduler=dict(mem=2000, queue=lsf_queue),
            requires=["insert-databases"]
        ),
        Task(
            fn=pdbe.publish,
            args=(df.pdbe, config["exchange"]["pdbe"]),
            name="publish-pdbe",
            scheduler=dict(queue=lsf_queue),
            requires=["export-pdbe"]
        ),

        # Export data for Elastic
        Task(
            fn=elastic.export_documents,
            args=(df.proteins, df.entries, df.proteomes, df.structures,
                  df.taxonomy, df.uniprot2ida, df.uniprot2matches,
                  df.uniprot2proteome, es_dirs, version),
            name="es-export",
            scheduler=dict(mem=16000, queue=lsf_queue),
            requires=["export-entries", "export-proteomes", "export-taxonomy"]
        ),

        # Export files for FTP
        Task(
            fn=ftp.flatfiles.export,
            args=(df.entries, df.uniprot2matches, pub_dir),
            name="export-flat-files",
            scheduler=dict(mem=4000, queue=lsf_queue),
            requires=["export-entries"]
        ),
        Task(
            fn=ftp.relnotes.export,
            args=(ipr_stg_url, pub_dir),
            name="export-release-notes",
            scheduler=dict(mem=1000, queue=lsf_queue),
            requires=["insert-release-notes", "insert-isoforms"]
        ),
        Task(
            fn=ftp.uniparc.export_matches,
            args=(ipr_pro_url, pub_dir),
            kwargs=dict(processes=8, tmpdir=tmp_dir),
            name="export-uniparc-xml",
            scheduler=dict(cpu=8, mem=8000, scratch=80000, queue=lsf_queue)
        ),
        Task(
            fn=ftp.xmlfiles.export_features_matches,
            args=(ipr_pro_url, df.proteins, df.uniprot2features, pub_dir),
            kwargs=dict(processes=8),
            name="export-features-xml",
            scheduler=dict(cpu=8, mem=8000, queue=lsf_queue),
            requires=["insert-databases", "export-proteins",
                      "uniprot2features"]
        ),
        Task(
            fn=ftp.xmlfiles.export_interpro,
            args=(ipr_stg_url, df.entries, df.entry2xrefs,
                  df.interpro2taxonomy, pub_dir),
            kwargs=dict(tmpdir=tmp_dir),
            name="export-interpro-xml",
            scheduler=dict(mem=8000, scratch=20000, queue=lsf_queue),
            requires=["insert-databases", "insert-entries", "insert-taxonomy"]
        ),
        Task(
            fn=ftp.xmlfiles.export_matches,
            args=(ipr_pro_url, ipr_stg_url, df.proteins, df.uniprot2matches,
                  pub_dir),
            kwargs=dict(processes=8),
            name="export-matches-xml",
            scheduler=dict(cpu=8, mem=24000, queue=lsf_queue),
            requires=["insert-databases", "export-proteins", "uniprot2matches"]
        ),
        Task(
            fn=ftp.xmlfiles.export_structure_matches,
            args=(ipr_pro_url, df.proteins, df.structures, pub_dir),
            name="export-structures-xml",
            scheduler=dict(mem=8000, queue=lsf_queue),
            requires=["export-proteins", "export-structures"]
        )
    ]

    # Indexing data in Elastic
    for cluster, nodes in config.items("elasticsearch"):
        hosts = [host.strip() for host in nodes.split(',') if host.strip()]

        if not hosts:
            continue

        hosts = list(set(hosts))

        tasks += [
            Task(
                fn=elastic.create_indices,
                args=(ipr_stg_url, hosts, version),
                name=f"es-init-{cluster}",
                scheduler=dict(mem=100, queue=lsf_queue),
                requires=["insert-databases", "export-entries",
                          "export-proteomes", "export-taxonomy"]
            ),
            Task(
                fn=elastic.index_documents,
                args=(hosts, os.path.join(df.elastic, cluster), version),
                kwargs=dict(threads=8),
                name=f"es-index-{cluster}",
                scheduler=dict(mem=16000, queue=lsf_queue),
                requires=[f"es-init-{cluster}"]
            ),
            Task(
                fn=elastic.publish,
                args=(hosts,),
                name=f"es-publish-{cluster}",
                scheduler=dict(mem=100, queue=lsf_queue),
                requires=["es-export", f"es-index-{cluster}"]
            )
        ]

    # Notify production unfreeze
    email_serv = config["email"]["server"]
    email_port = int(config["email"]["port"])
    email_addr = config["email"]["address"]
    tasks.append(
        Task(
            fn=email.notify_curators,
            args=(email_serv, email_port, email_addr),
            name="notify-curators",
            scheduler=dict(queue=lsf_queue),
            requires=["export-features-xml", "export-goa",
                      "export-matches-xml", "export-structures-xml",
                      "export-uniparc-xml", "insert-annotations",
                      "insert-residues"]
        )
    )

    return tasks


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
    workflow_dir = config["workflow"]["path"]

    tasks = gen_tasks(config)

    database = os.path.join(workflow_dir, f"{version}.sqlite")
    with Workflow(tasks, dir=workflow_dir, database=database) as workflow:
        workflow.run(args.tasks, dry_run=args.dry_run, monitor=not args.detach)


def test_database_links():
    parser = argparse.ArgumentParser(
        description="Test Oracle public database links"
    )
    parser.add_argument("config",
                        metavar="config.ini",
                        help="configuration file")
    args = parser.parse_args()

    if not os.path.isfile(args.config):
        parser.error(f"cannot open '{args.config}': no such file or directory")

    config = configparser.ConfigParser()
    config.read(args.config)
    ippro.test_db_links(config["databases"]["production"])


def drop_database():
    parser = argparse.ArgumentParser(
        description="Drop release/fallback MySQL database"
    )
    parser.add_argument("config",
                        metavar="config.ini",
                        help="configuration file")
    parser.add_argument("database", choices=("release", "fallback"))
    args = parser.parse_args()

    if not os.path.isfile(args.config):
        parser.error(f"cannot open '{args.config}': no such file or directory")

    config = configparser.ConfigParser()
    config.read(args.config)

    s = input(f"Do you want to drop the {args.database} database [y/N]? ")
    if s not in ('y', 'Y'):
        print("Aborted")
        return

    print("dropping database")
    staging.drop_database(config["databases"][args.database])
    print("done")


def traverse_bottom_up(tasks: Mapping[str, Task], name: str) -> set:
    result = {name}
    for r in tasks[name].requires:
        result |= traverse_bottom_up(tasks, r)

    return result


def find_leaves(filepath, arg=None, exclude=None) -> List[Task]:
    """
    Example:
        Find tasks using production DB, except insert-proteins, so we know
        that once all these tasks are complete, curators can start
        integrating again

    find_leaves("/path/to/config.ini", "user/password@production",
                exclude=["insert-proteins"])
    """
    config = configparser.ConfigParser()
    config.read(filepath)

    tasks = {}
    selection = set()
    for t in gen_tasks(config):
        tasks[t.name] = t
        if exclude is not None and t.name in exclude:
            continue
        elif arg is None or arg in t.args or arg in t.kwargs:
            selection.add(t.name)

    # Remove non-leaves (no other task depend on them)
    remove = set()
    for name in selection:
        t = tasks[name]
        for r in t.requires:
            remove |= traverse_bottom_up(tasks, r)

    return sorted(selection - remove)
