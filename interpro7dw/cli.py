#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import configparser
import os

from mundone import Task, Workflow

from interpro7dw import __version__
from interpro7dw.ebi.interpro import elastic
from interpro7dw.ebi.interpro import ftp
from interpro7dw.ebi.interpro import production as ippro
from interpro7dw.ebi.interpro import staging
from interpro7dw.ebi import ebisearch, goa, pdbe, uniprot


class DataFiles(object):
    def __init__(self, path: str):
        os.makedirs(path, exist_ok=True)

        self.clans = os.path.join(path, "sets")
        self.entries = os.path.join(path, "entries")
        self.entry2xrefs = os.path.join(path, "entry2xrefs")
        self.keys = os.path.join(path, "keys")
        self.overlapping = os.path.join(path, "overlapping")
        self.proteins = os.path.join(path, "proteins")
        self.proteomes = os.path.join(path, "proteomes")
        self.structures = os.path.join(path, "structures")
        self.taxonomy = os.path.join(path, "taxonomy")
        self.uniprot2comments = os.path.join(path, "uniprot2comments")
        self.uniprot2entries = os.path.join(path, "uniprot2entries")
        self.uniprot2evidence = os.path.join(path, "uniprot2evidence")
        self.uniprot2features = os.path.join(path, "uniprot2features")
        self.uniprot2ida = os.path.join(path, "uniprot2ida")
        self.uniprot2matches = os.path.join(path, "uniprot2matches")
        self.uniprot2name = os.path.join(path, "uniprot2name")
        self.uniprot2proteome = os.path.join(path, "uniprot2proteome")
        self.uniprot2residues = os.path.join(path, "uniprot2residues")
        self.uniprot2sequence = os.path.join(path, "uniprot2sequence")

        self.es_rel = os.path.join(path, "elastic", "rel")
        self.es_ida = os.path.join(path, "elastic", "ida")


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
    release_date = config["release"]["date"]
    data_dir = config["data"]["path"]
    tmp_dir = config["data"]["tmp"]
    ipr_pro_url = config["databases"]["production"]
    ipr_stg_url = config["databases"]["staging"]
    ipr_rel_url = config["databases"]["release"]
    pfam_url = config["databases"]["pfam"]
    lsf_queue = config["workflow"]["lsf_queue"]
    workflow_dir = config["workflow"]["path"]

    pub_dir = os.path.join(config["exchange"]["interpro"], version)
    os.makedirs(pub_dir, mode=0o775, exist_ok=True)

    df = DataFiles(data_dir)
    tasks = [
        # Populate 'independent' MySQL tables (do not rely on other tasks)
        Task(
            fn=staging.insert_databases,
            args=(ipr_pro_url, ipr_stg_url, version, release_date),
            name="insert-databases",
            scheduler=dict(mem=100, queue=lsf_queue)
        ),
        Task(
            fn=staging.insert_annotations,
            args=(pfam_url, ipr_stg_url),
            name="insert-annotations",
            kwargs=dict(dir=tmp_dir),
            scheduler=dict(cpu=2, mem=8000, scratch=40000, queue=lsf_queue)
        ),
        Task(
            fn=staging.init_clans,
            args=(ipr_pro_url, ipr_stg_url, df.clans),
            name="init-clans",
            scheduler=dict(mem=1000, queue=lsf_queue)
        ),

        # Chunk UniProt accessions
        Task(
            fn=ippro.chunk_proteins,
            args=(ipr_pro_url, df.keys),
            name="init-export",
            scheduler=dict(mem=16000, queue=lsf_queue)
        ),

        # Export data from InterPro Oracle database
        Task(
            fn=ippro.export_proteins,
            args=(ipr_pro_url, df.keys, df.proteins),
            kwargs=dict(dir=tmp_dir, processes=8),
            name="export-proteins",
            requires=["init-export"],
            scheduler=dict(cpu=8, mem=4000, scratch=2000, queue=lsf_queue)
        ),
        Task(
            fn=ippro.export_features,
            args=(ipr_pro_url, df.keys, df.uniprot2features),
            kwargs=dict(dir=tmp_dir, processes=8),
            name="uniprot2features",
            requires=["init-export"],
            scheduler=dict(cpu=8, mem=4000, scratch=8000, queue=lsf_queue)
        ),
        Task(
            fn=ippro.export_matches,
            args=(ipr_pro_url, df.keys, df.uniprot2matches),
            kwargs=dict(dir=tmp_dir, processes=8),
            name="uniprot2matches",
            requires=["init-export"],
            scheduler=dict(cpu=8, mem=8000, scratch=25000, queue=lsf_queue)
        ),
        Task(
            fn=ippro.export_residues,
            args=(ipr_pro_url, df.keys, df.uniprot2residues),
            kwargs=dict(dir=tmp_dir, processes=8),
            name="uniprot2residues",
            requires=["init-export"],
            scheduler=dict(cpu=8, mem=8000, scratch=16000, queue=lsf_queue)
        ),
        Task(
            fn=ippro.export_sequences,
            args=(ipr_pro_url, df.keys, df.uniprot2sequence),
            kwargs=dict(dir=tmp_dir, processes=8),
            name="uniprot2sequence",
            requires=["init-export"],
            scheduler=dict(cpu=8, mem=4000, scratch=40000, queue=lsf_queue)
        ),

        # Export data from UniProt Oracle database
        Task(
            fn=uniprot.export_comments,
            args=(ipr_pro_url, df.keys, df.uniprot2comments),
            kwargs=dict(dir=tmp_dir, processes=8),
            name="uniprot2comments",
            requires=["init-export"],
            scheduler=dict(cpu=8, mem=4000, scratch=2000, queue=lsf_queue)
        ),
        Task(
            fn=uniprot.export_name,
            args=(ipr_pro_url, df.keys, df.uniprot2name),
            kwargs=dict(dir=tmp_dir, processes=8),
            name="uniprot2name",
            requires=["init-export"],
            scheduler=dict(cpu=8, mem=4000, scratch=2000, queue=lsf_queue)
        ),
        Task(
            fn=uniprot.export_evidence,
            args=(ipr_pro_url, df.keys, df.uniprot2evidence),
            kwargs=dict(dir=tmp_dir, processes=8),
            name="uniprot2evidence",
            requires=["init-export"],
            scheduler=dict(cpu=8, mem=4000, scratch=1000, queue=lsf_queue)
        ),
        Task(
            fn=uniprot.export_proteome,
            args=(ipr_pro_url, df.keys, df.uniprot2proteome),
            kwargs=dict(dir=tmp_dir, processes=8),
            name="uniprot2proteome",
            requires=["init-export"],
            scheduler=dict(cpu=8, mem=4000, scratch=1000, queue=lsf_queue)
        ),
        Task(
            fn=ippro.export_taxonomy,
            args=(ipr_pro_url, df.taxonomy),
            name="export-taxonomy",
            scheduler=dict(mem=4000, queue=lsf_queue)
        ),
        Task(
            fn=ippro.export_entries,
            args=(ipr_pro_url, config["MetaCyc"]["path"],
                  df.clans, df.uniprot2matches, df.entries,
                  df.uniprot2entries),
            kwargs=dict(dir=tmp_dir, processes=8),
            name="export-entries",
            scheduler=dict(cpu=8, mem=8000, scratch=10000, queue=lsf_queue),
            requires=["init-clans", "uniprot2matches"]
        ),
        Task(
            fn=uniprot.export_proteomes,
            args=(ipr_pro_url, df.proteomes),
            name="export-proteomes",
            scheduler=dict(mem=100, queue=lsf_queue)
        ),

        # Export PDBe data
        Task(
            fn=pdbe.export_structures,
            args=(ipr_pro_url, df.structures),
            name="export-structures",
            scheduler=dict(mem=8000, queue=lsf_queue)
        ),

        # Export domain architectures
        Task(
            fn=staging.export_ida,
            args=(df.entries, df.uniprot2matches, df.uniprot2ida),
            kwargs=dict(dir=tmp_dir, processes=8),
            name="uniprot2ida",
            scheduler=dict(cpu=8, mem=8000, scratch=10000, queue=lsf_queue),
            requires=["export-entries", "uniprot2matches"]
        ),

        # MySQL tables
        Task(
            fn=staging.insert_entries,
            args=(ipr_stg_url, df.entries, df.proteins, df.structures,
                  df.uniprot2ida, df.uniprot2matches, df.uniprot2proteome,
                  df.entry2xrefs),
            kwargs=dict(dir=tmp_dir),
            # kwargs=dict(dir=tmp_dir, pro_url=ipr_pro_url),
            name="insert-entries",
            scheduler=dict(mem=8000, scratch=16000, queue=lsf_queue),
            requires=["export-entries", "export-proteins", "export-structures",
                      "uniprot2ida", "uniprot2matches", "uniprot2proteome"]
        ),
        Task(
            fn=staging.insert_clans,
            args=(ipr_stg_url, df.clans, df.entries, df.entry2xrefs),
            kwargs=dict(dir=tmp_dir),
            name="insert-clans",
            scheduler=dict(mem=16000, scratch=4000, queue=lsf_queue),
            requires=["init-clans", "export-entries", "insert-entries"]
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
            args=(df.proteins, df.structures, df.taxonomy, df.uniprot2comments,
                  df.uniprot2name, df.uniprot2entries,
                  df.uniprot2evidence, df.uniprot2features, df.uniprot2ida,
                  df.uniprot2proteome, df.uniprot2residues,
                  df.uniprot2sequence, ipr_pro_url, ipr_stg_url),
            name="insert-proteins",
            scheduler=dict(mem=8000, queue=lsf_queue),
            requires=["export-proteins", "export-structures",
                      "export-taxonomy", "uniprot2comments", "uniprot2name",
                      "export-entries", "uniprot2evidence", "uniprot2features",
                      "uniprot2ida", "uniprot2proteome", "uniprot2residues",
                      "uniprot2sequence", "insert-isoforms"]
        ),
        Task(
            fn=staging.insert_proteomes,
            args=(df.proteomes, df.structures, df.proteins, df.uniprot2ida,
                  df.uniprot2entries, df.uniprot2proteome, ipr_stg_url),
            name="insert-proteomes",
            scheduler=dict(mem=24000, queue=lsf_queue),
            requires=["export-proteomes", "export-structures",
                      "export-proteins", "uniprot2ida", "export-entries",
                      "uniprot2proteome"]
        ),
        Task(
            fn=staging.insert_structures,
            args=(df.entries, df.proteins, df.structures, df.uniprot2ida,
                  df.uniprot2matches, df.uniprot2proteome, ipr_stg_url),
            name="insert-structures",
            scheduler=dict(mem=8000, queue=lsf_queue),
            requires=["export-entries", "export-proteins", "export-structures",
                      "uniprot2ida", "uniprot2matches", "uniprot2proteome"]
        ),
        Task(
            fn=staging.insert_taxonomy,
            args=(df.entries, df.proteins, df.structures, df.taxonomy,
                  df.uniprot2matches, df.uniprot2proteome, ipr_stg_url),
            kwargs=dict(dir=tmp_dir),
            name="insert-taxonomy",
            scheduler=dict(mem=16000, scratch=20000, queue=lsf_queue),
            requires=["export-entries", "export-proteins", "export-structures",
                      "export-taxonomy", "uniprot2matches", "uniprot2proteome"]
        ),

        Task(
            fn=staging.make_release_notes,
            args=(df.entries, df.proteins, df.proteomes, df.structures,
                  df.taxonomy, df.uniprot2entries, df.uniprot2proteome,
                  ipr_rel_url, ipr_stg_url),
            name="insert-release-notes",
            scheduler=dict(mem=8000, queue=lsf_queue),
            requires=["export-entries", "export-proteins", "export-proteomes",
                      "export-structures", "export-taxonomy",
                      "uniprot2proteome"]
        ),

        # EBI Search
        Task(
            fn=ebisearch.export,
            args=(ipr_stg_url, df.entries, df.entry2xrefs,
                  config["ebisearch"]["staging"]),
            name="ebisearch",
            scheduler=dict(mem=12000, queue=lsf_queue),
            requires=["insert-databases", "insert-taxonomy", "export-entries",
                      "insert-entries"]
        ),
        Task(
            fn=ebisearch.publish,
            args=(config["ebisearch"]["staging"],
                  config["ebisearch"]["release"]),
            name="publish-ebisearch",
            scheduler=dict(queue=lsf_queue),
            requires=["ebisearch"]
        ),

        # Export data for Elastic
        Task(
            fn=elastic.relationship.dump_documents,
            args=(df.proteins, df.entries, df.proteomes, df.structures,
                  df.taxonomy, df.uniprot2ida, df.uniprot2matches,
                  df.uniprot2proteome, os.path.join(df.es_rel, "all"), version),
            name="es-rel",
            scheduler=dict(mem=16000, queue=lsf_queue),
            requires=["export-proteins", "export-entries", "export-proteomes",
                      "export-structures", "export-taxonomy", "uniprot2ida",
                      "uniprot2matches", "uniprot2proteome"]
        ),
        Task(
            fn=elastic.ida.dump_documents,
            args=(df.uniprot2ida, os.path.join(df.es_ida, "all"), version),
            name="es-ida",
            scheduler=dict(mem=2000, queue=lsf_queue),
            requires=["uniprot2ida"]
        ),

        # Export data for GOA
        Task(
            fn=goa.export,
            args=(ipr_pro_url, ipr_stg_url, config["exchange"]["goa"]),
            name="export-goa",
            scheduler=dict(mem=2000, queue=lsf_queue),
            requires=["insert-databases"]
        ),

        # Export files for FTP
        Task(
            fn=ftp.flatfiles.export,
            args=(df.entries, df.uniprot2matches, pub_dir),
            name="export-flat-files",
            scheduler=dict(mem=4000, queue=lsf_queue),
            requires=["export-entries", "uniprot2matches"]
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
            kwargs=dict(dir=tmp_dir, processes=8),
            name="export-uniparc-xml",
            scheduler=dict(cpu=8, mem=40000, scratch=130000, queue=lsf_queue)
        ),
        Task(
            fn=ftp.xmlfiles.export_features_matches,
            args=(ipr_pro_url, df.proteins, df.uniprot2features, pub_dir),
            kwargs=dict(processes=8),
            name="export-extra-xml",
            scheduler=dict(cpu=8, mem=24000, queue=lsf_queue),
            requires=["insert-databases", "export-proteins",
                      "uniprot2features"]
        ),
        Task(
            fn=ftp.xmlfiles.export_interpro,
            args=(ipr_stg_url, df.entries, df.entry2xrefs, pub_dir),
            kwargs=dict(dir=tmp_dir),
            name="export-interpro-xml",
            scheduler=dict(mem=8000, scratch=20000, queue=lsf_queue),
            requires=["insert-databases", "insert-entries", "insert-taxonomy"]
        ),
        Task(
            fn=ftp.xmlfiles.export_matches,
            args=(ipr_pro_url, ipr_stg_url, df.proteins, df.uniprot2matches,
                  pub_dir),
            kwargs=dict(processes=8),
            name="export-match-xml",
            scheduler=dict(cpu=8, mem=24000, queue=lsf_queue),
            requires=["insert-databases", "export-proteins", "uniprot2matches"]
        ),
        Task(
            fn=ftp.xmlfiles.export_structure_matches,
            args=(ipr_pro_url, df.proteins, df.structures, pub_dir),
            name="export-feature-xml",
            scheduler=dict(mem=8000, queue=lsf_queue),
            requires=["export-proteins", "export-structures"]
        ),

    ]

    # Indexing data in Elastic
    for cluster, nodes in config.items("elasticsearch"):
        hosts = [host.strip() for host in nodes.split(',') if host.strip()]

        if not hosts:
            continue

        hosts = list(set(hosts))

        tasks += [
            Task(
                fn=elastic.relationship.index_documents,
                args=(ipr_stg_url, hosts, os.path.join(df.es_rel, "all"),
                      version, os.path.join(df.es_rel, cluster)),
                name=f"es-rel-{cluster}",
                scheduler=dict(mem=4000, queue=lsf_queue),
                requires=["export-proteins", "export-entries",
                          "export-proteomes", "export-structures",
                          "export-taxonomy", "uniprot2ida", "uniprot2matches",
                          "uniprot2proteome"]
            ),
            Task(
                fn=elastic.ida.index_documents,
                args=(hosts, os.path.join(df.es_ida, "all"), version,
                      os.path.join(df.es_ida, cluster)),
                name=f"es-ida-{cluster}",
                scheduler=dict(mem=1000, queue=lsf_queue),
                requires=["uniprot2ida"]
            ),
            Task(
                fn=elastic.publish,
                args=(hosts,),
                name=f"publish-es-{cluster}",
                requires=["es-rel", f"es-rel-{cluster}",
                          "es-ida", f"es-ida-{cluster}"]
            ),
        ]

    database = os.path.join(workflow_dir, f"{version}.sqlite")
    with Workflow(tasks, dir=workflow_dir, database=database) as workflow:
        workflow.run(args.tasks, dry_run=args.dry_run, monitor=not args.detach)


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
