#!/usr/bin/env python

import argparse
import configparser
import os
import time
from typing import List, Mapping, Optional, Sequence

from mundone import Task, Workflow

from interpro7dw import __version__
from interpro7dw import interpro, pdbe, uniprot


class DataFiles:
    def __init__(self, root: str, pub_dir: str, create_dirs: bool = True):
        if create_dirs:
            os.makedirs(root, exist_ok=True)
            os.makedirs(pub_dir, mode=0o775, exist_ok=True)

        # Stores
        self.alignments = os.path.join(root, "alignments")
        self.proteins = os.path.join(root, "proteins")
        self.protein2domorg = os.path.join(root, "protein2domorg")
        self.protein2evidence = os.path.join(root, "protein2evidence")
        self.protein2features = os.path.join(root, "protein2features")
        self.protein2function = os.path.join(root, "protein2function")
        self.protein2matches = os.path.join(root, "protein2matches")
        self.protein2name = os.path.join(root, "protein2name")
        self.protein2proteome = os.path.join(root, "protein2proteome")
        self.protein2residues = os.path.join(root, "protein2residues")
        self.protein2sequence = os.path.join(root, "protein2sequence")

        # SimpleStores
        self.entryxrefs = os.path.join(root, "entryxrefs")
        self.isoforms = os.path.join(root, "isoforms")
        self.uniparc = os.path.join(root, "uniparc")

        # Data dumps
        self.clans = os.path.join(root, "clans")
        self.databases = os.path.join(root, "databases")
        self.overlapping_entries = os.path.join(root, "overlapping")
        self.proteomes = os.path.join(root, "proteomes")
        self.structures = os.path.join(root, "structures")
        self.taxa = os.path.join(root, "taxa")

        # Files for FTP
        self.pub_uniparc = os.path.join(pub_dir, "uniparc_match.tar.gz")


def gen_tasks(config: configparser.ConfigParser) -> List[Task]:
    release_version = config["release"]["version"]
    release_date = config["release"]["date"]
    data_dir = config["data"]["path"]
    temp_dir = config["data"]["tmp"]
    ipr_pro_url = config["databases"]["interpro_production"]
    pdbe_url = config["databases"]["pdbe"]
    pfam_url = config["databases"]["pfam"]
    uniprot_url = config["databases"]["uniprot"]
    pub_dir = config["exchange"]["interpro"]
    lsf_queue = config["workflow"]["lsf_queue"]

    df = DataFiles(data_dir, os.path.join(pub_dir, release_version))

    tasks = [
        # Data from InterPro
        Task(fn=interpro.oracle.entries.dump_databases,
             args=(ipr_pro_url, release_version, release_date, df.databases),
             kwargs=dict(update=config.getboolean("release", "update")),
             name="export-databases",
             scheduler=dict(mem=100, queue=lsf_queue)),
        Task(fn=interpro.oracle.taxa.export_taxa,
             args=(ipr_pro_url, df.taxa),
             name="export-taxa",
             scheduler=dict(mem=8000, queue=lsf_queue)),
        Task(fn=interpro.oracle.proteins.export_isoforms,
             args=(ipr_pro_url, df.isoforms),
             name="export-isoforms",
             # todo: review
             scheduler=dict(mem=8000, queue=lsf_queue)),
        Task(fn=interpro.oracle.proteins.export_proteins,
             args=(ipr_pro_url, df.proteins),
             kwargs=dict(tempdir=temp_dir, workers=8),
             name="export-proteins",
             scheduler=dict(cpu=8, mem=2000, scratch=10000, queue=lsf_queue)),
        Task(fn=interpro.oracle.proteins.export_features,
             args=(ipr_pro_url, df.proteins, df.protein2features),
             kwargs=dict(tempdir=temp_dir, workers=8),
             name="export-features",
             requires=["export-proteins"],
             scheduler=dict(cpu=8, mem=4000, scratch=10000, queue=lsf_queue)),
        Task(fn=interpro.oracle.proteins.export_matches,
             args=(ipr_pro_url, df.proteins, df.protein2matches),
             kwargs=dict(tempdir=temp_dir, workers=8),
             name="export-matches",
             requires=["export-proteins"],
             scheduler=dict(cpu=8, mem=4000, scratch=50000, queue=lsf_queue)),
        Task(fn=interpro.oracle.proteins.export_residues,
             args=(ipr_pro_url, df.proteins, df.protein2residues),
             kwargs=dict(tempdir=temp_dir, workers=8),
             name="export-residues",
             requires=["export-proteins"],
             scheduler=dict(cpu=8, mem=12000, scratch=20000, queue=lsf_queue)),
        Task(fn=interpro.oracle.proteins.export_sequences,
             args=(ipr_pro_url, df.proteins, df.protein2sequence),
             kwargs=dict(tempdir=temp_dir, workers=8),
             name="export-sequences",
             requires=["export-proteins"],
             scheduler=dict(cpu=8, mem=2000, scratch=50000, queue=lsf_queue)),
        Task(fn=interpro.oracle.entries.dump_domain_organisation,
             args=(ipr_pro_url, df.proteins, df.protein2matches,
                   df.protein2domorg),
             kwargs=dict(tempdir=temp_dir, workers=8),
             name="export-dom-orgs",
             requires=["export-matches"],
             scheduler=dict(cpu=8, mem=4000, scratch=20000, queue=lsf_queue)),
        Task(fn=interpro.oracle.entries.dump_similar_entries,
             args=(ipr_pro_url, df.protein2matches, df.overlapping_entries),
             name="export-sim-entries",
             requires=["export-matches"],
             scheduler=dict(mem=2000, queue=lsf_queue)),
        Task(fn=interpro.oracle.proteins.export_uniparc,
             args=(ipr_pro_url, df.uniparc),
             kwargs=dict(tempdir=temp_dir, workers=8),
             name="export-uniparc",
             # todo: review
             scheduler=dict(cpu=8, mem=24000, scratch=50000, queue=lsf_queue)),

        # Data from InterPro + other sources (Pfam, PDBe)
        Task(fn=interpro.oracle.clans.export_clans,
             args=(ipr_pro_url, pfam_url, df.clans, df.alignments),
             name="export-clans",
             scheduler=dict(mem=2000, queue=lsf_queue)),
        Task(fn=pdbe.export_structures,
             args=(ipr_pro_url, pdbe_url, df.structures),
             name="export-structures",
             scheduler=dict(mem=8000, queue=lsf_queue)),

        # Data from UniProt
        Task(fn=uniprot.proteomes.export_proteomes,
             args=(uniprot_url, df.proteomes),
             name="export-reference-proteomes",
             scheduler=dict(mem=100, queue=lsf_queue)),
        Task(fn=uniprot.proteins.export_entry2evidence,
             args=(uniprot_url, df.proteins, df.protein2evidence),
             kwargs=dict(tempdir=temp_dir, workers=8),
             name="export-evidences",
             requires=["export-proteins"],
             scheduler=dict(cpu=8, mem=1000, scratch=2000, queue=lsf_queue)),
        Task(fn=uniprot.proteins.export_entry2functions,
             args=(uniprot_url, df.proteins, df.protein2function),
             kwargs=dict(tempdir=temp_dir, workers=8),
             name="export-functions",
             requires=["export-proteins"],
             scheduler=dict(cpu=8, mem=2000, scratch=4000, queue=lsf_queue)),
        Task(fn=uniprot.proteins.export_entry2name,
             args=(uniprot_url, df.proteins, df.protein2name),
             kwargs=dict(tempdir=temp_dir, workers=8),
             name="export-names",
             requires=["export-proteins"],
             scheduler=dict(cpu=8, mem=1000, scratch=4000, queue=lsf_queue)),
        Task(fn=uniprot.proteins.export_entry2proteome,
             args=(uniprot_url, df.proteins, df.protein2proteome),
             kwargs=dict(tempdir=temp_dir, workers=8),
             name="export-proteomes",
             requires=["export-proteins"],
             scheduler=dict(cpu=8, mem=1000, scratch=50000, queue=lsf_queue)),

        # Cross-references
        Task(fn=interpro.xrefs.dump_entries,
             args=(uniprot_url, df.proteins, df.protein2matches,
                   df.protein2proteome, df.protein2domorg, df.structures,
                   config["data"]["metacyc"], df.entryxrefs),
             kwargs=dict(tempdir=temp_dir),
             name="export-xrefs",
             requires=["export-proteomes", "export-dom-orgs",
                       "export-structures"],
             # todo: review
             scheduler=dict(mem=32000, scratch=50000, queue=lsf_queue)),
    ]

    tasks += [
        # Add a "group" task, to include all export tasks
        Task(fn=time.sleep,
             args=(5,),
             name="export",
             requires=get_terminals(tasks)),
    ]

    tasks += [
        Task(fn=interpro.ftp.uniparc.archive_uniparc_matches,
             args=(df.uniparc, df.pub_uniparc),
             name="pub-uniparc",
             requires=["export-uniparc"],
             # todo: review
             scheduler=dict(mem=8000, queue=lsf_queue)),
    ]

    return tasks


def get_terminals(tasks: Sequence[Task],
                  targets: Optional[Sequence[str]] = None) -> List[Task]:
    """Returns a list of terminal/final tasks, i.e. tasks that are not
    dependencies for other tasks.

    :param tasks: A sequence of tasks to evaluate.
    :param targets: An optional sequence of task names.
        If provided, only target tasks thar are terminal nodes are returned.
    :return: A list of tasks.
    """

    # Create a dict of tasks (name -> task)
    tasks = {t.name: t for t in tasks}

    internal_nodes = set()
    for name in (targets or tasks):
        internal_nodes |= traverse_bottom_up(tasks, name)

    terminals = []

    for name in tasks:
        if name in internal_nodes:
            continue
        elif targets and name not in targets:
            continue
        else:
            terminals.append(tasks[name])

    return terminals


def traverse_bottom_up(tasks: Mapping[str, Task], name: str,
                       level: int = 0) -> set:
    internal_nodes = set()

    if level > 0:
        internal_nodes.add(name)

    for parent_name in tasks[name].requires:
        internal_nodes |= traverse_bottom_up(tasks, parent_name, level+1)

    return internal_nodes


def build():
    parser = argparse.ArgumentParser(
        description="Build InterPro7 data warehouse")
    parser.add_argument("config",
                        metavar="config.ini",
                        help="configuration file")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-t", "--tasks",
                       nargs="*",
                       default=None,
                       metavar="TASK",
                       help="tasks to run")
    group.add_argument("-a", "--all-except-completed",
                       action="store_true",
                       help="run all tasks while ignoring completed ones "
                            "(mutually exclusive with -t/--tasks)")
    parser.add_argument("--dry-run",
                        action="store_true",
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
        if args.all_except_completed:
            tasks = workflow.get_remaining_tasks()
        else:
            tasks = args.tasks

        workflow.run(tasks, dry_run=args.dry_run, monitor=not args.detach)
