#!/usr/bin/env python

import argparse
import configparser
import os
import time
from typing import List, Mapping, Optional, Sequence, Set

from mundone import Task, Workflow

from interpro7dw import __version__
from interpro7dw import interpro, pdbe, pfam, uniprot


class DataFiles:
    def __init__(self, root: str, pub_dir: str, create_dirs: bool = True):
        if create_dirs:
            os.makedirs(root, exist_ok=True)
            os.makedirs(pub_dir, mode=0o775, exist_ok=True)

        # Stores
        self.proteins = os.path.join(root, "proteins")
        self.protein2domorg = os.path.join(root, "protein2domorg")
        self.protein2evidence = os.path.join(root, "protein2evidence")
        self.protein2features = os.path.join(root, "protein2features")
        self.protein2functions = os.path.join(root, "protein2functions")
        self.protein2matches = os.path.join(root, "protein2matches")
        self.protein2name = os.path.join(root, "protein2name")
        self.protein2proteome = os.path.join(root, "protein2proteome")
        self.protein2residues = os.path.join(root, "protein2residues")
        self.protein2sequence = os.path.join(root, "protein2sequence")

        # SimpleStores
        self.clans_alignments = os.path.join(root, "clansalignments")
        self.clanxrefs = os.path.join(root, "clanxrefs")
        self.entryxrefs = os.path.join(root, "entryxrefs")
        self.hmms = os.path.join(root, "hmms")
        self.isoforms = os.path.join(root, "isoforms")
        self.pfam_alignments = os.path.join(root, "pfamalignments")
        self.proteomexrefs = os.path.join(root, "proteomexrefs")
        self.structmodels = os.path.join(root, "structmodels")
        self.structurexrefs = os.path.join(root, "structurexrefs")
        self.taxonxrefs = os.path.join(root, "taxonxrefs")
        self.uniparc = os.path.join(root, "uniparc")

        # Data dumps
        self.clans = os.path.join(root, "clans")
        self.databases = os.path.join(root, "databases")
        self.entries = os.path.join(root, "entries")
        self.overlapping_entries = os.path.join(root, "overlapping")
        self.proteomes = os.path.join(root, "proteomes")
        self.structures = os.path.join(root, "structures")
        self.taxa = os.path.join(root, "taxa")

        # Files for FTP
        self.pub_uniparc = os.path.join(pub_dir, "uniparc_match.tar.gz")
        self.relnotes = os.path.join(pub_dir, "announcements.txt")


def gen_tasks(config: configparser.ConfigParser,
              create_dirs: bool = True) -> List[Task]:
    release_version = config["release"]["version"]
    release_date = config["release"]["date"]
    data_dir = config["data"]["path"]
    temp_dir = config["data"]["tmp"]
    ipr_pro_url = config["databases"]["interpro_production"]
    ipr_stg_url = config["databases"]["interpro_staging"]
    ipr_rel_url = config["databases"]["interpro_fallback"]
    goa_url = config["databases"]["goa"]
    intact_url = config["databases"]["intact"]
    pdbe_url = config["databases"]["pdbe"]
    pfam_url = config["databases"]["pfam"]
    uniprot_url = config["databases"]["uniprot"]
    pub_dir = config["exchange"]["interpro"]
    lsf_queue = config["workflow"]["lsf_queue"]

    es_clusters = []
    es_root = os.path.join(data_dir, "elastic")
    es_dirs = [os.path.join(es_root, "default")]
    for cluster, nodes in config.items("elasticsearch"):
        hosts = [host.strip() for host in nodes.split(',') if host.strip()]

        if hosts:
            cluster_dir = os.path.join(es_root, cluster)
            es_clusters.append((cluster, list(set(hosts)), cluster_dir))
            es_dirs.append(cluster_dir)

    df = DataFiles(root=data_dir,
                   pub_dir=os.path.join(pub_dir, release_version),
                   create_dirs=create_dirs)

    tasks = [
        # Data from InterPro (not depending on other tasks)
        Task(fn=interpro.oracle.entries.export_databases,
             args=(ipr_pro_url, release_version, release_date, df.databases),
             kwargs=dict(update=config.getboolean("release", "update")),
             name="export-databases",
             scheduler=dict(mem=100, queue=lsf_queue)),
        Task(fn=interpro.oracle.proteins.export_isoforms,
             args=(ipr_pro_url, df.isoforms),
             name="export-isoforms",
             scheduler=dict(mem=4000, queue=lsf_queue)),
        Task(fn=interpro.oracle.proteins.export_proteins,
             args=(ipr_pro_url, df.proteins),
             kwargs=dict(tempdir=temp_dir),
             name="export-proteins",
             scheduler=dict(mem=4000, scratch=10000, queue=lsf_queue)),
        Task(fn=interpro.oracle.structures.export_structural_models,
             args=(ipr_pro_url, df.structmodels),
             name="export-struct-models",
             scheduler=dict(mem=10000, queue=lsf_queue)),
        Task(fn=interpro.oracle.taxa.export_taxa,
             args=(ipr_pro_url, df.taxa),
             name="export-taxa",
             scheduler=dict(mem=8000, queue=lsf_queue)),
        Task(fn=interpro.oracle.proteins.export_uniparc,
             args=(ipr_pro_url, df.uniparc),
             kwargs=dict(tempdir=temp_dir),
             name="export-uniparc",
             # todo: review
             scheduler=dict(mem=24000, scratch=50000, queue=lsf_queue)),

        # Data from InterPro and/or other sources (Pfam, PDBe)
        Task(fn=interpro.oracle.clans.export_clans,
             args=(ipr_pro_url, pfam_url, df.clans, df.clans_alignments),
             name="export-clans",
             scheduler=dict(mem=2000, queue=lsf_queue)),
        Task(fn=pfam.export_alignments,
             args=(pfam_url, df.pfam_alignments),
             name="export-pfam-alignments",
             scheduler=dict(mem=4000, queue=lsf_queue)),
        Task(fn=pdbe.export_structures,
             args=(ipr_pro_url, pdbe_url, df.structures),
             name="export-structures",
             scheduler=dict(mem=8000, queue=lsf_queue)),

        # Data from InterPro (after export-proteins)
        Task(fn=interpro.oracle.proteins.export_features,
             args=(ipr_pro_url, df.proteins, df.protein2features),
             kwargs=dict(tempdir=temp_dir),
             name="export-features",
             requires=["export-proteins"],
             scheduler=dict(mem=2000, scratch=10000, queue=lsf_queue)),
        Task(fn=interpro.oracle.proteins.export_matches,
             args=(ipr_pro_url, df.proteins, df.protein2matches),
             kwargs=dict(tempdir=temp_dir),
             name="export-matches",
             requires=["export-proteins"],
             scheduler=dict(mem=2000, scratch=50000, queue=lsf_queue)),
        Task(fn=interpro.oracle.proteins.export_residues,
             args=(ipr_pro_url, df.proteins, df.protein2residues),
             kwargs=dict(tempdir=temp_dir),
             name="export-residues",
             requires=["export-proteins"],
             scheduler=dict(mem=8000, scratch=20000, queue=lsf_queue)),
        Task(fn=interpro.oracle.proteins.export_sequences,
             args=(ipr_pro_url, df.proteins, df.protein2sequence),
             kwargs=dict(tempdir=temp_dir),
             name="export-sequences",
             requires=["export-proteins"],
             scheduler=dict(mem=2000, scratch=50000, queue=lsf_queue)),
        Task(fn=interpro.oracle.entries.dump_domain_organisation,
             args=(ipr_pro_url, df.proteins, df.protein2matches,
                   df.protein2domorg),
             kwargs=dict(tempdir=temp_dir),
             name="export-dom-orgs",
             requires=["export-matches"],
             scheduler=dict(mem=4000, scratch=20000, queue=lsf_queue)),
        Task(fn=interpro.oracle.hmms.export_hmms,
             args=(ipr_pro_url, df.protein2matches, df.hmms),
             name="export-hmms",
             requires=["export-matches"],
             scheduler=dict(mem=2000, queue=lsf_queue)),
        Task(fn=interpro.oracle.entries.dump_similar_entries,
             args=(ipr_pro_url, df.protein2matches, df.overlapping_entries),
             name="export-sim-entries",
             requires=["export-matches"],
             scheduler=dict(mem=2000, queue=lsf_queue)),

        # Data from UniProt
        Task(fn=uniprot.proteomes.export_proteomes,
             args=(uniprot_url, df.proteomes),
             name="export-reference-proteomes",
             scheduler=dict(mem=100, queue=lsf_queue)),
        Task(fn=uniprot.proteins.export_entry2evidence,
             args=(uniprot_url, df.proteins, df.protein2evidence),
             kwargs=dict(tempdir=temp_dir),
             name="export-evidences",
             requires=["export-proteins"],
             scheduler=dict(mem=2000, scratch=2000, queue=lsf_queue)),
        Task(fn=uniprot.proteins.export_entry2functions,
             args=(uniprot_url, df.proteins, df.protein2functions),
             kwargs=dict(tempdir=temp_dir),
             name="export-functions",
             requires=["export-proteins"],
             scheduler=dict(mem=2000, scratch=4000, queue=lsf_queue)),
        Task(fn=uniprot.proteins.export_entry2name,
             args=(uniprot_url, df.proteins, df.protein2name),
             kwargs=dict(tempdir=temp_dir),
             name="export-names",
             requires=["export-proteins"],
             scheduler=dict(mem=2000, scratch=4000, queue=lsf_queue)),
        Task(fn=uniprot.proteins.export_entry2proteome,
             args=(uniprot_url, df.proteins, df.protein2proteome),
             kwargs=dict(tempdir=temp_dir),
             name="export-proteomes",
             requires=["export-proteins"],
             scheduler=dict(mem=1000, scratch=1000, queue=lsf_queue)),

        # Exports entry cross-references (e.g entry-proteins, entry-taxa, etc.)
        Task(fn=interpro.xrefs.dump_entries,
             args=(uniprot_url, df.proteins, df.protein2matches,
                   df.protein2proteome, df.protein2domorg,
                   df.structures, df.taxa, config["data"]["metacyc"],
                   config["data"]["alphafold"], df.structmodels,
                   df.entryxrefs),
             kwargs=dict(tempdir=temp_dir),
             name="export-entry2xrefs",
             requires=["export-proteomes", "export-dom-orgs",
                       "export-structures", "export-taxa",
                       "export-struct-models"],
             scheduler=dict(mem=16000, scratch=100000, queue=lsf_queue)),

        # Exports entries (ready to be inserted into MySQL)
        Task(fn=interpro.oracle.entries.export_entries,
             args=(ipr_pro_url, goa_url, intact_url, df.clans,
                   df.overlapping_entries, df.entryxrefs, df.entries),
             kwargs=dict(update=config.getboolean("release", "update")),
             name="export-entries",
             requires=["export-clans", "export-sim-entries",
                       "export-entry2xrefs"],
             scheduler=dict(mem=4000, queue=lsf_queue)),

        # Exports cross-references for other entities (needed for counters)
        Task(fn=interpro.xrefs.dump_clans,
             args=(df.clans, df.proteins, df.protein2matches,
                   df.protein2proteome, df.protein2domorg, df.structures,
                   df.clanxrefs),
             kwargs=dict(tempdir=temp_dir),
             name="export-clan2xrefs",
             requires=["export-clans", "export-proteomes", "export-dom-orgs",
                       "export-structures"],
             scheduler=dict(mem=8000, scratch=20000, queue=lsf_queue)),
        Task(fn=interpro.xrefs.dump_proteomes,
             args=(df.proteins, df.protein2matches, df.protein2proteome,
                   df.protein2domorg, df.structures, df.entries,
                   df.proteomes, df.proteomexrefs),
             kwargs=dict(tempdir=temp_dir),
             name="export-proteome2xrefs",
             requires=["export-entries", "export-reference-proteomes"],
             # todo: review
             scheduler=dict(mem=16000, scratch=50000, queue=lsf_queue)),
        Task(fn=interpro.xrefs.dump_structures,
             args=(df.proteins, df.protein2matches, df.protein2proteome,
                   df.protein2domorg, df.structures, df.entries,
                   df.structurexrefs),
             name="export-structure2xrefs",
             requires=["export-entries"],
             # todo: review
             scheduler=dict(mem=16000, queue=lsf_queue)),
        Task(fn=interpro.xrefs.dump_taxa,
             args=(df.proteins, df.protein2matches, df.protein2proteome,
                   df.structures, df.entries, df.taxa, df.taxonxrefs),
             kwargs=dict(tempdir=temp_dir),
             name="export-taxon2xrefs",
             requires=["export-entries"],
             # todo: review
             scheduler=dict(mem=16000, scratch=50000, queue=lsf_queue)),
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

    insert_tasks = [
        Task(fn=interpro.mysql.entries.insert_annotations,
             args=(ipr_stg_url, df.hmms, df.pfam_alignments),
             name="insert-annotations",
             requires=["export-hmms", "export-pfam-alignments"],
             # todo: review
             scheduler=dict(mem=16000, queue=lsf_queue)),
        Task(fn=interpro.mysql.clans.insert_clans,
             args=(ipr_stg_url, df.clans, df.clanxrefs, df.clans_alignments),
             name="insert-clans",
             requires=["export-clan2xrefs"],
             scheduler=dict(mem=1000, queue=lsf_queue)),
        Task(fn=interpro.mysql.entries.insert_databases,
             args=(ipr_stg_url, df.databases),
             name="insert-databases",
             requires=["export-databases"],
             # todo: review
             scheduler=dict(mem=1000, queue=lsf_queue)),
        Task(fn=interpro.mysql.entries.insert_entries,
             args=(ipr_stg_url, pfam_url, df.entries, df.entryxrefs),
             name="insert-entries",
             requires=["export-entries"],
             # todo: review
             scheduler=dict(mem=16000, queue=lsf_queue)),
        Task(fn=interpro.mysql.proteins.insert_isoforms,
             args=(ipr_stg_url, df.entries, df.isoforms),
             name="insert-isoforms",
             requires=["export-entries", "export-isoforms"],
             # todo: review
             scheduler=dict(mem=4000, queue=lsf_queue)),
        Task(fn=interpro.mysql.proteins.insert_proteins,
             args=(ipr_stg_url, pdbe_url, df.entries, df.isoforms,
                   df.structures, df.taxa, df.proteins, df.protein2domorg,
                   df.protein2evidence, df.protein2functions,
                   df.protein2matches, df.protein2name, df.protein2proteome,
                   df.protein2sequence),
             name="insert-proteins",
             requires=["export-entries", "export-isoforms", "export-evidences",
                       "export-functions", "export-names", "export-sequences"],
             # todo: review
             scheduler=dict(mem=16000, queue=lsf_queue)),
        Task(fn=interpro.mysql.proteins.insert_protein_features,
             args=(ipr_stg_url, df.protein2features),
             name="insert-protein-features",
             requires=["export-features"],
             # todo: review
             scheduler=dict(mem=16000, queue=lsf_queue)),
        Task(fn=interpro.mysql.proteins.insert_protein_residues,
             args=(ipr_stg_url, df.protein2residues),
             name="insert-protein-residues",
             requires=["export-residues"],
             # todo: review
             scheduler=dict(mem=16000, queue=lsf_queue)),
        Task(fn=interpro.mysql.proteomes.insert_proteomes,
             args=(ipr_stg_url, df.proteomes, df.proteomexrefs),
             name="insert-proteomes",
             requires=["export-proteome2xrefs"],
             # todo: review
             scheduler=dict(mem=16000, queue=lsf_queue)),
        Task(fn=interpro.mysql.entries.insert_release_notes,
             args=(ipr_stg_url, ipr_rel_url, df.entries, df.proteomes,
                   df.structures, df.taxa, df.proteins, df.protein2matches,
                   df.protein2proteome, df.relnotes),
             name="insert-release-notes",
             # todo: review
             scheduler=dict(mem=12000, queue=lsf_queue),
             requires=["export-entries", "export-reference-proteomes",
                       "insert-databases"]),
        Task(fn=interpro.mysql.structures.insert_structures,
             args=(ipr_stg_url, df.structures, df.structurexrefs),
             name="insert-structures",
             requires=["export-structure2xrefs"],
             # todo: review
             scheduler=dict(mem=16000, queue=lsf_queue)),
        Task(fn=interpro.mysql.structures.insert_structural_models,
             args=(ipr_stg_url, df.entries, df.structmodels),
             name="insert-struct-models",
             requires=["export-entries", "export-struct-models"],
             # todo: review
             scheduler=dict(mem=8000, queue=lsf_queue)),
        Task(fn=interpro.mysql.taxa.insert_taxa,
             args=(ipr_stg_url, df.entries, df.taxa, df.taxonxrefs),
             name="insert-taxa",
             requires=["export-taxon2xrefs"],
             # todo: review
             scheduler=dict(mem=16000, queue=lsf_queue)),
    ]

    tasks += insert_tasks
    tasks += [
        Task(fn=time.sleep,
             args=(5,),
             name="insert",
             requires=get_terminals(tasks, [t.name for t in insert_tasks])),
    ]

    es_tasks = [
        Task(fn=interpro.elastic.export_documents,
             args=(df.proteins, df.protein2matches, df.protein2domorg,
                   df.protein2proteome, df.entries, df.proteomes,
                   df.structures, df.taxa, config["data"]["alphafold"],
                   es_dirs, release_version),
             name="es-export",
             requires=["export-entries", "export-reference-proteomes"],
             # todo: review
             scheduler=dict(mem=24000, queue=lsf_queue))
    ]

    for cluster, hosts, cluster_dir in es_clusters:
        es_tasks += [
            Task(
                fn=interpro.elastic.create_indices,
                args=(df.databases, hosts, release_version),
                name=f"es-init-{cluster}",
                scheduler=dict(mem=100, queue=lsf_queue),
                requires=["export-databases"] + list(es_tasks[0].requires)
            ),
            Task(
                fn=interpro.elastic.index_documents,
                args=(hosts, cluster_dir, release_version),
                kwargs=dict(threads=8),
                name=f"es-index-{cluster}",
                # todo: review
                scheduler=dict(mem=16000, queue=lsf_queue),
                requires=[f"es-init-{cluster}"]
            )
        ]

    tasks += es_tasks
    tasks += [
        Task(fn=time.sleep,
             args=(5,),
             name="elastic",
             requires=get_terminals(tasks, [t.name for t in es_tasks])),
    ]

    for cluster, hosts, cluster_dir in es_clusters:
        tasks += [
            Task(
                fn=interpro.elastic.publish,
                args=(hosts,),
                name=f"es-publish-{cluster}",
                scheduler=dict(mem=100, queue=lsf_queue),
                requires=["es-export", f"es-index-{cluster}"]
            )
        ]

    return tasks


def clean_deps(task: Task, tasks: Sequence[Task]) -> Set[str]:
    tasks = {t.name: t for t in tasks}

    direct_deps = set()
    all_deps = set()
    for parent_name in task.requires:
        direct_deps.add(parent_name)
        all_deps |= traverse_bottom_up(tasks, parent_name)

    return direct_deps - all_deps


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

    tasks = gen_tasks(config, create_dirs=not args.dry_run)

    database = os.path.join(workflow_dir, f"{version}.sqlite")
    with Workflow(tasks, dir=workflow_dir, database=database) as workflow:
        if args.all_except_completed:
            tasks = workflow.get_remaining_tasks()
        else:
            tasks = args.tasks

        workflow.run(tasks, dry_run=args.dry_run, monitor=not args.detach)
