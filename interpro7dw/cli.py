#!/usr/bin/env python

import argparse
import configparser
import os
import time
from typing import Optional

from mundone import Task, Workflow

from interpro7dw import __version__
from interpro7dw import alphafold, interpro, pdbe, pfam, uniprot


def wait(secs: int = 5):
    time.sleep(secs)


class DataFiles:
    def __init__(self, root: str):
        # BasicStores
        self.clans_alignments = os.path.join(root, "clan-alignments")
        self.entry2xrefs = os.path.join(root, "entry-xrefs")
        self.hmms = os.path.join(root, "hmms")
        self.isoforms = os.path.join(root, "isoforms")
        self.pdbematches = os.path.join(root, "pdbe-matches")
        self.pfam_alignments = os.path.join(root, "pfam-alignments")
        self.protein2features = os.path.join(root, "protein-features")
        self.protein2residues = os.path.join(root, "protein-residues")
        self.protein2uniparc = os.path.join(root, "protein-uniparc")
        self.structmodels = os.path.join(root, "structmodels")

        # KVStores
        self.proteins = os.path.join(root, "proteins")
        self.protein2alphafold = os.path.join(root, "protein-alphafold")
        self.protein2domorg = os.path.join(root, "protein-domorg")
        self.protein2evidence = os.path.join(root, "protein-evidence")
        self.protein2functions = os.path.join(root, "protein-functions")
        self.protein2matches = os.path.join(root, "protein-matches")
        self.protein2name = os.path.join(root, "protein-name")
        self.protein2proteome = os.path.join(root, "protein-proteome")
        self.protein2sequence = os.path.join(root, "protein-sequence")

        # Pickles
        self.clans = os.path.join(root, "clans")
        self.databases = os.path.join(root, "databases")
        self.overlapping = os.path.join(root, "overlapping")
        self.proteomes = os.path.join(root, "proteomes")
        self.structures = os.path.join(root, "structures")
        self.taxa = os.path.join(root, "taxa")

        # self.proteins = os.path.join(root, "proteins")
        # self.protein2domorg = os.path.join(root, "protein2domorg")
        # self.protein2evidence = os.path.join(root, "protein2evidence")
        # self.protein2features = os.path.join(root, "protein2features")
        # self.protein2functions = os.path.join(root, "protein2functions")
        # self.protein2matches = os.path.join(root, "protein2matches")
        # self.protein2name = os.path.join(root, "protein2name")
        # self.protein2proteome = os.path.join(root, "protein2proteome")
        # self.protein2sequence = os.path.join(root, "protein2sequence")
        #
        # # SimpleStores
        # self.clanxrefs = os.path.join(root, "clanxrefs")
        # self.entryxrefs = os.path.join(root, "entryxrefs")

        # self.proteomexrefs = os.path.join(root, "proteomexrefs")
        # self.structmodels = os.path.join(root, "structmodels")
        # self.structurexrefs = os.path.join(root, "structurexrefs")
        # self.taxonxrefs = os.path.join(root, "taxonxrefs")
        # self.uniparc = os.path.join(root, "uniparc")
        #
        # # Data dumps
        # self.databases = os.path.join(root, "databases")
        # self.entries = os.path.join(root, "entries")
        # self.overlapping_entries = os.path.join(root, "overlapping")
        # self.proteomes = os.path.join(root, "proteomes")
        # self.structures = os.path.join(root, "structures")
        # self.taxa = os.path.join(root, "taxa")


def gen_tasks(config: configparser.ConfigParser) -> list[Task]:
    release_version = config["release"]["version"]
    release_date = config["release"]["date"]
    update_db = config.getboolean("release", "update")
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
    pub_dir = os.path.join(config["exchange"]["interpro"], release_version)
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

    df = DataFiles(data_dir)

    tasks = [
        # Exports without dependencies
        Task(fn=interpro.oracle.clans.export_clans,
             args=(ipr_pro_url, pfam_url, df.clans, df.clans_alignments),
             name="export-clans",
             scheduler=dict(mem=2000, queue=lsf_queue)),
        Task(fn=interpro.oracle.databases.export,
             args=(ipr_pro_url, release_version, release_date, df.databases),
             kwargs=dict(update=update_db),
             name="export-databases",
             scheduler=dict(mem=100, queue=lsf_queue)),
        Task(fn=interpro.oracle.proteins.export_isoforms,
             args=(ipr_pro_url, df.isoforms),
             name="export-isoforms",
             scheduler=dict(mem=4000, queue=lsf_queue)),
        Task(fn=interpro.oracle.proteins.export_features,
             args=(ipr_pro_url, df.protein2features),
             name="export-features",
             scheduler=dict(mem=4000, queue=lsf_queue)),
        Task(fn=interpro.oracle.proteins.export_matches,
             args=(ipr_pro_url, df.protein2matches),
             name="export-matches",
             scheduler=dict(mem=4000, queue=lsf_queue)),
        Task(fn=interpro.oracle.proteins.export_proteins,
             args=(ipr_pro_url, df.proteins),
             name="export-proteins",
             scheduler=dict(mem=4000, queue=lsf_queue)),
        Task(fn=interpro.oracle.proteins.export_residues,
             args=(ipr_pro_url, df.protein2residues),
             name="export-residues",
             scheduler=dict(mem=4000, queue=lsf_queue)),
        Task(fn=interpro.oracle.proteins.export_uniparc,
             args=(ipr_pro_url, df.protein2uniparc),
             name="export-uniparc",
             scheduler=dict(mem=4000, queue=lsf_queue)),
        Task(fn=interpro.oracle.structures.export_pdbe_matches,
             args=(ipr_pro_url, df.pdbematches),
             name="export-pdbe-matches",
             scheduler=dict(mem=16000, queue=lsf_queue)),
        Task(fn=interpro.oracle.structures.export_structural_models,
             args=(ipr_pro_url, df.structmodels),
             name="export-struct-models",
             scheduler=dict(mem=16000, queue=lsf_queue)),
        Task(fn=interpro.oracle.taxa.export_taxa,
             args=(ipr_pro_url, df.taxa),
             name="export-taxa",
             scheduler=dict(mem=8000, queue=lsf_queue)),
        Task(fn=pdbe.export_structures,
             args=(ipr_pro_url, pdbe_url, df.structures),
             name="export-structures",
             scheduler=dict(mem=8000, queue=lsf_queue)),
        Task(fn=pfam.export_alignments,
             args=(pfam_url, df.pfam_alignments),
             name="export-pfam-alignments",
             scheduler=dict(mem=4000, queue=lsf_queue)),

        # Exports with dependencies
        Task(fn=alphafold.export,
             args=(config["data"]["alphafold"], df.proteins,
                   df.protein2alphafold),
             kwargs=dict(tempdir=temp_dir),
             name="export-alphafold",
             requires=["export-proteins"],
             scheduler=dict(mem=2000, tmp=50000, queue=lsf_queue)),
        Task(fn=interpro.oracle.hmms.export_hmms,
             args=(ipr_pro_url, df.protein2matches, df.hmms),
             name="export-hmms",
             requires=["export-matches"],
             scheduler=dict(mem=2000, queue=lsf_queue)),
        Task(fn=interpro.oracle.proteins.export_sequences,
             args=(ipr_pro_url, df.proteins, df.protein2sequence),
             kwargs=dict(tempdir=temp_dir),
             name="export-sequences",
             requires=["export-proteins"],
             scheduler=dict(mem=2000, tmp=50000, queue=lsf_queue)),
        Task(fn=uniprot.proteomes.export_proteomes,
             args=(uniprot_url, df.proteomes),
             name="export-reference-proteomes",
             scheduler=dict(mem=100, queue=lsf_queue)),
        Task(fn=uniprot.proteins.export_entry2evidence,
             args=(uniprot_url, df.proteins, df.protein2evidence),
             kwargs=dict(tempdir=temp_dir),
             name="export-evidences",
             requires=["export-proteins"],
             scheduler=dict(mem=4000, tmp=2000, queue=lsf_queue)),
        Task(fn=uniprot.proteins.export_entry2functions,
             args=(uniprot_url, df.proteins, df.protein2functions),
             kwargs=dict(tempdir=temp_dir),
             name="export-functions",
             requires=["export-proteins"],
             scheduler=dict(mem=2000, tmp=4000, queue=lsf_queue)),
        Task(fn=uniprot.proteins.export_entry2name,
             args=(uniprot_url, df.proteins, df.protein2name),
             kwargs=dict(tempdir=temp_dir),
             name="export-names",
             requires=["export-proteins"],
             scheduler=dict(mem=8000, tmp=4000, queue=lsf_queue)),
        Task(fn=uniprot.proteins.export_entry2proteome,
             args=(uniprot_url, df.proteins, df.protein2proteome),
             kwargs=dict(tempdir=temp_dir),
             name="export-proteomes",
             requires=["export-proteins"],
             scheduler=dict(mem=1000, tmp=1000, queue=lsf_queue)),
    ]

    tasks += [
        # Add a "group" task, to include all export tasks
        Task(fn=wait,
             name="export",
             requires=get_terminals(tasks)),
    ]

    tasks += [
        Task(fn=interpro.xrefs.domorgs.export,
             args=(df.proteins, df.protein2matches, df.protein2domorg),
             name="export-dom-orgs",
             requires=["export-matches"],
             scheduler=dict(mem=8000, queue=lsf_queue)),
        Task(fn=interpro.xrefs.entries.export_sim_entries,
             args=(df.protein2matches, df.overlapping),
             name="export-sim-entries",
             requires=["export-matches"],
             scheduler=dict(mem=2000, queue=lsf_queue)),
        Task(fn=interpro.xrefs.entries.export_xrefs,
             args=(uniprot_url, df.proteins, df.protein2matches,
                   df.protein2alphafold, df.protein2proteome,
                   df.protein2domorg, df.structmodels, df.structures,
                   df.taxa, config["data"]["metacyc"], df.entry2xrefs),
             kwargs=dict(interpro_uri=ipr_pro_url if update_db else None,
                         tempdir=temp_dir),
             name="export-entry2xrefs",
             requires=["export-proteomes", "export-dom-orgs",
                       "export-structures", "export-taxa",
                       "export-struct-models"],
             scheduler=dict(mem=16000, tmp=120000, queue=lsf_queue)),
    ]

    # tasks = [
    #
    #
    #
    #     # Exports entry cross-references (e.g entry-proteins, entry-taxa, etc.)
    #
    #
    #     # Exports entries (ready to be inserted into MySQL)
    #     Task(fn=interpro.oracle.entries.export_entries,
    #          args=(ipr_pro_url, goa_url, intact_url, df.clans,
    #                df.overlapping_entries, df.entryxrefs, df.entries),
    #          kwargs=dict(update=config.getboolean("release", "update")),
    #          name="export-entries",
    #          requires=["export-clans", "export-sim-entries",
    #                    "export-entry2xrefs"],
    #          scheduler=dict(mem=8000, queue=lsf_queue)),
    #
    #     # Exports cross-references for other entities (needed for counters)
    #     Task(fn=interpro.xrefs.dump_clans,
    #          args=(df.clans, df.proteins, df.protein2matches,
    #                df.protein2proteome, df.protein2domorg, df.structures,
    #                df.clanxrefs),
    #          kwargs=dict(tempdir=temp_dir),
    #          name="export-clan2xrefs",
    #          requires=["export-clans", "export-proteomes", "export-dom-orgs",
    #                    "export-structures"],
    #          scheduler=dict(mem=8000, tmp=20000, queue=lsf_queue)),
    #     Task(fn=interpro.xrefs.dump_proteomes,
    #          args=(df.proteins, df.protein2matches, df.protein2proteome,
    #                df.protein2domorg, df.structures, df.entries,
    #                df.proteomes, df.proteomexrefs),
    #          kwargs=dict(tempdir=temp_dir),
    #          name="export-proteome2xrefs",
    #          requires=["export-entries", "export-reference-proteomes"],
    #          scheduler=dict(mem=8000, tmp=5000, queue=lsf_queue)),
    #     Task(fn=interpro.xrefs.dump_structures,
    #          args=(df.proteins, df.protein2matches, df.protein2proteome,
    #                df.protein2domorg, df.structures, df.entries,
    #                df.structurexrefs),
    #          name="export-structure2xrefs",
    #          requires=["export-entries"],
    #          scheduler=dict(mem=8000, queue=lsf_queue)),
    #     Task(fn=interpro.xrefs.dump_taxa,
    #          args=(df.proteins, df.protein2matches, df.protein2proteome,
    #                df.structures, df.entries, df.taxa, df.taxonxrefs),
    #          kwargs=dict(tempdir=temp_dir),
    #          name="export-taxon2xrefs",
    #          requires=["export-entries"],
    #          scheduler=dict(mem=12000, tmp=60000, queue=lsf_queue)),
    #
    #     # UniParc matches (for FTP)
    #     Task(fn=interpro.oracle.proteins.export_uniparc,
    #          args=(ipr_pro_url, df.entries, df.uniparc),
    #          kwargs=dict(tempdir=temp_dir),
    #          name="export-uniparc",
    #          requires=["export-entries"],
    #          # TODO: review
    #          scheduler=dict(mem=16000, tmp=70000, queue=lsf_queue)),
    # ]
    #
    # tasks += [
    #     # Add a "group" task, to include all export tasks
    #     Task(fn=wait,
    #          name="export",
    #          requires=get_terminals(tasks)),
    # ]
    #
    # insert_tasks = [
    #     Task(fn=interpro.mysql.entries.insert_annotations,
    #          args=(ipr_stg_url, df.hmms, df.pfam_alignments),
    #          name="insert-annotations",
    #          requires=["export-hmms", "export-pfam-alignments"],
    #          scheduler=dict(mem=4000, queue=lsf_queue)),
    #     Task(fn=interpro.mysql.clans.insert_clans,
    #          args=(ipr_stg_url, df.clans, df.clanxrefs, df.clans_alignments),
    #          name="insert-clans",
    #          requires=["export-clan2xrefs"],
    #          scheduler=dict(mem=2000, queue=lsf_queue)),
    #     Task(fn=interpro.mysql.entries.insert_databases,
    #          args=(ipr_stg_url, df.databases),
    #          name="insert-databases",
    #          requires=["export-databases"],
    #          scheduler=dict(mem=1000, queue=lsf_queue)),
    #     Task(fn=interpro.mysql.entries.insert_entries,
    #          args=(ipr_stg_url, pfam_url, df.entries, df.entryxrefs),
    #          name="insert-entries",
    #          requires=["export-entries"],
    #          scheduler=dict(mem=12000, queue=lsf_queue)),
    #     Task(fn=interpro.mysql.proteins.insert_isoforms,
    #          args=(ipr_stg_url, df.entries, df.isoforms),
    #          name="insert-isoforms",
    #          requires=["export-entries", "export-isoforms"],
    #          scheduler=dict(mem=4000, queue=lsf_queue)),
    #     Task(fn=interpro.mysql.proteins.insert_proteins,
    #          args=(ipr_stg_url, pdbe_url, df.entries, df.isoforms,
    #                df.structures, df.taxa, df.proteins, df.protein2domorg,
    #                df.protein2evidence, df.protein2functions,
    #                df.protein2matches, df.protein2name, df.protein2proteome,
    #                df.protein2sequence),
    #          name="insert-proteins",
    #          requires=["export-entries", "export-isoforms", "export-evidences",
    #                    "export-functions", "export-names", "export-sequences"],
    #          scheduler=dict(mem=8000, queue=lsf_queue)),
    #     Task(fn=interpro.mysql.proteins.insert_protein_features,
    #          args=(ipr_stg_url, df.protein2features),
    #          name="insert-protein-features",
    #          requires=["export-features"],
    #          scheduler=dict(mem=1000, queue=lsf_queue)),
    #     Task(fn=interpro.mysql.proteins.insert_protein_residues,
    #          args=(ipr_stg_url, df.protein2residues),
    #          name="insert-protein-residues",
    #          requires=["export-residues"],
    #          scheduler=dict(mem=1000, queue=lsf_queue)),
    #     Task(fn=interpro.mysql.proteomes.insert_proteomes,
    #          args=(ipr_stg_url, df.proteomes, df.proteomexrefs),
    #          name="insert-proteomes",
    #          requires=["export-proteome2xrefs"],
    #          scheduler=dict(mem=1000, queue=lsf_queue)),
    #     Task(fn=interpro.mysql.entries.insert_release_notes,
    #          args=(ipr_stg_url, ipr_rel_url, df.entries, df.proteomes,
    #                df.structures, df.taxa, df.proteins, df.protein2matches,
    #                df.protein2proteome),
    #          name="insert-release-notes",
    #          scheduler=dict(mem=12000, queue=lsf_queue),
    #          requires=["export-entries", "export-reference-proteomes",
    #                    "insert-databases"]),
    #     Task(fn=interpro.mysql.structures.insert_structures,
    #          args=(ipr_stg_url, df.structures, df.structurexrefs),
    #          name="insert-structures",
    #          requires=["export-structure2xrefs"],
    #          scheduler=dict(mem=8000, queue=lsf_queue)),
    #     Task(fn=interpro.mysql.structures.insert_structural_models,
    #          args=(ipr_stg_url, df.entries, df.structmodels),
    #          name="insert-struct-models",
    #          requires=["export-entries", "export-struct-models"],
    #          scheduler=dict(mem=12000, queue=lsf_queue)),
    #     Task(fn=interpro.mysql.taxa.insert_taxa,
    #          args=(ipr_stg_url, df.entries, df.taxa, df.taxonxrefs),
    #          name="insert-taxa",
    #          requires=["export-taxon2xrefs"],
    #          scheduler=dict(mem=4000, queue=lsf_queue)),
    # ]
    #
    # tasks += insert_tasks
    # tasks += [
    #     Task(fn=wait,
    #          name="insert",
    #          requires=get_terminals(tasks, [t.name for t in insert_tasks])),
    # ]
    #
    # es_tasks = [
    #     Task(fn=interpro.elastic.export_documents,
    #          args=(df.proteins, df.protein2matches, df.protein2domorg,
    #                df.protein2proteome, df.entries, df.proteomes,
    #                df.structures, df.taxa, config["data"]["alphafold"],
    #                es_dirs, release_version),
    #          name="es-export",
    #          requires=["export-entries", "export-reference-proteomes"],
    #          scheduler=dict(mem=16000, queue=lsf_queue))
    # ]
    #
    # for cluster, hosts, cluster_dir in es_clusters:
    #     es_tasks += [
    #         Task(
    #             fn=interpro.elastic.create_indices,
    #             args=(df.databases, hosts, release_version),
    #             name=f"es-init-{cluster}",
    #             scheduler=dict(mem=100, queue=lsf_queue),
    #             requires=["export-databases"] + list(es_tasks[0].requires)
    #         ),
    #         Task(
    #             fn=interpro.elastic.index_documents,
    #             args=(hosts, cluster_dir, release_version),
    #             kwargs=dict(threads=8),
    #             name=f"es-index-{cluster}",
    #             # todo: review
    #             scheduler=dict(mem=16000, queue=lsf_queue),
    #             requires=[f"es-init-{cluster}"]
    #         )
    #     ]
    #
    # tasks += es_tasks
    # tasks += [
    #     Task(fn=wait,
    #          name="elastic",
    #          requires=get_terminals(tasks, [t.name for t in es_tasks])),
    # ]
    #
    # for cluster, hosts, cluster_dir in es_clusters:
    #     tasks += [
    #         Task(
    #             fn=interpro.elastic.publish,
    #             args=(hosts,),
    #             name=f"es-publish-{cluster}",
    #             scheduler=dict(mem=100, queue=lsf_queue),
    #             requires=["es-export", f"es-index-{cluster}"]
    #         )
    #     ]
    #
    # # Tasks for files to distribute to FTP
    # tasks += [
    #     Task(fn=interpro.ftp.flatfiles.export,
    #          args=(df.entries, df.protein2matches, pub_dir),
    #          name="ftp-flatfiles",
    #          requires=["export-entries"],
    #          # todo: review
    #          scheduler=dict(mem=16000, queue=lsf_queue)),
    #     Task(fn=interpro.ftp.relnotes.export,
    #          args=(ipr_stg_url, pub_dir),
    #          name="ftp-relnotes",
    #          requires=["insert-release-notes"],
    #          # todo: review
    #          scheduler=dict(mem=16000, queue=lsf_queue)),
    #     Task(fn=interpro.ftp.uniparc.archive_uniparc_matches,
    #          args=(df.uniparc, pub_dir),
    #          name="ftp-uniparc",
    #          requires=["export-uniparc"],
    #          # todo: review
    #          scheduler=dict(mem=8000, queue=lsf_queue)),
    #
    #     Task(fn=interpro.ftp.xmlfiles.export_interpro,
    #          args=(df.entries, df.entryxrefs, df.databases, df.taxa, pub_dir),
    #          name="ftp-interpro",
    #          requires=["export-entries", "export-databases"],
    #          # todo: review
    #          scheduler=dict(mem=16000, queue=lsf_queue)),
    #     Task(fn=interpro.ftp.xmlfiles.export_feature_matches,
    #          args=(df.databases, df.proteins, df.protein2features, pub_dir),
    #          name="ftp-features",
    #          requires=["export-databases", "export-features"],
    #          # todo: review
    #          scheduler=dict(mem=16000, queue=lsf_queue)),
    #     Task(fn=interpro.ftp.xmlfiles.export_matches,
    #          args=(df.databases, df.entries, df.isoforms, df.proteins,
    #                df.protein2matches, pub_dir),
    #          name="ftp-matches",
    #          requires=["export-databases", "export-entries",
    #                    "export-isoforms"],
    #          # todo: review
    #          scheduler=dict(mem=16000, queue=lsf_queue)),
    #     Task(
    #         fn=interpro.ftp.xmlfiles.export_structure_matches,
    #         args=(pdbe_url, df.proteins, df.structures, pub_dir),
    #         name="ftp-structures",
    #         # todo: review
    #         scheduler=dict(mem=8000, queue=lsf_queue),
    #         requires=["export-proteins", "export-structures"]
    #     ),
    #     Task(
    #         fn=wait,
    #         name="ftp",
    #         scheduler=dict(queue=lsf_queue),
    #         requires=["ftp-flatfiles", "ftp-relnotes", "ftp-uniparc",
    #                   "ftp-interpro", "ftp-features", "ftp-matches",
    #                   "ftp-structures"]
    #     )
    # ]
    #
    # # Tasks for other EMBL-EBI services
    # tasks += [
    #     Task(
    #         fn=ebisearch.export,
    #         args=(ipr_stg_url, df.entries, df.entryxrefs, df.taxa,
    #               os.path.join(data_dir, "ebisearch")),
    #         name="export-ebisearch",
    #         scheduler=dict(mem=12000, queue=lsf_queue),
    #         requires=["insert-databases", "export-entries"]
    #     ),
    #     Task(
    #         fn=ebisearch.publish,
    #         args=(os.path.join(data_dir, "ebisearch"),
    #               config["exchange"]["ebisearch"]),
    #         name="publish-ebisearch",
    #         scheduler=dict(queue=lsf_queue),
    #         requires=["export-ebisearch"]
    #     ),
    #     Task(
    #         fn=uniprot.goa.export,
    #         args=(ipr_pro_url, ipr_stg_url, pdbe_url, df.entries,
    #               df.entryxrefs, os.path.join(data_dir, "goa")),
    #         name="export-goa",
    #         # todo: review
    #         scheduler=dict(mem=12000, queue=lsf_queue),
    #         requires=["insert-databases", "export-entries"]
    #     ),
    #     Task(
    #         fn=uniprot.goa.publish,
    #         args=(os.path.join(data_dir, "goa"),
    #               config["exchange"]["goa"]),
    #         name="publish-goa",
    #         scheduler=dict(queue=lsf_queue),
    #         requires=["export-goa"]
    #     ),
    #     Task(
    #         fn=pdbe.export_pdb_matches,
    #         args=(ipr_pro_url, ipr_stg_url, df.entries,
    #               os.path.join(data_dir, "pdbe")),
    #         name="export-pdbe",
    #         # todo: review
    #         scheduler=dict(mem=12000, queue=lsf_queue),
    #         requires=["insert-databases", "export-entries"]
    #     ),
    #     Task(
    #         fn=pdbe.publish,
    #         args=(os.path.join(data_dir, "pdbe"), config["exchange"]["pdbe"]),
    #         name="publish-pdbe",
    #         scheduler=dict(queue=lsf_queue),
    #         requires=["export-pdbe"]
    #     ),
    #     Task(fn=wait,
    #          name="ebi-services",
    #          requires=["export-ebisearch", "export-goa", "export-pdbe"]),
    # ]

    return tasks


def clean_deps(task: Task, tasks: list[Task]) -> set[str]:
    tasks = {t.name: t for t in tasks}

    direct_deps = set()
    all_deps = set()
    for parent_name in task.requires:
        direct_deps.add(parent_name)
        all_deps |= traverse_bottom_up(tasks, parent_name)

    return direct_deps - all_deps


def get_terminals(tasks: list[Task],
                  targets: Optional[list[str]] = None) -> list[Task]:
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


def traverse_bottom_up(tasks: dict[str, Task], name: str,
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

    print(f"dropping database: {args.database}")
    uri = config["databases"][f"interpro_{args.database}"]
    interpro.mysql.utils.drop_database(uri)
    print("done")
