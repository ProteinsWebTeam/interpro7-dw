#!/usr/bin/env python

import argparse
import configparser
import os
import time
from typing import Optional

from mundone import Task, Workflow

from interpro7dw import __version__
from interpro7dw import alphafold, ebisearch, interpro, pdbe, pfam, uniprot


def wait(secs: int = 5):
    time.sleep(secs)


class DataFiles:
    def __init__(self, root: str):
        # BasicStores
        self.clan2xrefs = os.path.join(root, "clan-xrefs")
        self.clans_alignments = os.path.join(root, "clan-alignments")
        self.entry2xrefs = os.path.join(root, "entry-xrefs")
        self.hmms = os.path.join(root, "hmms")
        self.isoforms = os.path.join(root, "isoforms")
        self.pdbematches = os.path.join(root, "pdbe-matches")
        self.pfam_alignments = os.path.join(root, "pfam-alignments")
        self.protein2features = os.path.join(root, "protein-features")
        self.protein2residues = os.path.join(root, "protein-residues")
        self.proteome2xrefs = os.path.join(root, "proteome-xrefs")
        self.rosettafold = os.path.join(root, "structmodels")
        self.structure2xrefs = os.path.join(root, "structure-xfres")
        self.taxon2xrefs = os.path.join(root, "taxon-xrefs")

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
        self.uniparcmatches = os.path.join(root, "uniparc-matches")
        self.uniparcproteins = os.path.join(root, "uniparc-proteins")

        # Pickles
        self.clans = os.path.join(root, "clans")
        self.databases = os.path.join(root, "databases")
        self.entries = os.path.join(root, "entries")
        self.overlapping = os.path.join(root, "overlapping")
        self.proteomes = os.path.join(root, "proteomes")
        self.protein2structures = os.path.join(root, "protein-structures")
        self.structures = os.path.join(root, "structures")
        self.taxa = os.path.join(root, "taxa")


def gen_tasks(config: configparser.ConfigParser) -> list[Task]:
    release_version = config["release"]["version"]
    release_date = config["release"]["date"]
    update_db = config.getboolean("release", "update")
    data_dir = config["data"]["path"]
    temp_dir = config["data"]["tmp"]
    ipr_pro_uri = config["databases"]["interpro_production"]
    ipr_stg_uri = config["databases"]["interpro_staging"]
    ipr_rel_uri = config["databases"]["interpro_fallback"]
    ips_pro_uri = config["databases"]["iprscan_production"]
    goa_uri = config["databases"]["goa"]
    intact_uri = config["databases"]["intact"]
    pdbe_uri = config["databases"]["pdbe"]
    pfam_uri = config["databases"]["pfam"]
    uniprot_uri = config["databases"]["uniprot"]
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
             args=(ipr_pro_uri, pfam_uri, df.clans, df.clans_alignments),
             name="export-clans",
             scheduler=dict(mem=2000, queue=lsf_queue)),
        Task(fn=interpro.oracle.databases.export,
             args=(ipr_pro_uri, release_version, release_date, df.databases),
             kwargs=dict(update=update_db),
             name="export-databases",
             scheduler=dict(mem=100, queue=lsf_queue)),
        Task(fn=interpro.oracle.entries.export_entries,
             args=(ipr_pro_uri, goa_uri, intact_uri, df.entries),
             name="export-entries",
             scheduler=dict(mem=3000, queue=lsf_queue)),
        Task(fn=interpro.oracle.matches.export_isoforms,
             args=(ipr_pro_uri, df.isoforms),
             name="export-isoforms",
             scheduler=dict(mem=4000, queue=lsf_queue)),
        Task(fn=interpro.oracle.proteins.export_uniprot_proteins,
             args=(ipr_pro_uri, df.proteins),
             name="export-proteins",
             scheduler=dict(mem=1000, queue=lsf_queue)),
        Task(fn=interpro.oracle.matches.export_residues,
             args=(ipr_pro_uri, df.protein2residues),
             name="export-residues",
             scheduler=dict(mem=100, queue=lsf_queue)),
        Task(fn=interpro.oracle.proteins.export_uniparc_proteins,
             args=(ipr_pro_uri, df.uniparcproteins),
             name="export-uniparc-proteins",
             scheduler=dict(mem=1000, queue=lsf_queue)),
        Task(fn=interpro.oracle.structures.export_pdbe_matches,
             args=(ips_pro_uri, df.pdbematches),
             name="export-pdbe-matches",
             scheduler=dict(mem=2000, queue=lsf_queue)),
        Task(fn=interpro.oracle.structures.export_rosettafold,
             args=(ipr_pro_uri, df.rosettafold),
             name="export-rosettafold",
             scheduler=dict(mem=2000, queue=lsf_queue)),
        Task(fn=interpro.oracle.taxa.export_taxa,
             args=(ipr_pro_uri, df.taxa),
             name="export-taxa",
             scheduler=dict(mem=8000, queue=lsf_queue)),
        Task(fn=pdbe.export_structures,
             args=(pdbe_uri, df.structures),
             name="export-structures",
             scheduler=dict(mem=10000, queue=lsf_queue)),
        Task(fn=pdbe.export_segments,
             args=(pdbe_uri, df.protein2structures),
             name="export-structure-chains",
             scheduler=dict(mem=1000, queue=lsf_queue)),
        Task(fn=pfam.export_alignments,
             args=(pfam_uri, df.pfam_alignments),
             name="export-pfam-alignments",
             scheduler=dict(mem=4000, queue=lsf_queue)),

        # Exports with dependencies
        Task(fn=interpro.oracle.entries.export_pathways,
             args=(ipr_pro_uri, data_dir),
             name="export-pathways",
             requires=["export-entry2xrefs"],
             scheduler=dict(mem=3000, queue=lsf_queue)),
        Task(fn=alphafold.export,
             args=(config["data"]["alphafold"], df.proteins,
                   df.protein2alphafold),
             kwargs=dict(tempdir=temp_dir),
             name="export-alphafold",
             requires=["export-proteins"],
             scheduler=dict(mem=1000, queue=lsf_queue)),
        Task(fn=interpro.oracle.matches.export_features,
             args=(ipr_pro_uri, df.proteins, df.protein2features),
             kwargs=dict(processes=8, tempdir=temp_dir),
             name="export-features",
             requires=["export-proteins"],
             scheduler=dict(cpu=8, mem=4096, queue=lsf_queue)),
        Task(fn=interpro.oracle.matches.export_uniprot_matches,
             args=(ipr_pro_uri, df.proteins, df.protein2matches),
             kwargs=dict(processes=8, tempdir=temp_dir),
             name="export-matches",
             requires=["export-proteins"],
             scheduler=dict(cpu=8, mem=8000, queue=lsf_queue)),
        Task(fn=interpro.oracle.hmms.export_hmms,
             args=(ipr_pro_uri, df.protein2matches, df.hmms),
             name="export-hmms",
             requires=["export-matches"],
             scheduler=dict(mem=2000, queue=lsf_queue)),
        Task(fn=interpro.oracle.matches.export_uniparc_matches,
             args=(ipr_pro_uri, df.uniparcproteins, df.uniparcmatches),
             kwargs=dict(processes=8, tempdir=temp_dir),
             name="export-uniparc-matches",
             requires=["export-uniparc-proteins"],
             scheduler=dict(cpu=8, mem=16000, queue=lsf_queue)),
        Task(fn=interpro.oracle.proteins.export_uniprot_sequences,
             args=(ipr_pro_uri, df.proteins, df.protein2sequence),
             kwargs=dict(tempdir=temp_dir),
             name="export-sequences",
             requires=["export-proteins"],
             scheduler=dict(mem=2000, queue=lsf_queue)),
        Task(fn=uniprot.proteomes.export_proteomes,
             args=(uniprot_uri, df.proteomes),
             name="export-reference-proteomes",
             scheduler=dict(mem=100, queue=lsf_queue)),
        Task(fn=uniprot.proteins.export_entry2evidence,
             args=(uniprot_uri, df.proteins, df.protein2evidence),
             kwargs=dict(tempdir=temp_dir),
             name="export-evidences",
             requires=["export-proteins"],
             scheduler=dict(mem=1000, queue=lsf_queue)),
        Task(fn=uniprot.proteins.export_entry2functions,
             args=(uniprot_uri, df.proteins, df.protein2functions),
             kwargs=dict(tempdir=temp_dir),
             name="export-functions",
             requires=["export-proteins"],
             scheduler=dict(mem=2000, queue=lsf_queue)),
        Task(fn=uniprot.proteins.export_entry2name,
             args=(uniprot_uri, df.proteins, df.protein2name),
             kwargs=dict(tempdir=temp_dir),
             name="export-names",
             requires=["export-proteins"],
             scheduler=dict(mem=1000, queue=lsf_queue)),
        Task(fn=uniprot.proteins.export_entry2proteome,
             args=(uniprot_uri, df.proteins, df.protein2proteome),
             kwargs=dict(tempdir=temp_dir),
             name="export-proteomes",
             requires=["export-proteins"],
             scheduler=dict(mem=1000, queue=lsf_queue)),
    ]

    tasks += [
        # Add a "group" task, to include all export tasks
        Task(fn=wait,
             name="export",
             requires=get_terminals(tasks)),

        # Match lookup tables
        Task(fn=interpro.oracle.lookup.build_upi_md5_tbl,
             args=(ips_pro_uri,),
             name="lookup-md5",
             scheduler=dict(mem=4000, queue=lsf_queue)),
        Task(fn=interpro.oracle.lookup.build_lookup_tmp_tab,
             args=(ips_pro_uri,),
             name="lookup-insert-matches",
             requires=["lookup-md5"],
             scheduler=dict(mem=4000, queue=lsf_queue)),
        Task(fn=interpro.oracle.lookup.build_lookup_tmp_tab_idx,
             args=(ips_pro_uri,),
             name="lookup-index-matches",
             requires=["lookup-insert-matches"],
             scheduler=dict(mem=4000, queue=lsf_queue)),
        Task(fn=interpro.oracle.lookup.build_site_lookup_tmp_tab,
             args=(ips_pro_uri,),
             name="lookup-insert-sites",
             requires=["lookup-md5"],
             scheduler=dict(mem=4000, queue=lsf_queue)),
        Task(fn=interpro.oracle.lookup.build_site_lookup_tmp_tab_idx,
             args=(ips_pro_uri,),
             name="lookup-index-sites",
             requires=["lookup-insert-sites"],
             scheduler=dict(mem=4000, queue=lsf_queue)),
        Task(fn=wait,
             name="lookup",
             requires=["lookup-index-matches", "lookup-index-sites"]),
    ]

    xrefs_tasks = [
        Task(fn=interpro.xrefs.domorgs.export,
             args=(df.proteins, df.protein2matches, df.protein2domorg),
             kwargs=dict(processes=16, tempdir=temp_dir),
             name="export-dom-orgs",
             requires=["export-matches"],
             scheduler=dict(cpu=16, mem=16000, queue=lsf_queue)),
        Task(fn=interpro.xrefs.entries.export_sim_entries,
             args=(df.protein2matches, df.overlapping),
             name="export-sim-entries",
             requires=["export-matches"],
             scheduler=dict(mem=2000, queue=lsf_queue)),
        Task(fn=interpro.xrefs.entries.export_clan_xrefs,
             args=(df.clans, df.proteins, df.protein2matches,
                   df.protein2proteome, df.protein2domorg,
                   df.protein2structures, df.clan2xrefs),
             kwargs=dict(processes=8, tempdir=temp_dir),
             name="export-clan2xrefs",
             requires=["export-clans", "export-proteomes", "export-dom-orgs",
                       "export-structure-chains"],
             scheduler=dict(cpu=8, mem=10000, queue=lsf_queue)),
        Task(fn=interpro.xrefs.entries.export_xrefs,
             args=(uniprot_uri, df.proteins, df.protein2matches,
                   df.protein2alphafold, df.protein2proteome,
                   df.protein2domorg, df.rosettafold, df.protein2structures,
                   df.protein2evidence, df.taxa, config["data"]["metacyc"],
                   df.entry2xrefs),
             kwargs=dict(interpro_uri=ipr_pro_uri if update_db else None,
                         processes=16, tempdir=temp_dir),
             name="export-entry2xrefs",
             requires=["export-alphafold", "export-proteomes",
                       "export-dom-orgs", "export-structure-chains",
                       "export-taxa", "export-rosettafold",
                       "export-evidences"],
             scheduler=dict(cpu=16, mem=24000, queue=lsf_queue)),
        Task(fn=interpro.xrefs.proteomes.export_xrefs,
             args=(df.clans, df.proteins, df.protein2matches,
                   df.protein2proteome, df.protein2domorg,
                   df.protein2structures, df.proteomes, df.proteome2xrefs),
             kwargs=dict(processes=8, tempdir=temp_dir),
             name="export-proteome2xrefs",
             requires=["export-clans", "export-proteomes", "export-dom-orgs",
                       "export-structure-chains",
                       "export-reference-proteomes"],
             scheduler=dict(cpu=8, mem=10000, queue=lsf_queue)),
        Task(fn=interpro.xrefs.structures.export_xrefs,
             args=(df.clans, df.proteins, df.protein2matches,
                   df.protein2proteome, df.protein2domorg, df.structures,
                   df.protein2structures, df.structure2xrefs),
             name="export-structure2xrefs",
             requires=["export-clans", "export-proteomes", "export-dom-orgs",
                       "export-structures", "export-structure-chains"],
             scheduler=dict(mem=10000, queue=lsf_queue)),
        Task(fn=interpro.xrefs.taxa.export_xrefs,
             args=(df.proteins, df.protein2matches, df.protein2proteome,
                   df.protein2structures, df.taxa, df.taxon2xrefs),
             kwargs=dict(processes=16, tempdir=temp_dir),
             name="export-taxon2xrefs",
             requires=["export-matches", "export-proteomes",
                       "export-structure-chains", "export-taxa"],
             scheduler=dict(cpu=16, mem=24000, queue=lsf_queue)),
    ]

    tasks += xrefs_tasks
    tasks += [
        Task(fn=wait,
             name="xrefs",
             requires=get_terminals(tasks, [t.name for t in xrefs_tasks])),
        Task(fn=interpro.email.notify_curators,
             args=(config["email"]["server"],
                   config["email"]["from"],
                   config["email"]["to"]),
             name="notify-curators",
             requires=["export", "xrefs"])
    ]

    mysql_tasks = [
        Task(fn=interpro.mysql.entries.populate_annotations,
             args=(ipr_stg_uri, df.entries, df.hmms, df.pfam_alignments),
             name="insert-annotations",
             requires=["export-entries", "export-hmms",
                       "export-pfam-alignments"],
             scheduler=dict(mem=5000, queue=lsf_queue)),
        Task(fn=interpro.mysql.entries.index_annotations,
             args=(ipr_stg_uri,),
             name="index-annotations",
             requires=["insert-annotations"],
             scheduler=dict(queue=lsf_queue)),
        Task(fn=interpro.mysql.clans.populate,
             args=(ipr_stg_uri, df.clans, df.clan2xrefs, df.clans_alignments),
             name="insert-clans",
             requires=["export-clan2xrefs"],
             scheduler=dict(mem=2000, queue=lsf_queue)),
        Task(fn=interpro.mysql.databases.populate_databases,
             args=(ipr_stg_uri, df.databases),
             name="insert-databases",
             requires=["export-databases"],
             scheduler=dict(mem=100, queue=lsf_queue)),
        Task(fn=interpro.mysql.entries.populate_entries,
             args=(ipr_stg_uri, pfam_uri, df.clans, df.entries,
                   df.overlapping, df.entry2xrefs),
             name="insert-entries",
             requires=["export-clans", "export-entries",
                       "export-sim-entries", "export-entry2xrefs"],
             scheduler=dict(mem=10000, queue=lsf_queue)),
        Task(fn=interpro.mysql.entries.index_entries,
             args=(ipr_stg_uri,),
             name="index-entries",
             requires=["insert-entries"],
             scheduler=dict(queue=lsf_queue)),
        Task(fn=interpro.mysql.entries.populate_entry_taxa_distrib,
             args=(ipr_stg_uri, df.entries, df.entry2xrefs),
             name="insert-entries-taxa",
             requires=["export-entries", "export-entry2xrefs"],
             scheduler=dict(mem=8000, queue=lsf_queue)),
        Task(fn=interpro.mysql.proteins.populate_features,
             args=(ipr_stg_uri, df.protein2features),
             name="insert-features",
             requires=["export-features"],
             scheduler=dict(mem=1000, queue=lsf_queue)),
        Task(fn=interpro.mysql.proteins.index_features,
             args=(ipr_stg_uri,),
             name="index-features",
             requires=["insert-features"],
             scheduler=dict(queue=lsf_queue)),
        Task(fn=interpro.mysql.proteins.populate_isoforms,
             args=(ipr_stg_uri, df.isoforms),
             name="insert-isoforms",
             requires=["export-isoforms"],
             scheduler=dict(mem=1000, queue=lsf_queue)),
        Task(fn=interpro.mysql.proteins.populate_residues,
             args=(ipr_stg_uri, df.protein2residues),
             name="insert-residues",
             requires=["export-residues"],
             scheduler=dict(mem=1000, queue=lsf_queue)),
        Task(fn=interpro.mysql.proteins.index_residues,
             args=(ipr_stg_uri,),
             name="index-residues",
             requires=["insert-residues"],
             scheduler=dict(queue=lsf_queue)),
        Task(fn=interpro.mysql.proteins.populate_proteins,
             args=(ipr_stg_uri, df.clans, df.entries, df.isoforms,
                   df.structures, df.protein2structures, df.taxa, df.proteins,
                   df.protein2domorg, df.protein2evidence,
                   df.protein2functions, df.protein2matches, df.protein2name,
                   df.protein2proteome, df.protein2sequence),
             name="insert-proteins",
             requires=["export-clans", "export-entries", "export-isoforms",
                       "export-structures", "export-structure-chains",
                       "export-taxa", "export-dom-orgs", "export-evidences",
                       "export-functions", "export-names", "export-proteomes",
                       "export-sequences"],
             scheduler=dict(mem=8000, queue=lsf_queue)),
        Task(fn=interpro.mysql.proteins.index_proteins,
             args=(ipr_stg_uri,),
             name="index-proteins",
             requires=["insert-proteins"],
             scheduler=dict(queue=lsf_queue)),
        Task(fn=interpro.mysql.proteomes.populate,
             args=(ipr_stg_uri, df.proteomes, df.proteome2xrefs),
             name="insert-proteomes",
             requires=["export-proteome2xrefs"],
             scheduler=dict(mem=500, queue=lsf_queue)),
        Task(fn=interpro.mysql.structures.populate_rosettafold,
             args=(ipr_stg_uri, df.rosettafold),
             name="insert-rosettafold",
             requires=["export-rosettafold"],
             scheduler=dict(mem=500, queue=lsf_queue)),
        Task(fn=interpro.mysql.structures.populate_structures,
             args=(ipr_stg_uri, df.structures, df.protein2structures,
                   df.structure2xrefs),
             name="insert-structures",
             requires=["export-structure2xrefs"],
             scheduler=dict(mem=8000, queue=lsf_queue)),
        Task(fn=interpro.mysql.taxa.populate,
             args=(ipr_stg_uri, df.taxa, df.taxon2xrefs),
             name="insert-taxa",
             requires=["export-taxon2xrefs"],
             scheduler=dict(mem=4000, queue=lsf_queue)),
        Task(fn=interpro.mysql.taxa.index,
             args=(ipr_stg_uri,),
             name="index-taxa",
             requires=["insert-taxa"],
             scheduler=dict(queue=lsf_queue)),
        Task(fn=interpro.mysql.databases.populate_rel_notes,
             args=(ipr_stg_uri, ipr_rel_uri, df.clans, df.entries,
                   df.proteomes, df.protein2structures, df.structures,
                   df.taxa, df.proteins, df.protein2matches,
                   df.protein2proteome),
             name="insert-release-notes",
             scheduler=dict(mem=12000, queue=lsf_queue),
             requires=["export-clans", "export-entries", "export-matches",
                       "export-reference-proteomes", "export-structure-chains",
                       "export-structures", "export-taxa", "insert-databases"])
    ]

    tasks += mysql_tasks
    tasks += [
        Task(fn=wait,
             name="mysql",
             requires=get_terminals(tasks, [t.name for t in mysql_tasks])),
    ]

    # Tasks to index documents in Elasticsearch
    es_tasks = [
        Task(fn=interpro.elastic.export_documents,
             args=(df.proteins, df.protein2matches, df.protein2domorg,
                   df.protein2proteome, df.protein2structures,
                   df.protein2alphafold, df.proteomes, df.structures,
                   df.clans, df.entries, df.taxa, es_dirs, release_version),
             name="es-export",
             requires=["export-dom-orgs", "export-proteomes",
                       "export-structure-chains", "export-alphafold",
                       "export-reference-proteomes", "export-structures",
                       "export-clans", "export-entries", "export-taxa"],
             scheduler=dict(mem=20000, queue=lsf_queue))
    ]

    for cluster, hosts, cluster_dir in es_clusters:
        es_tasks += [
            Task(
                fn=interpro.elastic.create_indices,
                args=(df.databases, hosts, release_version),
                name=f"es-init-{cluster}",
                scheduler=dict(queue=lsf_queue),
                requires=["export-databases"] + list(es_tasks[0].requires)
            ),
            Task(
                fn=interpro.elastic.index_documents,
                args=(hosts, cluster_dir, release_version),
                kwargs=dict(threads=8),
                name=f"es-index-{cluster}",
                scheduler=dict(mem=8000, queue=lsf_queue),
                requires=[f"es-init-{cluster}"]
            )
        ]

    tasks += es_tasks
    tasks += [
        Task(fn=wait,
             name="elastic",
             requires=get_terminals(tasks, [t.name for t in es_tasks])),
    ]

    for cluster, hosts, cluster_dir in es_clusters:
        tasks += [
            Task(
                fn=interpro.elastic.publish,
                args=(hosts,),
                name=f"es-publish-{cluster}",
                scheduler=dict(queue=lsf_queue),
                requires=["es-export", f"es-index-{cluster}"]
            )
        ]

    # Task for other EMBL-EBI services
    service_tasks = [
        Task(
            fn=ebisearch.export,
            args=(df.clans, df.databases, df.entries, df.taxa,
                  df.entry2xrefs, os.path.join(data_dir, "ebisearch")),
            name="export-ebisearch",
            scheduler=dict(mem=20000, queue=lsf_queue),
            requires=["export-clans", "export-databases", "export-entries",
                      "export-entry2xrefs"]
        ),
        Task(
            fn=uniprot.goa.export,
            args=(df.databases, df.entries, df.structures, df.pdbematches,
                  df.entry2xrefs, os.path.join(data_dir, "goa")),
            name="export-goa",
            scheduler=dict(mem=8000, queue=lsf_queue),
            requires=["export-databases", "export-entries",
                      "export-structures", "export-pdbe-matches",
                      "export-entry2xrefs"]
        ),
        # Task(
        #     fn=pdbe.export_pdb_matches,
        #     args=(df.databases, df.pdbematches,
        #           os.path.join(data_dir, "pdbe")),
        #     name="export-pdbe",
        #     scheduler=dict(queue=lsf_queue),
        #     requires=["export-databases", "export-pdbe-matches"]
        # )
    ]

    # Add tasks for FTP
    exchange_tasks = service_tasks + [
        Task(fn=interpro.ftp.xmlfiles.export_feature_matches,
             args=(df.databases, df.proteins, df.protein2features, pub_dir),
             name="ftp-features",
             requires=["export-databases", "export-proteins",
                       "export-features"],
             scheduler=dict(mem=1000, queue=lsf_queue)),
        Task(fn=interpro.ftp.flatfiles.export,
             args=(df.entries, df.protein2matches, pub_dir),
             name="ftp-flatfiles",
             requires=["export-entries", "export-matches"],
             scheduler=dict(mem=8000, queue=lsf_queue)),
        Task(fn=interpro.ftp.xmlfiles.export_interpro,
             args=(df.entries, df.entry2xrefs, df.databases, df.taxa, pub_dir),
             name="ftp-interpro",
             requires=["export-entries", "export-entry2xrefs",
                       "export-databases"],
             scheduler=dict(mem=10000, queue=lsf_queue)),
        Task(fn=interpro.ftp.xmlfiles.export_matches,
             args=(df.databases, df.isoforms, df.proteins,
                   df.protein2matches, pub_dir),
             kwargs=dict(processes=8),
             name="ftp-matches",
             requires=["export-databases", "export-isoforms",
                       "export-matches"],
             scheduler=dict(cpu=8, mem=16000, queue=lsf_queue)),
        Task(fn=interpro.ftp.relnotes.export,
             args=(ipr_stg_uri, pub_dir),
             name="ftp-relnotes",
             requires=["insert-release-notes"],
             scheduler=dict(queue=lsf_queue)),
        Task(
            fn=interpro.ftp.xmlfiles.export_structure_matches,
            args=(df.structures, df.proteins, df.protein2structures, pub_dir),
            name="ftp-structures",
            scheduler=dict(mem=8000, queue=lsf_queue),
            requires=["export-structures", "export-proteins",
                      "export-structure-chains"]
        ),
        Task(fn=interpro.ftp.uniparc.archive_matches,
             args=(df.uniparcproteins, df.uniparcmatches, pub_dir),
             kwargs=dict(processes=8),
             name="ftp-uniparc",
             requires=["export-uniparc-matches"],
             scheduler=dict(cpu=8, mem=8000, queue=lsf_queue)),
    ]

    tasks += exchange_tasks
    tasks += [
        Task(fn=wait,
             name="exchange",
             requires=get_terminals(tasks, [t.name for t in exchange_tasks])),
    ]

    tasks += [
        Task(
            fn=ebisearch.publish,
            args=(os.path.join(data_dir, "ebisearch"),
                  config["exchange"]["ebisearch"]),
            name="publish-ebisearch",
            scheduler=dict(queue=lsf_queue),
            requires=["export-ebisearch"]
        ),
        Task(
            fn=uniprot.goa.publish,
            args=(os.path.join(data_dir, "goa"), config["exchange"]["goa"]),
            name="publish-goa",
            scheduler=dict(queue=lsf_queue),
            requires=["export-goa"]
        ),
        # Task(
        #     fn=pdbe.publish,
        #     args=(os.path.join(data_dir, "pdbe"), config["exchange"]["pdbe"]),
        #     name="publish-pdbe",
        #     scheduler=dict(queue=lsf_queue),
        #     requires=["export-pdbe"]
        # )
    ]

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
