import argparse
import importlib.metadata
import os
import tomllib
import time

from mundone import Task, Workflow, get_terminals

from interpro7dw import alphafold, ebisearch, interpro, pdbe, pfam, uniprot


def wait(secs: int = 5):
    time.sleep(secs)


class DataFiles:
    def __init__(self, root: str):
        # BasicStores
        self.clan2xrefs = os.path.join(root, "clan-xrefs")
        self.entry2xrefs = os.path.join(root, "entry-xrefs")
        self.hmms = os.path.join(root, "hmms")
        self.isoforms = os.path.join(root, "isoforms")
        self.pdbematches = os.path.join(root, "pdbe-matches")
        self.pfam_alignments = os.path.join(root, "pfam-alignments")
        self.protein2residues = os.path.join(root, "protein-residues")
        self.proteome2xrefs = os.path.join(root, "proteome-xrefs")
        self.structure2xrefs = os.path.join(root, "structure-xrefs")
        self.taxon2xrefs = os.path.join(root, "taxon-xrefs")

        # KVStores
        self.proteins = os.path.join(root, "proteins")
        self.protein2alphafold = os.path.join(root, "protein-alphafold")
        self.protein2domorg = os.path.join(root, "protein-domorg")
        self.protein2evidence = os.path.join(root, "protein-evidence")
        self.protein2features = os.path.join(root, "protein-features")
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
        self.uniprot2pdb = os.path.join(root, "protein-structures")
        self.structures = os.path.join(root, "structures")
        self.cath_scop = os.path.join(root, "cath-scop")
        self.taxa = os.path.join(root, "taxa")


def gen_tasks(config: dict) -> list[Task]:
    release_version = config["release"]["version"]
    release_date = config["release"]["date"]
    update_db = config["release"]["update"]
    data_dir = config["data"]["path"]
    data_src_dir = config["data"]["src"]
    temp_dir = config["data"]["tmp"]
    ipr_pro_uri = config["databases"]["interpro"]["production"]
    ipr_stg_uri = config["databases"]["interpro"]["staging"]
    ipr_rel_uri = config["databases"]["interpro"]["release"]
    ips_pro_uri = config["databases"]["iprscan"]["production"]
    goa_uri = config["databases"]["goa"]
    intact_uri = config["databases"]["intact"]
    pdbe_uri = config["databases"]["pdbe"]
    uniprot_uri = config["databases"]["uniprot"]
    pub_dir = os.path.join(config["exchange"]["interpro"], release_version)
    scheduler, queue = parse_scheduler(config["workflow"]["scheduler"])

    es_clusters = []
    es_dirs = [os.path.join(data_dir, "elastic", "default")]
    for cluster, properties in config["elasticsearch"].items():
        path = os.path.join(data_dir, "elastic", cluster)
        es_clusters.append({
            "id": cluster,
            "hosts": list(set(properties["nodes"])),
            "user": properties["user"],
            "password": properties["password"],
            "fingerprint": properties["fingerprint"],
            "path": path
        })
        es_dirs.append(path)

    df = DataFiles(data_dir)

    tasks = [
        # Exports without dependencies
        Task(fn=interpro.oracle.clans.export_clans,
             args=(ipr_pro_uri, df.clans),
             name="export-clans",
             scheduler=dict(type=scheduler, queue=queue, mem=2000, hours=1)),
        Task(fn=interpro.oracle.databases.export,
             args=(ipr_pro_uri, goa_uri, release_version, release_date,
                   df.databases),
             kwargs=dict(update=update_db),
             name="export-databases",
             scheduler=dict(type=scheduler, queue=queue, mem=1000, hours=1)),
        Task(fn=interpro.oracle.entries.export_entries,
             args=(ipr_pro_uri, goa_uri, intact_uri, df.entries),
             name="export-entries",
             scheduler=dict(type=scheduler, queue=queue, mem=3000, hours=1)),
        Task(fn=interpro.oracle.matches.export_isoforms,
             args=(ipr_pro_uri, df.isoforms),
             name="export-isoforms",
             scheduler=dict(type=scheduler, queue=queue, mem=4000, hours=1)),
        Task(fn=interpro.oracle.proteins.export_uniprot_proteins,
             args=(ipr_pro_uri, df.proteins),
             name="export-proteins",
             scheduler=dict(type=scheduler, queue=queue, mem=1000, hours=2)),
        Task(fn=interpro.oracle.matches.export_residues,
             args=(ipr_pro_uri, df.protein2residues),
             name="export-residues",
             scheduler=dict(type=scheduler, queue=queue, mem=800, hours=10)),
        Task(fn=interpro.oracle.structures.export_matches,
             args=(ipr_pro_uri, pdbe_uri, df.pdbematches),
             name="export-pdb-matches",
             scheduler=dict(type=scheduler, queue=queue, mem=3000, hours=36)),
        Task(fn=interpro.oracle.proteins.export_uniparc_proteins,
             args=(ipr_pro_uri, df.uniparcproteins),
             name="export-uniparc-proteins",
             scheduler=dict(type=scheduler, queue=queue, mem=1000, hours=3)),
        Task(fn=interpro.oracle.taxa.export_taxa,
             args=(ipr_pro_uri, df.taxa),
             name="export-taxa",
             scheduler=dict(type=scheduler, queue=queue, mem=8000, hours=1)),
        Task(fn=pdbe.export_entries,
             args=(pdbe_uri, df.structures),
             name="export-structures",
             scheduler=dict(type=scheduler, queue=queue, mem=10000, hours=1)),
        Task(fn=pdbe.export_cath_scop,
             args=(pdbe_uri, df.cath_scop),
             name="export-cath-scop",
             scheduler=dict(type=scheduler, queue=queue, mem=2000, hours=1)),
        Task(fn=pdbe.export_uniprot2pdb,
             args=(pdbe_uri, df.uniprot2pdb),
             name="export-uniprot2pdb",
             scheduler=dict(type=scheduler, queue=queue, mem=1000, hours=1)),
        Task(fn=pfam.export_alignments,
             args=(ipr_pro_uri, df.pfam_alignments),
             name="export-pfam-alignments",
             scheduler=dict(type=scheduler, queue=queue, mem=4000, hours=2)),

        # Exports with dependencies
        Task(fn=alphafold.export,
             args=(config["data"]["alphafold"], df.proteins,
                   df.protein2alphafold),
             kwargs=dict(tempdir=temp_dir),
             name="export-alphafold",
             requires=["export-proteins"],
             scheduler=dict(type=scheduler, queue=queue, mem=1000, hours=1)),
        Task(fn=interpro.oracle.matches.export_features,
             args=(ipr_pro_uri, df.proteins, df.protein2features),
             kwargs=dict(processes=8, tempdir=temp_dir),
             name="export-features",
             requires=["export-proteins"],
             scheduler=dict(type=scheduler, queue=queue, cpu=8, mem=4096,
                            hours=9)),
        Task(fn=interpro.oracle.matches.export_uniprot_matches,
             args=(ipr_pro_uri, df.proteins, df.protein2matches),
             kwargs=dict(processes=8, tempdir=temp_dir),
             name="export-matches",
             requires=["export-proteins"],
             scheduler=dict(type=scheduler, queue=queue, cpu=8, mem=10000,
                            hours=36)),
        Task(fn=interpro.oracle.hmms.export_hmms,
             args=(ipr_pro_uri, df.protein2matches, df.hmms),
             name="export-hmms",
             requires=["export-matches"],
             scheduler=dict(type=scheduler, queue=queue, mem=2000, hours=15)),
        Task(fn=interpro.oracle.matches.export_uniparc_matches,
             args=(ipr_pro_uri, df.uniparcproteins, df.uniparcmatches),
             kwargs=dict(processes=8, tempdir=temp_dir),
             name="export-uniparc-matches",
             requires=["export-uniparc-proteins"],
             scheduler=dict(type=scheduler, queue=queue, cpu=8, mem=16000,
                            hours=30)),
        Task(fn=interpro.oracle.proteins.export_uniprot_sequences,
             args=(ipr_pro_uri, df.proteins, df.protein2sequence),
             kwargs=dict(tempdir=temp_dir),
             name="export-sequences",
             requires=["export-proteins"],
             scheduler=dict(type=scheduler, queue=queue, mem=2000, hours=10)),
        Task(fn=uniprot.proteomes.export_proteomes,
             args=(uniprot_uri, df.proteomes),
             name="export-reference-proteomes",
             scheduler=dict(type=scheduler, queue=queue, mem=100, hours=1)),
        Task(fn=uniprot.proteins.export_entry2evidence,
             args=(uniprot_uri, df.proteins, df.protein2evidence),
             kwargs=dict(tempdir=temp_dir),
             name="export-evidences",
             requires=["export-proteins"],
             scheduler=dict(type=scheduler, queue=queue, mem=1000, hours=3)),
        Task(fn=uniprot.proteins.export_entry2functions,
             args=(uniprot_uri, df.proteins, df.protein2functions),
             kwargs=dict(tempdir=temp_dir),
             name="export-functions",
             requires=["export-proteins"],
             scheduler=dict(type=scheduler, queue=queue, mem=2000, hours=2)),
        Task(fn=uniprot.proteins.export_entry2name,
             args=(uniprot_uri, df.proteins, df.protein2name),
             kwargs=dict(tempdir=temp_dir),
             name="export-names",
             requires=["export-proteins"],
             scheduler=dict(type=scheduler, queue=queue, mem=1000, hours=3)),
        Task(fn=uniprot.proteins.export_entry2proteome,
             args=(uniprot_uri, df.proteins, df.protein2proteome),
             kwargs=dict(tempdir=temp_dir),
             name="export-proteomes",
             requires=["export-proteins"],
             scheduler=dict(type=scheduler, queue=queue, mem=1000, hours=1)),
    ]

    tasks += [
        # Add a "group" task, to include all export tasks
        Task(fn=wait,
             name="export",
             requires=get_terminals(tasks)),
    ]

    xrefs_tasks = [
        Task(fn=interpro.xrefs.domorgs.export,
             args=(df.proteins, df.protein2matches, df.protein2domorg),
             kwargs=dict(processes=16, tempdir=temp_dir),
             name="export-dom-orgs",
             requires=["export-matches"],
             scheduler=dict(type=scheduler, queue=queue, cpu=16, mem=16000,
                            hours=3)),
        Task(fn=interpro.xrefs.entries.export_sim_entries,
             args=(df.protein2matches, df.overlapping),
             name="export-sim-entries",
             requires=["export-matches"],
             scheduler=dict(type=scheduler, queue=queue, mem=2000, hours=14)),
        Task(fn=interpro.xrefs.clans.export_xrefs,
             args=(df.clans, df.proteins, df.protein2matches,
                   df.protein2proteome, df.protein2domorg,
                   df.pdbematches, df.clan2xrefs),
             kwargs=dict(processes=16, tempdir=temp_dir),
             name="export-clan2xrefs",
             requires=["export-clans", "export-proteomes", "export-dom-orgs",
                       "export-pdb-matches"],
             scheduler=dict(type=scheduler, queue=queue, cpu=16, mem=24000,
                            hours=3)),
        Task(fn=interpro.xrefs.entries.export_xrefs,
             args=(uniprot_uri, df.proteins, df.protein2matches,
                   df.protein2alphafold, df.protein2proteome,
                   df.protein2domorg, df.pdbematches,
                   df.protein2evidence, df.taxa, config["data"]["metacyc"],
                   df.entry2xrefs),
             kwargs=dict(interpro_uri=ipr_pro_uri if update_db else None,
                         processes=16, tempdir=temp_dir),
             name="export-entry2xrefs",
             requires=["export-alphafold", "export-proteomes",
                       "export-dom-orgs", "export-pdb-matches",
                       "export-taxa", "export-evidences"],
             scheduler=dict(type=scheduler, queue=queue, cpu=16, mem=30000,
                            hours=24)),
        Task(fn=interpro.xrefs.proteomes.export_xrefs,
             args=(df.proteins, df.protein2matches, df.protein2proteome,
                   df.structures, df.uniprot2pdb, df.pdbematches, df.proteomes,
                   df.proteome2xrefs),
             kwargs=dict(processes=16, tempdir=temp_dir),
             name="export-proteome2xrefs",
             requires=["export-matches", "export-proteomes",
                       "export-structures", "export-uniprot2pdb",
                       "export-pdb-matches", "export-reference-proteomes"],
             scheduler=dict(type=scheduler, queue=queue, cpu=16, mem=48000,
                            hours=6)),
        Task(fn=interpro.xrefs.structures.export_xrefs,
             args=(df.clans, df.proteins, df.protein2proteome,
                   df.protein2domorg, df.structures, df.pdbematches,
                   df.uniprot2pdb, df.structure2xrefs),
             name="export-structure2xrefs",
             requires=["export-clans", "export-proteomes", "export-dom-orgs",
                       "export-structures", "export-pdb-matches",
                       "export-uniprot2pdb"],
             scheduler=dict(type=scheduler, queue=queue, mem=10000, hours=1)),
        Task(fn=interpro.xrefs.taxa.export_xrefs,
             args=(df.proteins, df.protein2matches, df.protein2proteome,
                   df.structures, df.uniprot2pdb, df.pdbematches, df.taxa,
                   df.taxon2xrefs),
             kwargs=dict(processes=16, tempdir=temp_dir),
             name="export-taxon2xrefs",
             requires=["export-matches", "export-proteomes",
                       "export-structures", "export-uniprot2pdb",
                       "export-pdb-matches", "export-taxa"],
             scheduler=dict(type=scheduler, queue=queue, cpu=16, mem=48000,
                            hours=18)),
    ]
    tasks += xrefs_tasks
    tasks += [
        Task(fn=wait,
             name="xrefs",
             requires=get_terminals(tasks, [t.name for t in xrefs_tasks])),
    ]

    # InterProScan tasks
    tasks += [
        # Data files for InterProScan
        Task(fn=interpro.ftp.iprscan.package_data,
             args=(ipr_pro_uri, goa_uri, data_src_dir, release_version,
                   data_dir),
             name="interproscan",
             requires=["export-entry2xrefs"],
             scheduler=dict(type=scheduler, queue=queue, mem=4000, hours=6)),
        # Match lookup tables
        Task(fn=interpro.oracle.lookup.build_upi_md5_table,
             args=(ips_pro_uri,),
             name="lookup-md5",
             scheduler=dict(type=scheduler, queue=queue, mem=4000, hours=2)),
        Task(fn=interpro.oracle.lookup.build_matches_table,
             args=(ips_pro_uri,),
             name="lookup-matches",
             requires=["lookup-md5"],
             scheduler=dict(type=scheduler, queue=queue, mem=4000, hours=80)),
        Task(fn=interpro.oracle.lookup.build_site_table,
             args=(ips_pro_uri,),
             name="lookup-sites",
             requires=["lookup-matches"],
             scheduler=dict(type=scheduler, queue=queue, mem=4000, hours=96)),
        Task(fn=wait,
             name="lookup",
             requires=["lookup"]),
    ]

    mysql_tasks = [
        Task(fn=interpro.mysql.entries.populate_annotations,
             args=(ipr_stg_uri, df.entries, df.hmms, df.pfam_alignments),
             name="insert-annotations",
             requires=["export-entries", "export-hmms",
                       "export-pfam-alignments"],
             scheduler=dict(type=scheduler, queue=queue, mem=5000, hours=6)),
        Task(fn=interpro.mysql.entries.index_annotations,
             args=(ipr_stg_uri,),
             name="index-annotations",
             requires=["insert-annotations"],
             scheduler=dict(type=scheduler, queue=queue, mem=100, hours=3)),
        Task(fn=interpro.mysql.clans.populate,
             args=(ipr_stg_uri, df.clans, df.clan2xrefs),
             name="insert-clans",
             requires=["export-clan2xrefs"],
             scheduler=dict(type=scheduler, queue=queue, mem=2000, hours=6)),
        Task(fn=interpro.mysql.databases.populate_databases,
             args=(ipr_stg_uri, df.databases),
             name="insert-databases",
             requires=["export-databases"],
             scheduler=dict(type=scheduler, queue=queue, mem=100, hours=1)),
        Task(fn=interpro.mysql.entries.populate_entries,
             args=(ipr_pro_uri, ipr_stg_uri, df.clans, df.entries,
                   df.overlapping, df.entry2xrefs, df.structures),
             name="insert-entries",
             requires=["export-clans", "export-entries",
                       "export-sim-entries", "export-entry2xrefs",
                       "export-structures"],
             scheduler=dict(type=scheduler, queue=queue, mem=10000, hours=3)),
        Task(fn=interpro.mysql.entries.index_entries,
             args=(ipr_stg_uri,),
             name="index-entries",
             requires=["insert-entries"],
             scheduler=dict(type=scheduler, queue=queue, mem=100, hours=1)),
        Task(fn=interpro.mysql.entries.populate_entry_taxa_distrib,
             args=(ipr_stg_uri, df.entries, df.entry2xrefs),
             name="insert-entries-taxa",
             requires=["export-entries", "export-entry2xrefs"],
             scheduler=dict(type=scheduler, queue=queue, mem=8000, hours=8)),
        Task(fn=interpro.mysql.proteins.populate_features,
             args=(ipr_stg_uri, df.protein2features),
             name="insert-features",
             requires=["export-features"],
             scheduler=dict(type=scheduler, queue=queue, mem=1000, hours=14)),
        Task(fn=interpro.mysql.proteins.index_features,
             args=(ipr_stg_uri,),
             name="index-features",
             requires=["insert-features"],
             scheduler=dict(type=scheduler, queue=queue, mem=100, hours=10)),
        Task(fn=interpro.mysql.proteins.populate_isoforms,
             args=(ipr_stg_uri, df.isoforms),
             name="insert-isoforms",
             requires=["export-isoforms"],
             scheduler=dict(type=scheduler, queue=queue, mem=1000, hours=1)),
        Task(fn=interpro.mysql.proteins.populate_residues,
             args=(ipr_stg_uri, df.protein2residues),
             name="insert-residues",
             requires=["export-residues"],
             scheduler=dict(type=scheduler, queue=queue, mem=1000, hours=8)),
        Task(fn=interpro.mysql.proteins.index_residues,
             args=(ipr_stg_uri,),
             name="index-residues",
             requires=["insert-residues"],
             scheduler=dict(type=scheduler, queue=queue, mem=100, hours=2)),
        Task(fn=interpro.mysql.proteins.populate_proteins,
             args=(ipr_stg_uri, df.clans, df.entries, df.isoforms,
                   df.cath_scop, df.uniprot2pdb, df.taxa, df.proteins,
                   df.protein2domorg, df.protein2evidence,
                   df.protein2functions, df.protein2matches, df.protein2name,
                   df.protein2proteome, df.protein2sequence,
                   df.protein2alphafold),
             name="insert-proteins",
             requires=["export-clans", "export-entries", "export-isoforms",
                       "export-cath-scop", "export-uniprot2pdb",
                       "export-taxa", "export-dom-orgs", "export-evidences",
                       "export-functions", "export-names", "export-proteomes",
                       "export-sequences", "export-alphafold"],
             scheduler=dict(type=scheduler, queue=queue, mem=8000, hours=35)),
        Task(fn=interpro.mysql.proteins.index_proteins,
             args=(ipr_stg_uri,),
             name="index-proteins",
             requires=["insert-proteins"],
             scheduler=dict(type=scheduler, queue=queue, mem=100, hours=36)),
        Task(fn=interpro.mysql.proteomes.populate,
             args=(ipr_stg_uri, df.proteomes, df.proteome2xrefs),
             name="insert-proteomes",
             requires=["export-proteome2xrefs"],
             scheduler=dict(type=scheduler, queue=queue, mem=500, hours=2)),
        Task(fn=interpro.mysql.proteomes.index,
             args=(ipr_stg_uri,),
             name="index-proteomes",
             requires=["insert-proteomes"],
             scheduler=dict(type=scheduler, queue=queue, mem=100, hours=4)),
        Task(fn=interpro.mysql.structures.populate_structures,
             args=(ipr_stg_uri, df.structures, df.uniprot2pdb, df.pdbematches,
                   df.structure2xrefs),
             name="insert-structures",
             requires=["export-structure2xrefs"],
             scheduler=dict(type=scheduler, queue=queue, mem=8000, hours=3)),
        Task(fn=interpro.mysql.taxa.populate,
             args=(ipr_stg_uri, df.taxa, df.taxon2xrefs),
             name="insert-taxa",
             requires=["export-taxon2xrefs"],
             scheduler=dict(type=scheduler, queue=queue, mem=4000, hours=8)),
        Task(fn=interpro.mysql.taxa.index,
             args=(ipr_stg_uri,),
             name="index-taxa",
             requires=["insert-taxa"],
             scheduler=dict(type=scheduler, queue=queue, mem=100, hours=12)),
        Task(fn=interpro.mysql.databases.populate_rel_notes,
             args=(ipr_stg_uri, ipr_rel_uri, df.clans, df.entries,
                   df.proteomes, df.structures, df.structure2xrefs,
                   df.taxa, df.proteins, df.protein2matches,
                   df.protein2proteome),
             name="insert-release-notes",
             scheduler=dict(type=scheduler, queue=queue, mem=12000, hours=16),
             requires=["export-clans", "export-entries", "export-matches",
                       "export-reference-proteomes", "export-structure2xrefs",
                       "export-taxa", "insert-databases"])
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
                   df.protein2name, df.protein2proteome, df.uniprot2pdb,
                   df.pdbematches, df.protein2alphafold, df.proteomes,
                   df.structures, df.clans, df.entries, df.taxa, es_dirs,
                   release_version),
             kwargs=dict(processes=8, tempdir=temp_dir),
             name="es-export",
             requires=["export-dom-orgs", "export-proteomes", "export-names",
                       "export-uniprot2pdb", "export-pdb-matches",
                       "export-alphafold", "export-reference-proteomes",
                       "export-structures", "export-clans", "export-entries",
                       "export-taxa"],
             scheduler=dict(type=scheduler, queue=queue, cpu=8, mem=48000,
                            hours=24)),
    ]

    for cluster in es_clusters:
        es_tasks += [
            Task(
                fn=interpro.elastic.create_indices,
                args=(df.databases, cluster["hosts"], cluster["user"],
                      cluster["password"], cluster["fingerprint"],
                      release_version),
                name=f"es-init-{cluster['id']}",
                scheduler=dict(type=scheduler, queue=queue, mem=100, hours=1),
                requires=["export-databases"] + list(es_tasks[0].requires)
            ),
            Task(
                fn=interpro.elastic.index_documents,
                args=(cluster["hosts"], cluster["user"], cluster["password"],
                      cluster["fingerprint"], cluster["path"], release_version),
                kwargs=dict(processes=8),
                name=f"es-index-{cluster['id']}",
                scheduler=dict(type=scheduler, queue=queue, cpu=8, mem=16000,
                               hours=36),
                requires=[f"es-init-{cluster['id']}"]
            )
        ]

    tasks += es_tasks
    tasks += [
        Task(fn=wait,
             name="elastic",
             requires=get_terminals(tasks, [t.name for t in es_tasks])),
    ]

    for cluster in es_clusters:
        tasks += [
            Task(
                fn=interpro.elastic.publish,
                args=(cluster["hosts"], cluster["user"], cluster["password"],
                      cluster["fingerprint"]),
                name=f"es-publish-{cluster['id']}",
                scheduler=dict(type=scheduler, queue=queue, mem=100, hours=1),
                requires=["es-export", f"es-index-{cluster['id']}"]
            )
        ]

    # Task for other EMBL-EBI services
    ebi_files_tasks = [
        Task(
            fn=ebisearch.export,
            args=(df.clans, df.databases, df.entries, df.taxa,
                  df.entry2xrefs, os.path.join(data_dir, "ebisearch")),
            name="export-ebisearch",
            scheduler=dict(type=scheduler, queue=queue, mem=20000, hours=25),
            requires=["export-clans", "export-databases", "export-entries",
                      "export-entry2xrefs"]
        ),
        Task(
            fn=uniprot.goa.export,
            args=(ipr_pro_uri, df.databases, df.entries, df.structures,
                  df.pdbematches, df.uniprot2pdb, df.entry2xrefs,
                  os.path.join(data_dir, "goa")),
            name="export-goa",
            scheduler=dict(type=scheduler, queue=queue, mem=8000, hours=16),
            requires=["export-databases", "export-entries",
                      "export-structures", "export-pdb-matches",
                      "export-uniprot2pdb", "export-entry2xrefs"]
        ),
    ]

    # Add tasks for FTP
    ftp_files_tasks = [
        Task(fn=interpro.ftp.xmlfiles.export_feature_matches,
             args=(df.databases, df.proteins, df.protein2features, pub_dir),
             kwargs=dict(processes=8),
             name="ftp-features",
             requires=["export-databases", "export-proteins",
                       "export-features"],
             scheduler=dict(type=scheduler, queue=queue, cpu=8, mem=4000,
                            hours=18)),
        Task(fn=interpro.ftp.flatfiles.export,
             args=(df.entries, df.protein2matches, pub_dir),
             name="ftp-flatfiles",
             requires=["export-entries", "export-matches"],
             scheduler=dict(type=scheduler, queue=queue, mem=8000, hours=18)),
        Task(fn=interpro.ftp.xmlfiles.export_interpro,
             args=(df.entries, df.entry2xrefs, df.databases, df.taxa, pub_dir),
             name="ftp-interpro",
             requires=["export-entries", "export-entry2xrefs",
                       "export-databases"],
             scheduler=dict(type=scheduler, queue=queue, mem=10000, hours=3)),
        Task(fn=interpro.ftp.xmlfiles.export_matches,
             args=(df.databases, df.isoforms, df.proteins,
                   df.protein2matches, pub_dir),
             kwargs=dict(processes=8),
             name="ftp-matches",
             requires=["export-databases", "export-isoforms",
                       "export-matches"],
             scheduler=dict(type=scheduler, queue=queue, cpu=8, mem=24000,
                            hours=20)),
        Task(fn=interpro.ftp.relnotes.export,
             args=(ipr_stg_uri, pub_dir),
             name="ftp-relnotes",
             requires=["insert-release-notes"],
             scheduler=dict(type=scheduler, queue=queue, mem=100, hours=1)),
        # Task(
        #     fn=interpro.ftp.xmlfiles.export_structure_matches,
        #     args=(df.structures, df.proteins, df.uniprot2pdb, pub_dir),
        #     name="ftp-structures",
        #     scheduler=dict(type=scheduler, queue=queue,mem=8000, hours=),
        #     requires=["export-structures", "export-proteins",
        #               "export-structure-chains"]
        # ),
        Task(fn=interpro.ftp.uniparc.archive_matches,
             args=(df.uniparcproteins, df.uniparcmatches, pub_dir),
             kwargs=dict(processes=8),
             name="ftp-uniparc",
             requires=["export-uniparc-matches"],
             scheduler=dict(type=scheduler, queue=queue, cpu=8, mem=90000,
                            hours=33)),
    ]

    tasks += ebi_files_tasks + ftp_files_tasks
    tasks += [
        Task(fn=wait,
             name="ebi-files",
             requires=get_terminals(tasks, [t.name for t in ebi_files_tasks])),
        Task(fn=wait,
             name="ftp-files",
             requires=get_terminals(tasks, [t.name for t in ftp_files_tasks])),
    ]

    tasks += [
        Task(
            fn=ebisearch.publish,
            args=(os.path.join(data_dir, "ebisearch"),
                  config["exchange"]["ebisearch"]),
            name="publish-ebisearch",
            scheduler=dict(type=scheduler, queue=queue, mem=500, hours=6),
            requires=["export-ebisearch"]
        ),
        Task(
            fn=uniprot.goa.publish,
            args=(os.path.join(data_dir, "goa"), config["exchange"]["goa"]),
            name="publish-goa",
            scheduler=dict(type=scheduler, queue=queue, mem=500, hours=1),
            requires=["export-goa"]
        ),
    ]

    return tasks


def build():
    parser = argparse.ArgumentParser(
        description="Build InterPro7 data warehouse")
    parser.add_argument("config",
                        metavar="config.toml",
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
    try:
        pkg_version = importlib.metadata.version("interpro7-dw")
    except importlib.metadata.PackageNotFoundError:
        pass
    else:
        parser.add_argument("-v", "--version", action="version",
                            version=f"%(prog)s {pkg_version}",
                            help="show the version and exit")
    args = parser.parse_args()

    if not os.path.isfile(args.config):
        parser.error(f"cannot open '{args.config}': no such file or directory")

    with open(args.config, "rb") as fh:
        config = tomllib.load(fh)

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


def parse_scheduler(value: str) -> tuple[str, str | None]:
    values = value.split(":")
    if len(values) == 2:
        return values[0], values[1]
    elif len(values) == 1:
        return values[0], None
    raise ValueError(value)
