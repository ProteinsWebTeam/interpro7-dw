#!/usr/bin/env python

import configparser
import os
from typing import List

from mundone import Task, Workflow

from interpro7dw import interpro, pdbe, uniprot


class DataFiles:
    def __init__(self, root: str):
        os.makedirs(root, exist_ok=True)

        # Stores
        self.alignments = os.path.join(root, "alignments.store")
        self.proteins = os.path.join(root, "proteins.store")
        self.protein2domorg = os.path.join(root, "protein2domorg.store")
        self.protein2evidence = os.path.join(root, "protein2evidence.store")
        self.protein2features = os.path.join(root, "protein2features.store")
        self.protein2function = os.path.join(root, "protein2function.store")
        self.protein2matches = os.path.join(root, "protein2matches.store")
        self.protein2name = os.path.join(root, "protein2name.store")
        self.protein2proteome = os.path.join(root, "protein2proteome.store")
        self.protein2residues = os.path.join(root, "protein2residues.store")
        self.protein2sequence = os.path.join(root, "protein2sequence.store")

        # Data dumps
        self.clans = os.path.join(root, "clans.pickle")
        self.overlapping_entries = os.path.join(root, "overlapping.pickle")
        self.proteomes = os.path.join(root, "proteomes.pickle")
        self.structures = os.path.join(root, "structures.pickle")
        self.taxa = os.path.join(root, "taxa.pickle")


def gen_tasks(config: configparser.ConfigParser) -> List[Task]:
    ipr_pro_url = config["databases"]["interpro_production"]
    pdbe_url = config["databases"]["pdbe"]
    pfam_url = config["databases"]["pfam"]
    uniprot_url = config["databases"]["uniprot"]
    data_dir = config["data"]["path"]
    temp_dir = config["data"]["tmp"]
    lsf_queue = config["workflow"]["lsf_queue"]

    df = DataFiles(data_dir)

    tasks = [
        # Data from InterPro
        Task(fn=interpro.oracle.taxa.export_taxa,
             args=(ipr_pro_url, df.taxa),
             name="export-taxa",
             scheduler=dict(mem=8000, queue=lsf_queue)),
        Task(fn=interpro.oracle.proteins.export_proteins,
             args=(ipr_pro_url, df.proteins),
             kwargs=dict(tempdir=temp_dir, workers=8),
             name="export-proteins",
             scheduler=dict(cpu=8, mem=16000, scratch=50000, queue=lsf_queue)),
        Task(fn=interpro.oracle.proteins.export_features,
             args=(ipr_pro_url, df.proteins, df.protein2features),
             kwargs=dict(tempdir=temp_dir, workers=8),
             name="export-features",
             requires=["export-proteins"],
             scheduler=dict(cpu=8, mem=16000, scratch=50000, queue=lsf_queue)),
        Task(fn=interpro.oracle.proteins.export_matches,
             args=(ipr_pro_url, df.proteins, df.protein2matches),
             kwargs=dict(tempdir=temp_dir, workers=8),
             name="export-matches",
             requires=["export-proteins"],
             scheduler=dict(cpu=8, mem=16000, scratch=50000, queue=lsf_queue)),
        Task(fn=interpro.oracle.proteins.export_residues,
             args=(ipr_pro_url, df.proteins, df.protein2residues),
             kwargs=dict(tempdir=temp_dir, workers=8),
             name="export-residues",
             requires=["export-proteins"],
             scheduler=dict(cpu=8, mem=16000, scratch=50000, queue=lsf_queue)),
        Task(fn=interpro.oracle.proteins.export_sequences,
             args=(ipr_pro_url, df.proteins, df.protein2sequence),
             kwargs=dict(tempdir=temp_dir, workers=8),
             name="export-sequences",
             requires=["export-proteins"],
             scheduler=dict(cpu=8, mem=16000, scratch=50000, queue=lsf_queue)),
        Task(fn=interpro.oracle.entries.dump_domain_organisation,
             args=(ipr_pro_url, df.proteins, df.protein2matches,
                   df.protein2domorg),
             kwargs=dict(tempdir=temp_dir, workers=8),
             name="export-dom-orgs",
             requires=["export-matches"],
             scheduler=dict(cpu=8, mem=16000, scratch=50000, queue=lsf_queue)),
        Task(fn=interpro.oracle.entries.dump_similar_entries,
             args=(ipr_pro_url, df.protein2matches, df.overlapping_entries),
             name="export-sim-entries",
             requires=["export-matches"],
             scheduler=dict(mem=16000, queue=lsf_queue)),

        # Data from InterPro + other sources (Pfam, PDBe)
        Task(fn=interpro.oracle.clans.export_clans,
             args=(ipr_pro_url, pfam_url, df.clans, df.alignments),
             name="export-clans",
             scheduler=dict(mem=8000, queue=lsf_queue)),
        Task(fn=pdbe.export_structures,
             args=(ipr_pro_url, pdbe_url, df.structures),
             name="export-structures",
             scheduler=dict(mem=8000, queue=lsf_queue)),

        # Data from UniProt
        Task(fn=uniprot.proteomes.export_proteomes,
             args=(uniprot_url, df.proteomes),
             name="export-reference-proteomes",
             scheduler=dict(mem=4000, queue=lsf_queue)),
        Task(fn=uniprot.proteins.export_entry2evidence,
             args=(uniprot_url, df.proteins, df.protein2evidence),
             kwargs=dict(tempdir=temp_dir, workers=8),
             name="export-evidences",
             requires=["export-proteins"],
             scheduler=dict(cpu=8, mem=16000, scratch=50000, queue=lsf_queue)),
        Task(fn=uniprot.proteins.export_entry2functions,
             args=(uniprot_url, df.proteins, df.protein2function),
             kwargs=dict(tempdir=temp_dir, workers=8),
             name="export-functions",
             requires=["export-proteins"],
             scheduler=dict(cpu=8, mem=16000, scratch=50000, queue=lsf_queue)),
        Task(fn=uniprot.proteins.export_entry2name,
             args=(uniprot_url, df.proteins, df.protein2name),
             kwargs=dict(tempdir=temp_dir, workers=8),
             name="export-names",
             requires=["export-proteins"],
             scheduler=dict(cpu=8, mem=16000, scratch=50000, queue=lsf_queue)),
        Task(fn=uniprot.proteins.export_entry2proteome,
             args=(uniprot_url, df.proteins, df.protein2proteome),
             kwargs=dict(tempdir=temp_dir, workers=8),
             name="export-proteomes",
             requires=["export-proteins"],
             scheduler=dict(cpu=8, mem=16000, scratch=50000, queue=lsf_queue)),
    ]

    return tasks
