#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import configparser
import os

from mundone import Task, Workflow

from i7dw import elastic, mysql
from i7dw.ebi import ebisearch, goa, interpro, uniprot


def parse_elastic_hosts(str_hosts):
    hosts = []
    for host in str_hosts.split(','):
        host = host.strip()

        if not host:
            continue

        pair = host.split(':')
        if len(pair) == 2:
            host = pair[0]
            port = int(pair[1])
        else:
            host = host
            port = 9200

        hosts.append({'host': host, 'port': port})

    return hosts


def cli():
    parser = argparse.ArgumentParser(description='Build InterPro7 data warehouse')
    parser.add_argument('config', metavar='config.ini', help='configuration file')
    parser.add_argument('-t', '--tasks', nargs='*', help='tasks to run')
    parser.add_argument('-l', '--list', action='store_true', default=False,
                        help='list steps that would be processed, but do not process them')
    parser.add_argument('--nodep', action='store_true', default=False,
                        help='do not include dependencies (run only the requested tasks)')
    parser.add_argument('-d', '--detach', action='store_true', default=False,
                        help='do not wait for tasks to complete (only tasks without dependencies are run)')
    args = parser.parse_args()

    if not os.path.isfile(args.config):
        raise RuntimeError("cannot open '{}': no such file or directory".format(args.config))

    config = configparser.ConfigParser()
    config.read(args.config)

    export_dir = config['export']['dir']
    if not os.path.isdir(export_dir):
        os.makedirs(export_dir)

    interpro_oracle = config['databases']['interpro_oracle']
    interpro_mysql = config['databases']['interpro_mysql']
    pfam_mysql = config['databases']['pfam_mysql']

    tasks = [
        Task(
            name='export_comments',
            fn=uniprot.export_protein_comments,
            args=(interpro_oracle, os.path.join(export_dir, 'comments.bs')),
            lsf=dict(queue=config['workflow']['queue'], mem=1000)
        ),
        Task(
            name='export_descriptions',
            fn=uniprot.export_protein_descriptions,
            args=(interpro_oracle, os.path.join(export_dir, 'descriptions.bs')),
            lsf=dict(queue=config['workflow']['queue'], mem=2000)
        ),
        Task(
            name='export_evidences',
            fn=uniprot.export_protein_evidence,
            args=(interpro_oracle, os.path.join(export_dir, 'evidences.bs')),
            lsf=dict(queue=config['workflow']['queue'], mem=1000)
        ),
        Task(
            name='export_genes',
            fn=uniprot.export_protein_gene,
            args=(interpro_oracle, os.path.join(export_dir, 'genes.bs')),
            lsf=dict(queue=config['workflow']['queue'], mem=1000)
        ),
        Task(
            name='export_proteomes',
            fn=uniprot.export_protein_proteomes,
            args=(interpro_oracle, os.path.join(export_dir, 'proteomes.bs')),
            lsf=dict(queue=config['workflow']['queue'], mem=1000)
        ),
        Task(
            name='export_annotations',
            fn=goa.export_annotations,
            args=(interpro_oracle, os.path.join(export_dir, 'annotations.bs')),
            lsf=dict(queue=config['workflow']['queue'], mem=8000)
        ),
        Task(
            name='export_struct_matches',
            fn=interpro.export_struct_matches,
            args=(interpro_oracle, os.path.join(export_dir, 'struct_matches.bs')),
            lsf=dict(queue=config['workflow']['queue'], mem=4000)
        ),
        Task(
            name='export_prot_matches',
            fn=interpro.export_prot_matches,
            args=(interpro_oracle, os.path.join(export_dir, 'prot_matches.bs')),
            lsf=dict(queue=config['workflow']['queue'], mem=16000)
        ),
        Task(
            name='export_residues',
            fn=interpro.export_residues,
            args=(interpro_oracle, os.path.join(export_dir, 'residues.bs')),
            lsf=dict(queue=config['workflow']['queue'], mem=4000)
        ),
        Task(
            name='export_proteins',
            fn=interpro.export_proteins,
            args=(interpro_oracle, os.path.join(export_dir, 'proteins.bs')),
            lsf=dict(queue=config['workflow']['queue'], mem=4000)
        ),

        Task(
            name='init_mysql',
            fn=mysql.init,
            args=(interpro_mysql,),
            lsf=dict(queue=config['workflow']['queue'])
        ),
        Task(
            name='insert_taxa',
            fn=mysql.insert_taxa,
            args=(interpro_oracle, interpro_mysql),
            lsf=dict(queue=config['workflow']['queue'], mem=4000),
            requires=['init_mysql'],
        ),
        Task(
            name='insert_proteomes',
            fn=mysql.insert_proteomes,
            args=(interpro_oracle, interpro_mysql),
            lsf=dict(queue=config['workflow']['queue'], mem=1000),
            requires=['insert_taxa'],
        ),

        Task(
            name='insert_databases',
            fn=mysql.insert_databases,
            args=(interpro_oracle, interpro_mysql),
            lsf=dict(queue=config['workflow']['queue']),
            requires=['init_mysql'],
        ),

        Task(
            name='insert_entries',
            fn=mysql.insert_entries,
            args=(interpro_oracle, pfam_mysql, interpro_mysql),
            lsf=dict(queue=config['workflow']['queue'], mem=16000),
            requires=['insert_databases'],
        ),
        Task(
            name='insert_annotations',
            fn=mysql.insert_annotations,
            args=(pfam_mysql, interpro_mysql),
            lsf=dict(queue=config['workflow']['queue'], mem=16000),
            requires=['insert_entries'],
        ),

        Task(
            name='insert_structures',
            fn=mysql.insert_structures,
            args=(interpro_oracle, interpro_mysql),
            lsf=dict(queue=config['workflow']['queue'], mem=16000),
            requires=['insert_databases'],
        ),

        Task(
            name='insert_sets',
            fn=mysql.insert_sets,
            args=(pfam_mysql, interpro_mysql),
            lsf=dict(queue=config['workflow']['queue']),
            requires=['insert_databases'],
        ),

        Task(
            name='insert_proteins',
            fn=mysql.insert_proteins,
            args=(
                interpro_mysql,
                os.path.join(export_dir, 'proteins.bs'),
                os.path.join(export_dir, 'evidences.bs'),
                os.path.join(export_dir, 'descriptions.bs'),
                os.path.join(export_dir, 'comments.bs'),
                os.path.join(export_dir, 'proteomes.bs'),
                os.path.join(export_dir, 'genes.bs'),
                os.path.join(export_dir, 'annotations.bs'),
                os.path.join(export_dir, 'residues.bs'),
                os.path.join(export_dir, 'struct_matches.bs')
            ),
            lsf=dict(queue=config['workflow']['queue'], mem=16000),
            requires=[
                'insert_taxa', 'export_proteins', 'export_evidences', 'export_descriptions', 'export_comments',
                'export_proteomes', 'export_genes', 'export_annotations', 'export_residues', 'export_struct_matches'
            ],
        ),

        Task(
            name='index_relationships',
            fn=elastic.index_relationships,
            args=(
                interpro_oracle,
                interpro_mysql,
                os.path.join(export_dir, 'proteins.bs'),
                os.path.join(export_dir, 'descriptions.bs'),
                os.path.join(export_dir, 'comments.bs'),
                os.path.join(export_dir, 'proteomes.bs'),
                os.path.join(export_dir, 'prot_matches.bs'),
                config['elastic']['dir']
            ),
            kwargs=dict(
                producers=3,
                loaders=4,
                supermatches=True,
                hosts=parse_elastic_hosts(config['elastic']['hosts']),
                doc_type=config['elastic']['type'],
                suffix=config['meta']['release'],
                indices_json=config['elastic']['indices'],
                properties_json=config['elastic']['properties'],
                shards=config.getint('elastic', 'shards')
            ),
            lsf=dict(queue=config['workflow']['queue'], mem=48000, cpu=9),
            requires=[
                'insert_entries', 'insert_sets', 'insert_taxa',
                'export_proteins', 'export_descriptions', 'export_comments', 'export_proteomes', 'export_prot_matches'
            ],
        ),

        Task(
            name='collect',
            fn=elastic.collect,
            args=(
                interpro_mysql,
                parse_elastic_hosts(config['elastic']['hosts']),
                config['elastic']['type'],
                config['elastic']['dir']
            ),
            kwargs=dict(
                suffix=config['meta']['release'],
                processes=4,
                gzip=True
            ),
            lsf=dict(queue=config['workflow']['queue'], mem=16000, cpu=5),
            requires=['index_relationships'],
        ),

        Task(
            name='update_alias',
            fn=elastic.update_alias,
            args=(
                interpro_mysql,
                parse_elastic_hosts(config['elastic']['hosts']),
                config['elastic']['alias'],
                config['meta']['release']
            ),
            lsf=dict(queue=config['workflow']['queue']),
            requires=['collect'],
        ),

        Task(
            name='ebisearch',
            fn=ebisearch.export,
            args=(
                interpro_mysql,
                os.path.join(export_dir, 'proteins.bs'),
                os.path.join(export_dir, 'prot_matches.bs'),
                os.path.join(export_dir, 'struct_matches.bs'),
                os.path.join(export_dir, 'proteomes.bs'),
                config['meta']['name'],
                config['meta']['release'],
                config['meta']['release_date'],
                config['ebisearch']['dir'],
            ),
            lsf=dict(queue=config['workflow']['queue'], mem=24000, tmp=10000),
            requires=[
                'insert_entries', 'export_proteins', 'export_prot_matches', 'export_struct_matches',
                'export_proteomes',
            ],
        ),
    ]

    if args.detach:
        secs = 0
        cascade_kill = False
    else:
        secs = 10
        cascade_kill = True

    w = Workflow(tasks, dir=config['workflow']['dir'], db=config['workflow']['db'], cascade_kill=cascade_kill)
    w.run(args.tasks, process=(not args.list), incdep=(not args.nodep), secs=secs)
