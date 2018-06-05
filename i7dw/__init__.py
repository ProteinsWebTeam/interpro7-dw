#!/usr/bin/env python
# -*- coding: utf-8 -*-

__version__ = '0.2.0'

import logging
import os

from . import disk
from .ebi import pdbe


logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s: %(levelname)s: %(message)s',
    datefmt='%y-%m-%d %H:%M:%S'
)


def export(uri, comments, descriptions, evidences, genes, proteomes, annotations, struct_matches, prot_matches,
           residues, proteins, dst, chunk_size=100000):
    logging.info('getting PDBe structures')
    structures = pdbe.get_structures(uri, citations=True, fragments=True, by_protein=True)

    _comments = disk.Store(comments)
    _descriptions = disk.Store(descriptions)
    _evidences = disk.Store(evidences)
    _genes = disk.Store(genes)
    _proteomes = disk.Store(proteomes)
    _annotations = disk.Store(annotations)
    _struct_matches = disk.Store(struct_matches)
    _prot_matches = disk.Store(prot_matches)
    _residues = disk.Store(residues)
    _proteins = disk.Store(proteins)

    logging.info('merging data')

    with disk.File(dst) as fh:
        proteins = {}
        cnt = 1
        for accession, protein in _proteins.iter():
            # Enrich with data from other stores
            protein.update({
                'comments': _comments.get(accession, []),
                'descriptions': _descriptions.get(accession, (None, None)),
                'evidence': _evidences.get(accession),
                'gene': _genes.get(accession),
                'proteomes': _proteomes.get(accession, []),
                'goa': _annotations.get(accession, []),
                'struct_matches': _struct_matches.get(accession, {}),
                'prot_matches': _prot_matches.get(accession, []),
                'residues': _residues.get(accession, {}),
                'structures': structures.get(accession, {})
            })

            proteins[accession.lower()] = protein
            cnt += 1

            if len(proteins) == chunk_size:
                fh.write(proteins)
                proteins = {}

            if not cnt % 1000000:
                logging.info('{:>12}'.format(cnt))

        fh.write(proteins)
        logging.info('{:>12}'.format(cnt))

    for f in (comments, descriptions, evidences, genes, proteomes, annotations, struct_matches, prot_matches,
              residues, proteins):
        os.unlink(f)
