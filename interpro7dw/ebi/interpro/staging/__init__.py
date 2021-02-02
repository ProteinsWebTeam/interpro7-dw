# -*- coding: utf-8 -*-

from .annotations import insert_annotations
from .database import insert_databases, insert_release_notes
from .clan import insert_clans
from .entry import insert_entries
from .protein import insert_isoforms
from .protein import insert_extra_features, insert_proteins, insert_residues
from .proteome import insert_proteomes
from .structure import insert_structures, insert_structural_models
from .taxonomy import insert_taxonomy
from .utils import drop_database
