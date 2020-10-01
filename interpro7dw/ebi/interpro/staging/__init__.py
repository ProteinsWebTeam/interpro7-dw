# -*- coding: utf-8 -*-

from .annotations import insert_annotations
from .database import insert_databases, make_release_notes
from .clan import init_clans, insert_clans
from .entry import insert_entries
from .protein import export_ida
from .protein import export_uniprot2entries
from .protein import insert_isoforms
from .protein import insert_proteins
from .proteome import insert_proteomes
from .structure import insert_structures
from .taxonomy import insert_taxonomy
from .utils import drop_database
