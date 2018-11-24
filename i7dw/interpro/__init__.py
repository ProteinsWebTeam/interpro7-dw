#!/usr/bin/env python
# -*- coding: utf-8 -*-

from .elastic import (
    init_dir as init_elastic,
    create_documents,
    index_documents,
    update_alias
)

from .export import (
    chunk_proteins,
    export_protein2matches,
    export_protein2features,
    export_protein2residues,
    export_proteins,
    export_sequences
)

from .mysql import (
    init_tables,
    insert_taxa,
    insert_proteomes,
    insert_databases,
    insert_entries,
    insert_annotations,
    insert_structures,
    insert_sets,
    insert_proteins,
    make_release_notes,
    get_entries,
    get_entry_databases,
    get_sets,
    update_counts
)

from .supermatch import calculate_relationships

from .xref import export as export_xrefs
