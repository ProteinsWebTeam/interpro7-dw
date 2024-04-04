import os
import pickle
import shelve
import shutil
from copy import deepcopy
from multiprocessing import Process, Queue
from tempfile import mkdtemp, mkstemp

from interpro7dw.interpro.oracle.entries import Entry
from interpro7dw.utils import logger
from interpro7dw.utils.store import Directory, KVStore
from . import config


def export_documents(proteins_file: str, matches_file: str, domorgs_file: str,
                     protein2proteome_file: str, uniprot2pdb_file: str,
                     pdbmatches_file: str, alphafold_file: str,
                     proteomes_file: str, structures_file: str,
                     clans_file: str, entries_file: str, taxa_file: str,
                     outdirs: list[str], version: str,
                     processes: int = 8, tempdir: str | None = None):
    logger.info("starting")
    if tempdir:
        # Ensure the directory exists
        os.makedirs(tempdir, exist_ok=True)

    # Create a subdirectory where all temporary files will be created
    tempdir = mkdtemp(dir=tempdir)

    inqueue = Queue()
    outqueue = Queue()
    workers = []
    for _ in range(max(1, processes - 1)):
        p = Process(
            target=gen_rel_docs,
            args=(proteins_file, matches_file, domorgs_file,
                  protein2proteome_file, alphafold_file, proteomes_file,
                  uniprot2pdb_file, pdbmatches_file, version,
                  inqueue, outqueue, tempdir)
        )
        p.start()
        workers.append(p)

    logger.info("creating directories")
    directories = []
    for path in outdirs:
        if os.path.isdir(path):
            shutil.rmtree(path)

        os.makedirs(path, mode=0o775)
        directories.append(Directory(root=path))
        open(os.path.join(path, f"{version}{config.LOAD_SUFFIX}"), "w")

    logger.info("dump entries/taxonomy/structures properties")
    props_file = prepare_props(taxa_file, entries_file, clans_file,
                               structures_file, tempdir)
    for _ in workers:
        inqueue.put(props_file)

    # Submit tasks
    n_tasks = 0
    with KVStore(proteins_file) as proteins:
        keys = proteins.get_keys()
        for i, start in enumerate(keys):
            try:
                stop = keys[i+1]
            except IndexError:
                stop = None

            inqueue.put((start, stop))
            n_tasks += 1

        for _ in range(len(workers)):
            inqueue.put(None)

    logger.info(f"{n_tasks:,} tasks submitted")
    running = len(workers)
    n_done = n_documents = 0
    step = milestone = 5
    seen_entries = set()
    seen_structures = set()
    seen_taxa = set()
    while running:
        done, file, n = outqueue.get()
        if done:
            running -= 1

            # Load structures/taxa observed by the child process
            with open(file, "rb") as fh:
                seen_entries |= pickle.load(fh)
                seen_structures |= pickle.load(fh)
                seen_taxa |= pickle.load(fh)

            os.unlink(file)
            continue

        move_docs_file(file, directories)

        n_done += 1
        n_documents += n
        progress = n_done / n_tasks * 100
        if progress >= milestone:
            logger.info(f"{progress:>6.1f}% ({n_documents:,} documents)")
            milestone += step

    for p in workers:
        p.join()

    """
    Add additional entry-structure pairs
    i.e. entries with PDB matches where the PDB structure is not associated
    to a protein in UniProtKB
    """
    logger.info("creating additional relationship documents")
    with open(entries_file, "rb") as fh:
        entries = pickle.load(fh)

    member2clan = {}
    with open(clans_file, "rb") as fh:
        for clan in pickle.load(fh).values():
            for entry_acc, _, _, _, _ in clan["members"]:
                member2clan[entry_acc] = (clan["accession"], clan["name"])

    with open(structures_file, "rb") as fh:
        structures = pickle.load(fh)

    with open(taxa_file, "rb") as fh:
        taxa = pickle.load(fh)

    documents = []
    with shelve.open(pdbmatches_file, writeback=False) as pdb_matches:
        for pdb_chain, pdb_entry in pdb_matches.items():
            if pdb_chain in seen_structures:
                continue

            pdb_id, chain_id = pdb_chain.split("_")

            try:
                structure = structures[pdb_id]
            except KeyError:
                continue

            for entry_acc, match in pdb_entry["matches"].items():
                entry = entries[entry_acc]

                if entry.deletion_date or not entry.public:
                    continue

                if entry.integrated_in:
                    integrated_in = entry.integrated_in.lower()
                else:
                    integrated_in = None

                seen_entries.add(entry_acc)
                doc = init_rel_doc()
                doc.update({
                    "entry_acc": entry.accession.lower(),
                    "entry_db": entry.database.lower(),
                    "entry_type": entry.type.lower(),
                    "entry_date": entry.creation_date.strftime("%Y-%m-%d"),
                    "entry_protein_locations": [],
                    "entry_go_terms": [t["identifier"] for t in entry.go_terms],
                    "entry_integrated": integrated_in,
                    "text_entry": join(entry.accession, entry.short_name,
                                       entry.name,
                                       entry.type.lower(), integrated_in),
                    "structure_acc": pdb_id.lower(),
                    "structure_resolution": structure["resolution"],
                    "structure_date": structure["date"],
                    "structure_evidence": structure["evidence"],
                    "text_structure": join(pdb_id,
                                           structure["evidence"],
                                           structure["name"]),
                    "structure_chain_acc": chain_id,
                    "structure_chain": pdb_chain,
                    "structure_protein_length": pdb_entry["length"],
                    "entry_structure_locations": match["locations"]
                })

                if entry_acc in member2clan:
                    clan_acc, clan_name = member2clan[entry_acc]
                    doc.update({
                        "set_acc": clan_acc.lower(),
                        "set_db": entry.database.lower(),
                        "text_set": join(clan_acc, clan_name),
                    })

                taxon_id = structure["taxonomy"].get(chain_id)
                if taxon_id and taxon_id in taxa:
                    taxon = taxa[taxon_id]
                    seen_taxa.add(taxon_id)

                    doc.update({
                        "tax_id": taxon_id,
                        "tax_name": taxon["sci_name"],
                        "tax_lineage": taxon["lineage"],
                        "tax_rank": taxon["rank"],
                        "text_taxonomy": join(taxon_id, taxon["full_name"],
                                              taxon["rank"])
                    })

                documents.append((
                    get_rel_doc_index(doc),
                    get_rel_doc_id(doc),
                    doc
                ))

    # Entry documents
    for entry in entries.values():
        if (entry.accession in seen_entries or
                entry.deletion_date or
                not entry.public):
            continue

        if entry.integrated_in:
            integrated_in = entry.integrated_in.lower()
        else:
            integrated_in = None

        doc = init_rel_doc()
        doc.update({
            "entry_acc": entry.accession.lower(),
            "entry_db": entry.database.lower(),
            "entry_type": entry.type.lower(),
            "entry_date": entry.creation_date.strftime("%Y-%m-%d"),
            "entry_protein_locations": [],
            "entry_go_terms": [t["identifier"] for t in entry.go_terms],
            "entry_integrated": integrated_in,
            "text_entry": join(entry.accession, entry.short_name,
                               entry.name,
                               entry.type.lower(), integrated_in),
        })

        if entry.accession in member2clan:
            clan_acc, clan_name = member2clan[entry.accession]
            doc.update({
                "set_acc": clan_acc.lower(),
                "set_db": entry.database.lower(),
                "text_set": join(clan_acc, clan_name),
            })

        documents.append((
            get_rel_doc_index(doc),
            get_rel_doc_id(doc),
            doc
        ))

    # Taxon documents
    for taxon_id, taxon in taxa.items():
        if taxon_id in seen_taxa:
            continue

        doc = init_rel_doc()
        doc.update({
            "tax_id": taxon["id"],
            "tax_name": taxon["full_name"],
            "tax_lineage": taxon["lineage"],
            "tax_rank": taxon["rank"],
            "text_taxonomy": join(taxon["id"], taxon["full_name"],
                                  taxon["rank"])
        })

        documents.append((
            get_rel_doc_index(doc),
            get_rel_doc_id(doc),
            doc
        ))

    n_documents += len(documents)
    for i in range(0, len(documents), 100000):
        fd, file = mkstemp(dir=tempdir)
        with open(fd, "wb") as fh:
            pickle.dump(documents[i:i+100000], fh)

        move_docs_file(file, directories)

    documents.clear()

    logger.info("creating domain architecture documents")
    documents = gen_ida_docs(domorgs_file, entries, version)
    n_documents += len(documents)
    for i in range(0, len(documents), 100000):
        fd, file = mkstemp(dir=tempdir)
        with open(fd, "wb") as fh:
            pickle.dump(documents[i:i+100000], fh)

        move_docs_file(file, directories)

    documents.clear()
    shutil.rmtree(tempdir)

    for path in outdirs:
        open(os.path.join(path, f"{version}{config.DONE_SUFFIX}"), "w").close()

    logger.info(f"done ({n_documents:,} documents)")


def prepare_props(taxa_file: str, entries_file: str, clans_file: str,
                  structures_file: str, outdir: str) -> str:
    taxa = {}
    with open(taxa_file, "rb") as fh:
        for taxon_id, taxon in pickle.load(fh).items():
            taxa[taxon_id] = {
                "tax_id": taxon_id,
                "tax_lineage": taxon["lineage"],
                "tax_name": taxon["sci_name"],
                "tax_rank": taxon["rank"],
                "text_taxonomy": join(taxon_id, taxon["full_name"],
                                      taxon["rank"])
            }

    member2clan = {}
    with open(clans_file, "rb") as fh:
        for clan in pickle.load(fh).values():
            for entry_acc, _, _, _, _ in clan["members"]:
                member2clan[entry_acc] = (clan["accession"], clan["name"])

    entries = {}
    with open(entries_file, "rb") as fh:
        for entry_acc, entry in pickle.load(fh).items():
            entries[entry_acc] = {
                "entry_acc": entry_acc.lower(),
                "entry_db": entry.database.lower(),
                "entry_type": entry.type.lower(),
                "entry_date": entry.creation_date.strftime("%Y-%m-%d"),
                "entry_go_terms": [t["identifier"] for t in
                                   entry.go_terms],
                "entry_integrated": (entry.integrated_in.lower()
                                     if entry.integrated_in else None),
                "text_entry": join(entry_acc, entry.short_name, entry.name,
                                   entry.type.lower(), entry.integrated_in),
            }

            if entry_acc in member2clan:
                clan_acc, clan_name = member2clan[entry_acc]
                entries[entry_acc].update({
                    "set_acc": clan_acc.lower(),
                    "set_db": entry.database.lower(),
                    "text_set": join(clan_acc, clan_name),
                })

    structures = {}
    with open(structures_file, "rb") as fh:
        for pdb_id, entry in pickle.load(fh).items():
            structures[pdb_id] = {
                "structure_acc": pdb_id.lower(),
                "structure_resolution": entry["resolution"],
                "structure_date": entry["date"],
                "structure_evidence": entry["evidence"],
                "text_structure": join(pdb_id, entry["evidence"], entry["name"])
            }

    fd, file = mkstemp(dir=outdir)
    with open(fd, "wb") as fh:
        pickle.dump(entries, fh)
        pickle.dump(structures, fh)
        pickle.dump(taxa, fh)

    return file


def init_rel_doc() -> dict:
    return {key: None for key in config.REL_BODY["mappings"]["properties"]}


def join(*args, separator: str = " ") -> str:
    items = []

    for item in args:
        if item is None:
            continue
        elif isinstance(item, (int, float)):
            item = str(item)
        elif isinstance(item, (list, set, tuple)):
            item = separator.join(map(str, item))
        elif isinstance(item, dict):
            item = separator.join(map(str, item.values()))
        elif not isinstance(item, str):
            continue

        items.append(item)

    return separator.join(items)


def get_rel_doc_id(doc: dict) -> str:
    return join(doc["protein_acc"],
                doc["proteome_acc"],
                doc["entry_acc"],
                doc["set_acc"],
                doc["structure_acc"],
                doc["structure_chain_acc"],
                doc["tax_id"], separator="-")


def get_rel_doc_index(doc: dict) -> str:
    return doc["entry_db"] or config.REL_DEFAULT_INDEX


def gen_ida_docs(domorgs_file: str, entries: dict[str, Entry],
                 version: str) -> list[tuple[str, str, dict]]:
    documents = []
    seen_domains = set()
    with KVStore(domorgs_file) as domains:
        for protein_acc, domain in domains.items():
            domain_id = domain["id"]
            if domain_id not in seen_domains:
                seen_domains.add(domain_id)
                locations = []
                for loc in domain["locations"]:
                    entry_acc = loc["pfam"]
                    locations.append({
                        "accession": entry_acc,
                        "name": entries[entry_acc].short_name,
                        "coordinates": [{
                            "fragments": [{
                                "start": loc["start"],
                                "end": loc["end"]
                            }]
                        }]
                    })

                    if entry_acc := loc["interpro"]:
                        locations.append({
                            "accession": entry_acc,
                            "name": entries[entry_acc].short_name,
                            "coordinates": [{
                                "fragments": [{
                                    "start": loc["start"],
                                    "end": loc["end"]
                                }]
                            }]
                        })

                documents.append((
                    config.IDA_INDEX + version,
                    domain["id"],
                    {
                        "ida_id": domain["id"],
                        "ida": domain["key"],
                        "representative": {
                            "accession": domain["protein"],
                            "length": domain["length"],
                            "domains": locations
                        },
                        "counts": domain["count"]
                    }
                ))

    return documents


def gen_rel_docs(proteins_file: str, matches_file: str,
                 domorgs_file: str, protein2proteome_file: str,
                 alphafold_file: str, proteomes_file: str,
                 uniprot2pdb_file: str, pdb_matches_file: str, version: str,
                 inqueue: Queue, outqueue: Queue, outdir: str):
    props_file = inqueue.get()
    with open(props_file, "rb") as fh:
        entries = pickle.load(fh)
        structures = pickle.load(fh)
        taxa = pickle.load(fh)

    with open(proteomes_file, "rb") as fh:
        proteomes = pickle.load(fh)

    with open(uniprot2pdb_file, "rb") as fh:
        uniprot2pdb = pickle.load(fh)

    pdb_matches = shelve.open(pdb_matches_file, writeback=False)
    proteins_store = KVStore(proteins_file)
    matches_store = KVStore(matches_file)
    proteomes_store = KVStore(protein2proteome_file)
    alphafold_store = KVStore(alphafold_file)
    domorgs_store = KVStore(domorgs_file)
    seen_entries = set()
    seen_structures = set()
    seen_taxa = set()
    for start, stop in iter(inqueue.get, None):
        documents = []
        for protein_acc, protein in proteins_store.range(start, stop):
            taxon_id = protein["taxid"]
            taxon = taxa[taxon_id]
            seen_taxa.add(taxon_id)
            protein_docs = 0

            af_models = alphafold_store.get(protein_acc, [])
            if af_models:
                # af_models: sorted list of tuples (AFDB ID, pLDDT score)
                # sorted by pLDDT score (ascending order)
                af_score = af_models[-1][1]
            else:
                af_score = -1

            # Creates an empty document (all properties set to None)
            base_doc = init_rel_doc()

            # Set protein-related properties
            base_doc.update({
                # Protein
                "protein_acc": protein_acc.lower(),
                "protein_length": protein["length"],
                "protein_is_fragment": protein["fragment"],
                "protein_af_score": af_score,
                "protein_db": ("reviewed" if protein["reviewed"]
                               else "unreviewed"),
                "text_protein": join(protein_acc,
                                     protein["identifier"],
                                     taxon["tax_name"]),
                # Taxonomy
                **taxa[taxon_id]
            })

            proteome_id = proteomes_store.get(protein_acc)
            if proteome_id:
                # Adds proteome-related properties
                proteome = proteomes[proteome_id]
                base_doc.update({
                    "proteome_acc": proteome_id.lower(),
                    "proteome_name": proteome["name"],
                    "proteome_is_reference": proteome["is_reference"],
                    "text_proteome": join(proteome_id,
                                          proteome["name"],
                                          proteome["assembly"],
                                          proteome["taxon_id"],
                                          proteome["strain"]),
                })

            # Get domain architecture, if any
            try:
                domain = domorgs_store[protein_acc]
            except KeyError:
                domain_id = domain_str = None
                domain_members = set()
            else:
                domain_id = domain["id"]
                domain_str = domain["key"]
                domain_members = domain["members"]

            # Get structures mapped the protein, and structural matches
            protein_structures = uniprot2pdb.get(protein_acc, {})
            pdb_props = {}         # Structure-related properties
            struct_entries = {}    # Structural entry matches
            structs_w_rel = set()  # structures with relationships
            for pdb_chain, segments in protein_structures.items():
                pdb_id, chain_id = pdb_chain.split("_")
                try:
                    structure = structures[pdb_id]
                    pdb_entry = pdb_matches[pdb_chain]
                except KeyError:
                    continue

                locations = []
                for segment in segments:
                    locations.append({
                        "fragments": [{
                            # Coord. of UniProt entry on the PDB seq
                            "start": segment["structure_start"],
                            "end": segment["structure_end"],
                            # Coord. of PDB entry on the UniProt seq
                            "protein_start": segment["protein_start"],
                            "protein_end": segment["protein_end"],
                        }]
                    })

                pdb_props[pdb_chain] = {
                    # Structure/chain info
                    **structure,
                    # Used to map a protein to a structure
                    "structure_protein_acc": protein_acc.lower(),
                    "structure_protein_length": pdb_entry["length"],
                    # SIFTS mapping coordinates
                    "structure_protein_locations": locations,
                }

                for entry_acc, match in pdb_entry["matches"].items():
                    locations = match["locations"]

                    try:
                        struct_entries[entry_acc][pdb_chain] = locations
                    except KeyError:
                        struct_entries[entry_acc] = {pdb_chain: locations}

            # Get (UniProt) protein matches
            s_matches, e_matches = matches_store.get(protein_acc, ({}, {}))
            prot_entries = {**s_matches, **e_matches}

            accessions = set(prot_entries.keys()) | set(struct_entries.keys())
            for entry_acc in accessions:
                seen_entries.add(entry_acc)
                entry_doc = deepcopy(base_doc)
                entry_doc.update(entries[entry_acc])

                if entry_acc in domain_members:
                    entry_doc.update({
                        "ida_id": domain_id,
                        "ida": domain_str,
                    })

                if entry_acc in prot_entries:
                    match = prot_entries[entry_acc]
                    locations = match["locations"]
                    entry_doc["entry_protein_locations"] = locations
                else:
                    # Clear props used for protein-entry relationships
                    entry_doc.update({
                        "protein_acc": None,
                        "protein_is_fragment": None,
                        "protein_af_score": None,
                        "text_protein": None
                    })

                if entry_acc in struct_entries:
                    # Entry has structural matches: one doc per struct/chain
                    for pdb_chain, locations in struct_entries[entry_acc].items():
                        structs_w_rel.add(pdb_chain)
                        doc = deepcopy(entry_doc)
                        doc.update({
                            **pdb_props[pdb_chain],
                            "entry_structure_locations": locations
                        })
                        documents.append((
                            get_rel_doc_index(doc) + version,
                            get_rel_doc_id(doc),
                            doc
                        ))
                        protein_docs += 1
                else:
                    # No entry-structure relationship
                    documents.append((
                        get_rel_doc_index(entry_doc) + version,
                        get_rel_doc_id(entry_doc),
                        entry_doc
                    ))
                    protein_docs += 1

            # Add protein-structure relationships
            for pdb_chain, props in pdb_props.items():
                seen_structures.add(pdb_chain)
                if pdb_chain in structs_w_rel:
                    continue

                struct_doc = deepcopy(base_doc)
                struct_doc.update({
                    **props,
                    # Clear props used for protein-entry relationships
                    "protein_acc": None,
                    "protein_is_fragment": None,
                    "protein_af_score": None,
                    "text_protein": None
                })
                documents.append((
                    get_rel_doc_index(struct_doc) + version,
                    get_rel_doc_id(struct_doc),
                    struct_doc
                ))
                protein_docs += 1

            if not protein_docs:
                # Not protein matches: simple protein document
                documents.append((
                    get_rel_doc_index(base_doc) + version,
                    get_rel_doc_id(base_doc),
                    base_doc
                ))

        fd, file = mkstemp(dir=outdir)
        with open(fd, "wb") as fh:
            pickle.dump(documents, fh)

        outqueue.put((False, file, len(documents)))

    pdb_matches.close()
    proteins_store.close()
    matches_store.close()
    proteomes_store.close()
    alphafold_store.close()
    domorgs_store.close()

    # Send seen structures and taxa to parent process
    fd, file = mkstemp(dir=outdir)
    with open(fd, "wb") as fh:
        pickle.dump(seen_entries, fh)
        pickle.dump(seen_structures, fh)
        pickle.dump(seen_taxa, fh)

    outqueue.put((True, file, 0))


def move_docs_file(src: str, directories: list[Directory]):
    for dir_obj in directories:
        dst = dir_obj.mktemp(createfile=False)
        shutil.copy(src, dst)
        os.rename(dst, f"{dst}{config.EXTENSION}")

    os.unlink(src)
