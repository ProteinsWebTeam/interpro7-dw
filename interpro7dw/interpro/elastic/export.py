import os
import pickle
import shelve
import shutil
from copy import deepcopy
from multiprocessing import Process, Queue
from tempfile import mkstemp

from interpro7dw.utils import logger
from interpro7dw.utils.store import Directory, KVStore
from . import config


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


def gen_ida_docs(
        domorgs_file: str,
        entries_file: str,
        version: str
) -> list[tuple[str, str, dict]]:
    short_names = {}
    with open(entries_file, "rb") as fh:
        for acc, entry in pickle.load(fh).items():
            short_names[acc] = entry.short_name

    domorgs_store = KVStore(domorgs_file)
    seen_domains = set()
    documents = []
    for protein_acc, domain in domorgs_store.items():
        domain_id = domain["id"]
        if domain_id not in seen_domains:
            seen_domains.add(domain_id)
            locations = []
            for loc in domain["locations"]:
                entry_acc = loc["pfam"]
                locations.append({
                    "accession": entry_acc,
                    "name": short_names[entry_acc],
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
                        "name": short_names[entry_acc],
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

    domorgs_store.close()
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

                    if match["database"].lower() == "panther":
                        """
                        PANTHER: remove the node ID 
                        (other databases do not have a subfamily property)
                        """
                        for loc in locations:
                            try:
                                del loc["subfamily"]["node"]
                            except KeyError:
                                continue  # No subfamily annotation

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
                    "set_db": entry.database,
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


def export_mp(proteins_file: str, matches_file: str, domorgs_file: str,
              protein2proteome_file: str, uniprot2pdb_file: str,
              pdbmatches_file: str, alphafold_file: str, proteomes_file: str,
              structures_file: str, clans_file: str, entries_file: str,
              taxa_file: str, outdirs: list[str], tmpdir: str,
              version: str, processes: int = 8):
    logger.info("starting")
    inqueue = Queue()
    outqueue = Queue()
    workers = []
    for _ in range(max(1, processes - 1)):
        p = Process(
            target=gen_rel_docs,
            args=(proteins_file, matches_file, domorgs_file,
                  protein2proteome_file, alphafold_file, proteomes_file,
                  uniprot2pdb_file, pdbmatches_file, version,
                  inqueue, outqueue, tmpdir)
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

    if os.path.isdir(tmpdir):
        shutil.rmtree(tmpdir)

    os.makedirs(tmpdir, mode=0o775)

    logger.info("dump entries/taxonomy/structures properties")
    props_file = prepare_props(taxa_file, entries_file, clans_file,
                               structures_file, tmpdir)
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

    logger.info("waiting for tasks")
    running = len(workers)
    n_done = n_documents = 0
    step = milestone = 1
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
    logger.info("creating additional entry-structure documents")
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
                    entries[entry_acc].update({
                        "set_acc": clan_acc.lower(),
                        "set_db": entry.database,
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

    logger.info("creating additional entry documents")
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
                "set_db": entry.database,
                "text_set": join(clan_acc, clan_name),
            })

        documents.append((
            get_rel_doc_index(doc),
            get_rel_doc_id(doc),
            doc
        ))

    # TODO: Do we need these? Find case and decide
    logger.info("creating additional taxon documents")
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

    fd, file = mkstemp(dir=tmpdir)
    with open(fd, "wb") as fh:
        pickle.dump(documents, fh)

    n_documents += len(documents)
    move_docs_file(file, directories)

    shutil.rmtree(tmpdir)

    for path in outdirs:
        open(os.path.join(path, f"{version}{config.DONE_SUFFIX}"), "w").close()

    logger.info(f"done ({n_documents:,} documents)")


def move_docs_file(src: str, directories: list[Directory]):
    for dir_obj in directories:
        dst = dir_obj.mktemp(createfile=False)
        shutil.copy(src, dst)
        os.rename(dst, f"{dst}{config.EXTENSION}")

    os.unlink(src)


def export_documents(proteins_file: str, matches_file: str, domorgs_file: str,
                     protein2proteome_file: str, uniprot2pdb_file: str,
                     pdbmatches_file: str, alphafold_file: str,
                     proteomes_file: str, structures_file: str, clans_file: str,
                     entries_file: str, taxa_file: str, outdirs: list[str],
                     version: str, cachesize: int = 100000):
    directories = []
    for path in outdirs:
        try:
            shutil.rmtree(path)
        except FileNotFoundError:
            pass

        os.makedirs(path, mode=0o775)
        directories.append(Directory(root=path))
        open(os.path.join(path, f"{version}{config.LOAD_SUFFIX}"), "w").close()

    logger.info("loading PDBe data")
    with open(structures_file, "rb") as fh:
        # memory: 4.2 GB
        structures = pickle.load(fh)

    with open(uniprot2pdb_file, "rb") as fh:
        # memory: 0.5 GB
        uniprot2pdb = pickle.load(fh)

    # memory: 9.1 GB
    pdb2entry = {}
    pdb2seqlen = {}
    with shelve.open(pdbmatches_file, writeback=False) as d:
        for pdb_chain, pdb_entry in d.items():
            pdb2entry[pdb_chain] = {}
            pdb2seqlen[pdb_chain] = pdb_entry["length"]
            for entry_acc, match in pdb_entry["matches"].items():
                pdb2entry[pdb_chain][entry_acc] = match["locations"]

    logger.info("loading proteomes")
    with open(proteomes_file, "rb") as fh:
        # memory: 0.5 GB
        proteomes = pickle.load(fh)

    logger.info("loading taxonomy")
    with open(taxa_file, "rb") as fh:
        # memory: 2.6 GB
        taxa = pickle.load(fh)

    logger.info("loading entries and clans")
    with open(entries_file, "rb") as fh:
        # memory: 2.1 GB
        entries = pickle.load(fh)

    # memory: 0.5 GB
    member2clan = {}
    with open(clans_file, "rb") as fh:
        for clan in pickle.load(fh).values():
            for entry_acc, _, _, _, _ in clan["members"]:
                member2clan[entry_acc] = (clan["accession"], clan["name"])

    logger.info("writing protein-based documents")
    proteins_store = KVStore(proteins_file)
    matches_store = KVStore(matches_file)
    proteomes_store = KVStore(protein2proteome_file)
    alphafold_store = KVStore(alphafold_file)
    domorgs_store = KVStore(domorgs_file)

    i = 0
    documents = []
    num_documents = 0
    seen_domains = set()
    seen_entries = set()
    seen_structures = set()
    seen_taxa = set()
    for i, (protein_acc, protein) in enumerate(proteins_store.items()):
        taxon_id = protein["taxid"]
        taxon = taxa[taxon_id]
        seen_taxa.add(taxon_id)

        try:
            domain = domorgs_store[protein_acc]
        except KeyError:
            domain_id = domain_str = None
            domain_members = set()
        else:
            domain_id = domain["id"]
            domain_str = domain["key"]
            domain_members = domain["members"]
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

        af_models = alphafold_store.get(protein_acc, [])
        if af_models:
            # af_models: sorted list of tuples (AFDB ID, pLDDT score)
            # sorted by pLDDT score (ascending order)
            af_score = af_models[-1][1]
        else:
            af_score = -1

        # Creates an empty document (all properties set to None)
        doc = init_rel_doc()
        doc.update({
            "protein_acc": protein_acc.lower(),
            "protein_length": protein["length"],
            "protein_is_fragment": protein["fragment"],
            "protein_af_score": af_score,
            "protein_db": "reviewed" if protein["reviewed"] else "unreviewed",
            "text_protein": join(protein_acc,
                                 protein["identifier"],
                                 taxon["sci_name"]),

            # Taxonomy
            "tax_id": taxon_id,
            "tax_name": taxon["sci_name"],
            "tax_lineage": taxon["lineage"],
            "tax_rank": taxon["rank"],
            "text_taxonomy": join(taxon_id, taxon["full_name"], taxon["rank"])
        })

        proteome_id = proteomes_store.get(protein_acc)
        if proteome_id:
            # Adds proteome
            proteome = proteomes[proteome_id]
            doc.update({
                "proteome_acc": proteome_id.lower(),
                "proteome_name": proteome["name"],
                "proteome_is_reference": proteome["is_reference"],
                "text_proteome": join(proteome_id,
                                      proteome["name"],
                                      proteome["assembly"],
                                      proteome["taxon_id"],
                                      proteome["strain"]),
            })

        # PDBe structures and chains (without entries)
        pdb_documents = {}
        # Entries matching structure chains associated to the protein
        pdb_entries = {}
        for pdb_chain, segments in uniprot2pdb.get(protein_acc, {}).items():
            pdb_id, chain_id = pdb_chain.split("_")

            try:
                structure = structures[pdb_id]
            except KeyError:
                continue

            try:
                chain_seq_length = pdb2seqlen[pdb_chain]
            except KeyError:
                continue

            locations = []
            for segment in segments:
                locations.append({
                    "fragments": [{
                        # Coordinates of UniProt entry on the PDB sequence
                        "start": segment["structure_start"],
                        "end": segment["structure_end"],
                        # Coordinates of PDB entry on the UniProt sequence
                        "protein_start": segment["protein_start"],
                        "protein_end": segment["protein_end"],
                    }]
                })

            pdb_doc = deepcopy(doc)
            pdb_doc.update({
                "structure_acc": pdb_id.lower(),
                "structure_resolution": structure["resolution"],
                "structure_date": structure["date"],
                "structure_evidence": structure["evidence"],
                "text_structure": join(pdb_id,
                                       structure["evidence"],
                                       structure["name"]),
                "structure_chain_acc": chain_id,
                "structure_chain": pdb_chain,
                "structure_protein_acc": protein_acc.lower(),
                "structure_protein_length": chain_seq_length,
                "structure_protein_locations": locations,
            })
            pdb_documents[pdb_chain] = pdb_doc

            # Adds entries matching the sequence of this PDBe chain
            for entry_acc in pdb2entry.get(pdb_chain, {}):
                try:
                    pdb_entries[entry_acc].append(pdb_chain)
                except KeyError:
                    pdb_entries[entry_acc] = [pdb_chain]

        # Entries matching the protein (UniProt) sequence
        s_matches, e_matches = matches_store.get(protein_acc, ({}, {}))
        protein_entries = {**s_matches, **e_matches}

        # structures/chains associated to at least one entry
        structures_with_entries = set()
        # number of documents for this protein
        protein_docs = 0

        all_entries = set(protein_entries.keys()) | set(pdb_entries.keys())
        for entry_acc in all_entries:
            seen_entries.add(entry_acc)
            entry = entries[entry_acc]
            entry_database = entry.database.lower()

            if entry.integrated_in:
                integrated_in = entry.integrated_in.lower()
            else:
                integrated_in = None

            entry_doc = deepcopy(doc)
            entry_doc.update({
                "entry_acc": entry_acc.lower(),
                "entry_db": entry_database,
                "entry_type": entry.type.lower(),
                "entry_date": entry.creation_date.strftime("%Y-%m-%d"),
                "entry_go_terms": [t["identifier"] for t in
                                   entry.go_terms],
                "entry_integrated": integrated_in,
                "text_entry": join(entry_acc, entry.short_name, entry.name,
                                   entry.type.lower(), integrated_in),
            })

            if entry_acc in member2clan:
                clan_acc, clan_name = member2clan[entry_acc]
                entry_doc.update({
                    "set_acc": clan_acc.lower(),
                    "set_db": entry_database,
                    "text_set": join(clan_acc, clan_name),
                })

            if entry_acc in domain_members:
                entry_doc.update({
                    "ida_id": domain_id,
                    "ida": domain_str,
                })

            if entry_acc in protein_entries:
                locations = protein_entries[entry_acc]["locations"]

                if entry_database == "panther":
                    """
                    PANTHER: remove the node ID 
                    (other databases do not have a subfamily property)
                    """
                    for loc in locations:
                        try:
                            del loc["subfamily"]["node"]
                        except KeyError:
                            continue  # No subfamily annotation

                entry_doc["entry_protein_locations"] = locations
            else:
                entry_doc["protein_acc"] = None
                entry_doc["protein_is_fragment"] = None
                entry_doc["protein_af_score"] = None
                entry_doc["text_protein"] = None

            if entry_acc in pdb_entries:
                for pdb_chain in pdb_entries[entry_acc]:
                    structures_with_entries.add(pdb_chain)
                    locations = pdb2entry[pdb_chain][entry_acc]
                    _entry_doc = deepcopy(entry_doc)

                    for k, v in pdb_documents[pdb_chain].items():
                        if "structure" in k and v is not None:
                            _entry_doc[k] = v

                    _entry_doc["entry_structure_locations"] = locations
                    documents.append((
                        entry_database + version,
                        get_rel_doc_id(_entry_doc),
                        _entry_doc
                    ))
                    protein_docs += 1
            else:
                documents.append((
                    entry_database + version,
                    get_rel_doc_id(entry_doc),
                    entry_doc
                ))
                protein_docs += 1

        for pdb_chain, pdb_doc in pdb_documents.items():
            seen_structures.add(pdb_chain)
            if pdb_chain in structures_with_entries:
                # Already associated to an entry (therefore to the protein)
                continue

            # TODO: should we do this? Find case and decide.
            pdb_doc.update({
                "ida_id": domain_id,
                "ida": domain_str
            })

            documents.append((
                # Not overlapping any entry -> not associated to a member DB
                config.REL_DEFAULT_INDEX + version,
                get_rel_doc_id(pdb_doc),
                pdb_doc
            ))
            protein_docs += 1

        if protein_docs == 0:
            # No relationships for this protein: fallback to protein doc
            documents.append((
                config.REL_DEFAULT_INDEX + version,
                get_rel_doc_id(doc),
                doc
            ))

        while len(documents) >= cachesize:
            for dir_obj in directories:
                filepath = dir_obj.mktemp()

                with open(filepath, "wb") as fh:
                    pickle.dump(documents[:cachesize], fh)

                os.rename(filepath, f"{filepath}{config.EXTENSION}")

            del documents[:cachesize]
            num_documents += cachesize

        if (i + 1) % 1e7 == 0:
            logger.info(f"{i + 1:>15,}")

    proteins_store.close()
    matches_store.close()
    proteomes_store.close()
    alphafold_store.close()
    domorgs_store.close()
    logger.info(f"{i + 1:>15,}")

    """
    Add additional entry-structure pairs
    i.e. entries with PDB matches where the PDB structure is not associated
    to a protein in UniProtKB
    """
    logger.info("writing entry-structure documents")
    for pdb_chain, structure_entries in pdb2entry.items():
        if pdb_chain in seen_structures:
            continue

        pdb_id, chain_id = pdb_chain.split("_")

        try:
            structure = structures[pdb_id]
        except KeyError:
            continue

        for entry_acc, locations in structure_entries.items():
            entry = entries[entry_acc]

            if entry.deletion_date or not entry.public:
                continue

            database = entry.database.lower()
            if entry.integrated_in:
                integrated_in = entry.integrated_in.lower()
            else:
                integrated_in = None

            doc = init_rel_doc()
            doc.update({
                "entry_acc": entry.accession.lower(),
                "entry_db": database,
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
                "structure_protein_length": pdb2seqlen[pdb_chain],
                "entry_structure_locations": locations
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

            if entry_acc in member2clan:
                clan_acc, clan_name = member2clan[entry.accession]
                doc.update({
                    "set_acc": clan_acc.lower(),
                    "set_db": database,
                    "text_set": join(clan_acc, clan_name),
                })

            documents.append((
                database + version,
                get_rel_doc_id(doc),
                doc
            ))
            seen_entries.add(entry_acc)

    # Adds unseen entries
    logger.info("writing entry documents")
    for entry in entries.values():
        if (entry.accession in seen_entries or
                entry.deletion_date or
                not entry.public):
            continue

        if entry.integrated_in:
            integrated_in = entry.integrated_in.lower()
        else:
            integrated_in = None

        database = entry.database.lower()
        doc = init_rel_doc()
        doc.update({
            "entry_acc": entry.accession.lower(),
            "entry_db": database,
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
                "set_db": database,
                "text_set": join(clan_acc, clan_name),
            })

        documents.append((
            database + version,
            get_rel_doc_id(doc),
            doc
        ))

    # TODO: Do we need these? Find case and decide
    # Adds unseen taxa
    logger.info("writing taxon documents")
    for taxon in taxa.values():
        if taxon["id"] in seen_taxa:
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
            config.REL_DEFAULT_INDEX + version,
            get_rel_doc_id(doc),
            doc
        ))

    num_documents += len(documents)
    while documents:
        for dir_obj in directories:
            filepath = dir_obj.mktemp()

            with open(filepath, "wb") as fh:
                pickle.dump(documents[:cachesize], fh)

            os.rename(filepath, f"{filepath}{config.EXTENSION}")

        del documents[:cachesize]

    # Adds done sentinel file
    for path in outdirs:
        open(os.path.join(path, f"{version}{config.DONE_SUFFIX}"), "w").close()

    logger.info(f"complete: {num_documents:,} documents")
