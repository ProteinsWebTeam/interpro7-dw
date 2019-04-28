import hashlib
import time
from multiprocessing import Process, Queue
from tempfile import mkdtemp
from typing import List

from . import index, set_ready
from .. import mysql
from ... import logger, pdbe
from ...io import JsonFileOrganizer, Store


NODB_INDEX = "others"


class DocumentProducer(Process):
    def __init__(self, ora_ipr: str, my_ipr: str, task_queue: Queue,
                 done_queue: Queue, outdir: str, min_overlap: int=20,
                 docs_per_file: int=10000):
        super().__init__()
        self.ora_ipr = ora_ipr
        self.my_ipr = my_ipr
        self.task_queue = task_queue
        self.done_queue = done_queue
        self.min_overlap = min_overlap
        self.docs_per_file = docs_per_file
        self.organizer = JsonFileOrganizer(outdir, docs_per_file)

        self.entries = {}
        self.integrated = {}
        self.entry2set = {}
        self.proteomes = {}
        self.pfam = set()
        self.structures = {}
        self.protein2pdb = {}

    def run(self):
        # Get PDBe structures, entries, sets, and proteomes
        self.structures = pdbe.get_structures(self.ora_ipr)

        for pdb_id, s in self.structures.items():
            for acc in s["proteins"]:
                if acc in self.protein2pdb:
                    self.protein2pdb[acc].add(pdb_id)
                else:
                    self.protein2pdb[acc] = {pdb_id}

        self.entries = mysql.entry.get_entries(self.my_ipr)
        self.integrated = {
            acc: e["integrated"]
            for acc, e in self.entries.items()
            if e["integrated"]
        }
        self.entry2set = {
            entry_ac: (set_ac, s["database"])
            for set_ac, s in mysql.entry.get_sets(self.my_ipr).items()
            for entry_ac in s["members"]
        }
        self.proteomes = mysql.proteome.get_proteomes(self.my_ipr)

        # List Pfam entries (for IDA)
        self.pfam = {
            e["accession"]
            for e in self.entries.values()
            if e["database"] == "pfam"
        }

        cnt = 0
        types = {
            "protein": self.process_protein,
            "entry": self.process_entry,
            "taxonomy": self.process_taxonomy
        }

        for _type, chunk in iter(self.task_queue.get, None):
            fn = types[_type]

            for args in chunk:
                for doc in fn(*args):
                    self.organizer.add(doc)
                    cnt += 1

        self.organizer.flush()
        self.done_queue.put(cnt)

    @staticmethod
    def repr_frag(f: dict) -> tuple:
        return f["start"], f["end"]

    def process_protein(self, accession: str, identifier: str, name: str,
                        database: str, length: int, comments: list,
                        matches: list, proteome_id: str, taxon: dict) -> list:
        entry_matches = {}
        condensed_entries = {}
        dom_arch = []
        dom_entries = set()
        for m in matches:
            method_ac = m["method_ac"]
            if method_ac in entry_matches:
                e = entry_matches[method_ac]
            else:
                e = entry_matches[method_ac] = []

            e.append({
                "fragments": m["fragments"],
                "model_acc": m["model_ac"],
                "seq_feature": m["seq_feature"]
            })

            entry_ac = self.integrated.get(method_ac)
            if method_ac in self.pfam:
                dom_entries.add(method_ac)
                if entry_ac:
                    dom_entries.add(entry_ac)
                    dom_arch.append("{}:{}".format(method_ac, entry_ac))
                else:
                    dom_arch.append("{}".format(method_ac))

            if entry_ac:
                if entry_ac in condensed_entries:
                    condensed_entries[entry_ac] += m["fragments"]
                else:
                    condensed_entries[entry_ac] = m["fragments"]

        for entry_ac, frags in condensed_entries.items():
            fragments = []
            start = end = None
            for f in sorted(frags, key=self.repr_frag):
                if start is None:
                    # Leftmost fragment
                    start = f["start"]
                    end = f["end"]
                elif f["start"] > end:
                    """
                            end
                        ----]
                              [----
                              s
                    -> new fragment

                       but if end=34 and s=35, we do not want to merge:

                           end s
                          ----][----
                    """
                    fragments.append((start, end))
                    start = f["start"]
                    end = f["end"]
                elif f["end"] > end:
                    """
                            end
                        ----]
                          ------]
                                e
                    -> extend
                    """
                    end = f["end"]

            fragments.append((start, end))
            e = entry_matches[entry_ac] = []
            for start, end in fragments:
                e.append({
                    "fragments": [{
                        "start": start,
                        "end": end,
                        "dc-status": "CONTINUOUS"
                    }],
                    "seq_feature": None,
                    "model_acc": None
                })

        if dom_arch:
            dom_arch = '-'.join(dom_arch)
            dom_arch_id = hashlib.sha1(dom_arch.encode("utf-8")).hexdigest()
        else:
            dom_arch = dom_arch_id = None

        if length <= 100:
            size = "small"
        elif length <= 1000:
            size = "medium"
        else:
            size = "large"

        doc = self.init_document()

        # Add protein info
        doc.update({
            "protein_acc": accession.lower(),
            "protein_length": length,
            "protein_size": size,
            "protein_db": database,
            "text_protein": self._join(
                accession, identifier, name, database, comments
            ),

            "tax_id": taxon["taxId"],
            "tax_name": taxon["scientificName"],
            "tax_lineage": taxon["lineage"].strip().split(),
            "tax_rank": taxon["rank"],
            "text_taxonomy": self._join(
                taxon["taxId"], taxon["fullName"], taxon["rank"]
            ),
        })

        if proteome_id:
            p = self.proteomes[proteome_id]
            doc.update({
                "proteome_acc": proteome_id.lower(),
                "proteome_name": p["name"],
                "proteome_is_reference": p["is_reference"],
                "text_proteome": self._join(proteome_id, *list(p.values()))
            })

        documents = []
        if entry_matches:
            # Add entries
            for entry_ac in entry_matches:
                entry = self.entries[entry_ac]
                if entry["integrated"]:
                    entry["integrated"] = entry["integrated"].lower()

                go_terms = [t["identifier"] for t in entry["go_terms"]]
                _doc = doc.copy()
                _doc.update({
                    "entry_acc": entry["accession"].lower(),
                    "entry_db": entry["database"],
                    "entry_type": entry["type"],
                    "entry_date": entry["date"].strftime("%Y-%m-%d"),
                    "entry_integrated": entry["integrated"],
                    "text_entry": self._join(
                        entry["accession"], entry["name"],
                        entry["type"], entry["descriptions"], *go_terms
                    ),
                    "entry_protein_locations": entry_matches[entry_ac],
                    "entry_go_terms": go_terms
                })

                if entry["accession"] in dom_entries:
                    _doc.update({
                        "ida_id": dom_arch_id,
                        "ida": dom_arch
                    })

                _set = self.entry2set.get(entry_ac)
                if _set:
                    set_ac, set_db = _set
                    _doc.update({
                        "set_acc": set_ac.lower(),
                        "set_db": set_db,
                        # todo: implement set integration (e.g. pathways)
                        "set_integrated": [],
                        "text_set": self._join(set_ac, set_db)
                    })

                documents.append(_doc)
        else:
            documents.append(doc)

        # Add PDBe structures (and chains)
        _documents = []
        for pdbe_id in self.protein2pdb.get(accession, []):
            structure = self.structures[pdbe_id]
            text = self._join(
                pdbe_id,
                structure["evidence"],
                structure["name"],
                ' '.join([
                    pub["title"]
                    for pub in structure["citations"].values()
                    if pub.get("title")
                ])
            )

            for doc in documents:
                chains = structure["proteins"][accession]

                _doc = doc.copy()
                _doc.update({
                    #"structure_acc": pdbe_id.lower(),
                    #"structure_acc": pdbe_id,
                    "structure_resolution": structure["resolution"],
                    "structure_date": structure["date"].strftime("%Y-%m-%d"),
                    "structure_evidence": structure["evidence"],
                    "protein_structure": chains
                })

                for chain_id, chain in chains.items():
                    _doc_chain = _doc.copy()
                    _doc_chain.update({
                        "structure_chain_acc": chain_id,
                        "structure_protein_locations": [
                            {
                                "fragments": [
                                    {
                                        "start": fragment["protein_start"],
                                        "end": fragment["protein_end"]
                                    }
                                ]
                            } for fragment in chain
                        ],
                        "structure_chain": "{} - {}".format(
                            pdbe_id, chain_id
                        ),
                        "text_structure": "{} {}".format(chain_id, text)
                    })
                    _documents.append(_doc_chain)

        if _documents:
            documents = _documents

        for doc in documents:
            doc["id"] = self._join(
                doc["protein_acc"], doc["proteome_acc"], doc["entry_acc"],
                doc["set_acc"], doc["structure_acc"],
                doc["structure_chain_acc"],
                separator='-'
            )

        return documents

    def process_entry(self, accession: str) -> list:
        entry = self.entries[accession]
        if entry["integrated"]:
            entry["integrated"] = entry["integrated"].lower()

        go_terms = [t["identifier"] for t in entry["go_terms"]]
        doc = self.init_document()
        doc.update({
            "entry_acc": entry["accession"].lower(),
            "entry_db": entry["database"],
            "entry_type": entry["type"],
            "entry_date": entry["date"].strftime("%Y-%m-%d"),
            "entry_integrated": entry["integrated"],
            "text_entry": self._join(
                entry["accession"], entry["name"],
                entry["type"], entry["descriptions"], *go_terms
            ),
            "entry_protein_locations": [],
            "entry_go_terms": go_terms
        })

        _set = self.entry2set.get(accession)
        if _set:
            set_ac, set_db = _set
            doc.update({
                "set_acc": set_ac.lower(),
                "set_db": set_db,
                # todo: implement set integration (e.g. pathways)
                "set_integrated": [],
                "text_set": self._join(set_ac, set_db)
            })
            doc["id"] = doc["entry_acc"] + '-' + doc["set_acc"]
        else:
            doc["id"] = doc["entry_acc"]

        return [doc]

    def process_taxonomy(self, taxon: dict) -> list:
        doc = self.init_document()
        doc.update({
            "tax_id": taxon["taxId"],
            "tax_name": taxon["scientificName"],
            "tax_lineage": taxon["lineage"].strip().split(),
            "tax_rank": taxon["rank"],
            "text_taxonomy": self._join(
                taxon["taxId"], taxon["fullName"], taxon["rank"]
            ),
            "id": taxon["taxId"]
        })
        return [doc]

    @staticmethod
    def init_document() -> dict:
        return {
            "id": None,

            # Protein
            "protein_acc": None,
            "protein_length": None,
            "protein_size": None,
            "protein_db": None,
            "text_protein": None,

            # Taxonomy
            "tax_id": None,
            "tax_name": None,
            "tax_lineage": None,
            "tax_rank": None,
            "text_taxonomy": None,

            # Proteome
            "proteome_acc": None,
            "proteome_name": None,
            "proteome_is_reference": None,
            "text_proteome": None,

            # Entry
            "entry_acc": None,
            "entry_db": None,
            "entry_type": None,
            "entry_date": None,
            "entry_protein_locations": None,
            "entry_go_terms": None,
            "entry_integrated": None,
            "text_entry": None,

            # Set
            "set_acc": None,
            "set_db": None,
            "set_integrated": None,
            "text_set": None,

            # Structure
            "structure_acc": None,
            "structure_resolution": None,
            "structure_date": None,
            "structure_evidence": None,

            # Chain
            "structure_chain_acc": None,
            "structure_protein_locations": None,
            "protein_structure_chain": None,
            "structure_chain": None,
            "text_structure": None,

            # Domain architecture
            "ida_id": None,
            "ida": None
        }

    @staticmethod
    def _join(*args, separator: str=' ') -> str:
        # Use underscore to NOT override Process.join()
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

        return separator.join(items).lower()


def create_documents(ora_ipr: str, my_ipr: str, src_proteins: str,
                     src_names: str, src_comments: str, src_proteomes: str,
                     src_matches: str, outdir: str, **kwargs):
    processes = kwargs.get("processes", 1)
    chunk_size = kwargs.get("chunk_size", 100000)
    limit = kwargs.get("limit", 0)

    processes = max(1, processes-1)  # minus one for parent process
    task_queue = Queue(processes)
    done_queue = Queue()
    workers = []
    for _ in range(processes):
        p = DocumentProducer(ora_ipr, my_ipr, task_queue, done_queue,
                             mkdtemp(dir=outdir))
        p.start()
        workers.append(p)

    # MySQL data
    logger.info("loading data from MySQL")
    taxa = mysql.taxonomy.get_taxa(my_ipr, lineage=True)
    integrated = {}
    entry_accessions = set()
    for entry_ac, e in mysql.entry.get_entries(my_ipr).items():
        entry_accessions.add(entry_ac)
        if e["integrated"]:
            integrated[entry_ac] = e["integrated"]

    # Open stores
    proteins = Store(src_proteins)
    protein2names = Store(src_names)
    protein2comments = Store(src_comments)
    protein2proteome = Store(src_proteomes)
    protein2matches = Store(src_matches)

    tax_ids = set(taxa.keys())

    logger.info("starting")
    n_proteins = 0
    chunk = []
    entries_with_matches = set()
    enqueue_time = 0
    ts = time.time()

    for acc, protein in proteins:
        tax_id = protein["taxon"]
        taxon = taxa[tax_id]

        name, other_names = protein2names.get(acc, (None, None))
        matches = protein2matches.get(acc, [])

        # Enqueue protein
        chunk.append((
            acc,
            protein["identifier"],
            name,
            "reviewed" if protein["isReviewed"] else "unreviewed",
            protein["length"],
            protein2comments.get(acc, []),
            matches,
            protein2proteome.get(acc),
            taxon
        ))

        if len(chunk) == chunk_size:
            t = time.time()
            task_queue.put(("protein", chunk))
            enqueue_time += time.time() - t
            chunk = []

        # Keep track of taxa associated to at least one protein
        try:
            tax_ids.remove(tax_id)
        except KeyError:
            pass

        # Keep track of entries with protein matches
        for m in matches:
            method_ac = m["method_ac"]
            entries_with_matches.add(method_ac)

            if method_ac in integrated:
                entries_with_matches.add(integrated[method_ac])

        n_proteins += 1
        if n_proteins == limit:
            break
        elif not n_proteins % 10000000:
            logger.info("{:>12,} ({:.0f} proteins/sec)".format(
                n_proteins, n_proteins / (time.time() - ts)
            ))

    if chunk:
        t = time.time()
        task_queue.put(("protein", chunk))
        enqueue_time += time.time() - t

    logger.info("{:>12,} ({:.0f} proteins/sec)".format(
        n_proteins, n_proteins / (time.time() - ts)
    ))

    # Add entries without matches
    chunk = [
        (entry_ac,)
        for entry_ac in entry_accessions - entries_with_matches
    ]
    for i in range(0, len(chunk), chunk_size):
        t = time.time()
        task_queue.put(("entry", chunk[i:i+chunk_size]))
        enqueue_time += time.time() - t

    # Add taxa without proteins
    chunk = [(taxa[tax_id],) for tax_id in tax_ids]
    for i in range(0, len(chunk), chunk_size):
        t = time.time()
        task_queue.put(("taxonomy", chunk[i:i+chunk_size]))
        enqueue_time += time.time() - t

    logger.info("enqueue time: {:>10.0f} seconds".format(enqueue_time))

    # Poison pill
    for _ in workers:
        task_queue.put(None)

    # Closing stores
    proteins.close()
    protein2names.close()
    protein2comments.close()
    protein2proteome.close()
    protein2matches.close()

    n_docs = sum([done_queue.get() for _ in workers])

    # Wait for workers to finish
    for p in workers:
        p.join()

    # Delete loading file so Loaders know that all files are generated
    set_ready(outdir)

    logger.info("complete: {:,} documents".format(n_docs))


class DocumentController(index.DocumentController):
    def __init__(self, **kwargs):
        self.suffix = kwargs.get("suffix", "")

    def dump(self, doc: dict) -> dict:
        if doc["entry_db"]:
            idx = doc["entry_db"] + self.suffix
        else:
            idx = NODB_INDEX + self.suffix

        return {
            "_op_type": "index",
            "_index": idx,
            "_type": "relationship",
            "_id": doc["id"],
            "_source": doc
        }

    def parse(self, item: dict):
        try:
            del item["index"]["data"]
        except KeyError:
            pass


def index_documents(my_ipr: str, hosts: List[str], src: str, **kwargs):
    # Load databases (base of indices)
    indices = list(mysql.database.get_databases(my_ipr).keys())
    indices.append(NODB_INDEX)

    if kwargs.get("body_path"):
        # Create indices
        index.create_indices(hosts=hosts,
                             indices=indices,
                             body_path=kwargs.pop("body_path"),
                             doc_type="relationship",
                             **kwargs
                             )

    controller = DocumentController(**kwargs)
    alias = kwargs.get("alias")
    if index.index_documents(hosts, controller, src, **kwargs) and alias:
        index.update_alias(hosts, indices, alias, **kwargs)


def update_alias(my_ipr: str, hosts: List[str], alias: str, **kwargs):
    indices = list(mysql.database.get_databases(my_ipr).keys())
    indices.append(NODB_INDEX)
    index.update_alias(hosts, indices, alias, **kwargs)
