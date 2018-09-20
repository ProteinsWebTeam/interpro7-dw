#!/usr/bin/env python
# -*- coding: utf-8 -*-

import glob
import gzip
import hashlib
import json
import logging
import os
import shutil
import tempfile
import time
from multiprocessing import Process, Queue

from elasticsearch import Elasticsearch, helpers, exceptions

from . import dbms, disk, mysql
from .ebi import interpro, pdbe

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s: %(levelname)s: %(message)s',
    datefmt='%y-%m-%d %H:%M:%S'
)

LOADING_FILE = "loading"
EXTRA_INDEX = "others"


def init_dir(path):
    if os.path.isdir(path):
        shutil.rmtree(path)

    os.makedirs(path)
    open(os.path.join(path, LOADING_FILE), "w").close()


def parse_host(host: str) -> dict:
    host = host.split(':')
    if len(host) == 2:
        return {
            "host": host[0],
            "port": int(host[1])
        }
    else:
        return {
            "host": host,
            "port": 9200
        }


class DocumentProducer(Process):
    def __init__(self, ora_ippro: str, my_ippro: str, queue_in: Queue,
                 queue_out: Queue, outdir: str, **kwargs):
        super().__init__()
        self.ora_ippro = ora_ippro
        self.my_ippro = my_ippro
        self.queue_in = queue_in
        self.queue_out = queue_out
        self.outdir = outdir
        self.min_overlap = kwargs.get("min_overlap", 20)
        self.max_size = kwargs.get("max_size", 1000000)
        self.chunk_size = kwargs.get("chunk_size", 10000)

        self.entries = None
        self.sets = None
        self.proteomes = None
        self.pfam = set()
        self.structures = None

    def run(self):
        logging.info("{} ({}) started".format(self.name, os.getpid()))

        # Get PDBe structures, entries, sets, and proteomes
        self.structures = pdbe.get_structures(self.ora_ippro,
                                              citations=True,
                                              fragments=True,
                                              by_protein=True)
        self.entries = mysql.get_entries(self.my_ippro)
        self.sets = mysql.get_sets(self.my_ippro)
        self.proteomes = mysql.get_proteomes(self.my_ippro)

        # List Pfam entries (for IDA)
        self.pfam = {
            e["accession"]
            for e in self.entries.values()
            if e["database"] == "pfam"
        }

        documents = []
        cnt = 0
        types = {
            "protein": self.process_protein,
            "entry": self.process_entry,
            "taxonomy": self.process_taxonomy
        }

        while True:
            task = self.queue_in.get()
            if task is None:
                break

            _type, chunk = task
            fn = types[_type]

            for args in chunk:
                documents += fn(*args)

                if len(documents) >= self.max_size:
                    cnt += len(documents)
                    self.dump(documents, self.outdir, self.chunk_size)
                    documents = []

        if documents:
            cnt += len(documents)
            self.dump(documents, self.outdir, self.chunk_size)
            documents = []

        logging.info("{} ({}) terminated ({} documents)".format(
            self.name, os.getpid(), cnt)
        )

    def process_protein(self, accession: str, identifier: str, name: str,
                        database: str, length: int, comments: list,
                        matches: list, proteomes: list, taxon: dict) -> list:
        # Prepare matches/supermatches
        entry_matches = {}
        supermatches = []
        dom_arch = []
        dom_entries = set()
        for m in matches:
            entry_ac = m["entry_ac"]
            method_ac = m["method_ac"]
            if m["model_ac"] and m["model_ac"] != method_ac:
                model_ac = m["model_ac"]
            else:
                model_ac = None

            # todo: remove when I5 bug fixed
            fragments = [f for f in m["fragments"] if f["start"] <= f["end"]]

            if method_ac in entry_matches:
                e = entry_matches[method_ac]
            else:
                e = entry_matches[method_ac] = []

            e.append({
                "fragments": fragments,
                "model_acc": model_ac,
                "seq_feature": m["seq_feature"]
            })

            if method_ac in self.pfam:
                dom_entries.add(method_ac)
                if entry_ac:
                    dom_entries.add(entry_ac)
                    dom_arch.append("{}:{}".format(method_ac, entry_ac))
                else:
                    dom_arch.append("{}".format(method_ac))

            if entry_ac:
                try:
                    entry = self.entries[entry_ac]
                except KeyError:
                    continue  # todo: log error

                supermatches.append(
                    interpro.Supermatch(
                        entry_ac,
                        entry["root"],
                        min([f["start"] for f in fragments]),
                        max([f["end"] for f in fragments])
                    )
                )

        if dom_arch:
            dom_arch = '-'.join(dom_arch)
            dom_arch_id = hashlib.sha1(dom_arch.encode("utf-8")).hexdigest()
        else:
            dom_arch = dom_arch_id = None

        # Merge overlapping supermatches
        prot_supermatches = []
        sets = interpro.merge_supermatches(supermatches, self.min_overlap)
        for s in sets:
            for sm in s.supermatches:
                for entry_ac in sm.get_entries():
                    # Supermatch rows
                    prot_supermatches.append((
                        entry_ac.upper(),
                        sm.start,
                        sm.end
                    ))

                    # Add supermatches to Elastic matches
                    if entry_ac in entry_matches:
                        e = entry_matches[entry_ac]
                    else:
                        e = entry_matches[entry_ac] = []

                    e.append({
                        "fragments": [{'start': sm.start, 'end': sm.end}],
                        "model_acc": None,
                        "seq_feature": None
                    })

        self.queue_out.put((
            accession,
            'S' if database == "reviewed" else 'T',
            prot_supermatches
        ))

        if length <= 100:
            size = "small"
        elif length <= 1000:
            size = "medium"
        else:
            size = "large"

        doc = self.init_document()

        # Accession in lower case
        accession_lc = accession.lower()

        # Add protein info
        doc.update({
            "protein_acc": accession_lc,
            "protein_length": length,
            "protein_size": size,
            "protein_db": database,
            "text_protein": self._join(
                accession_lc, identifier, name, database, comments
            ),

            "tax_id": taxon["taxId"],
            "tax_name": taxon["scientificName"],
            "tax_lineage": taxon["lineage"].strip().split(),
            "tax_rank": taxon["rank"],
            "text_taxonomy": self._join(
                taxon["taxId"], taxon["fullName"], taxon["rank"]
            ),
        })

        # Add proteomes
        documents = [doc]
        _documents = []
        for upid in proteomes:
            try:
                p = self.proteomes[upid]
            except KeyError:
                continue  # todo: log error

            for doc in documents:
                _doc = doc.copy()
                _doc.update({
                    "proteome_acc": upid,
                    "proteome_name": p["name"],
                    "proteome_is_reference": p["is_reference"],
                    "text_proteome": self._join(upid, *list(p.values()))
                })
                _documents.append(_doc)

        if _documents:
            documents = _documents
            _documents = []

        # Add entries
        for entry_ac in entry_matches:
            try:
                entry = self.entries[entry_ac]
            except KeyError:
                continue  # todo: log error

            sets = self.sets.get(entry_ac)
            go_terms = [t["identifier"] for t in entry["go_terms"]]
            for doc in documents:
                _doc = doc.copy()
                _doc.update({
                    "entry_acc": entry["accession"],
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

                if sets:
                    for set_ac, set_db in sets.items():
                        _doc_set = _doc.copy()
                        _doc_set.update({
                            "set_acc": set_ac,
                            "set_db": set_db,
                            # todo: implement set integration (e.g. pathways)
                            "set_integrated": [],
                            "text_set": self._join(set_ac, set_db)
                        })
                        _documents.append(_doc_set)
                else:
                    _documents.append(_doc)

        if _documents:
            documents = _documents
            _documents = []

        # Add PDBe structures (and chains)
        for structure in self.structures.get(accession, {}).values():
            text = self._join(
                structure["accession"],
                structure["evidence"],
                structure["name"],
                ' '.join([
                    pub["title"]
                    for pub in structure["citations"].values()
                    if pub.get("title")
                ])
            )

            for doc in documents:
                _doc = doc.copy()
                _doc.update({
                    "structure_acc": structure["accession"],
                    "structure_resolution": structure["resolution"],
                    "structure_date": structure["date"].strftime("%Y-%m-%d"),
                    "structure_evidence": structure["evidence"]
                })

                for chain in structure["proteins"][accession]:
                    fragments = structure["proteins"][accession][chain]
                    _doc_chain = _doc.copy()
                    _doc_chain.update({
                        "structure_chain_acc": chain,
                        "protein_structure_locations": [
                            {"fragments": [
                                {"start": m["start"], "end": m["end"]}
                            ]} for m in fragments
                        ],
                        "structure_chain": "{} - {}".format(
                            structure["accession"], chain
                        ),
                        "text_structure": "{} {}".format(chain, text)
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
        try:
            entry = self.entries[accession]
        except KeyError:
            return []  # todo: log error

        go_terms = [t["identifier"] for t in entry["go_terms"]]
        doc = self.init_document()
        doc.update({
            "entry_acc": entry["accession"],
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

        sets = self.sets.get(accession)
        if sets:
            documents = []
            for set_ac, set_db in sets.items():
                _doc = doc.copy()
                _doc.update({
                    "set_acc": set_ac,
                    "set_db": set_db,
                    # todo: implement set integration (e.g. pathways)
                    "set_integrated": [],
                    "text_set": self._join(set_ac, set_db),
                    "id": self._join(
                        entry["accession"], set_ac, separator='-'
                    )
                })
                documents.append(_doc)
            return documents
        else:
            doc["id"] = entry["accession"]
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
    def dump(documents: list, outdir: str, chunk_size: int):
        if len(documents) > chunk_size:
            """
            Too many documents for one single file: 
            create a directory and write files inside
            """
            outdir = tempfile.mkdtemp(dir=outdir)

            for i in range(0, len(documents), chunk_size):
                fd, path = tempfile.mkstemp(dir=outdir)
                os.close(fd)

                with gzip.open(path, "wt") as fh:
                    json.dump(documents[i:i+chunk_size], fh)

                os.rename(path, path + ".json.gz")
        else:
            # All documents fit in a single file
            fd, path = tempfile.mkstemp(dir=outdir)
            os.close(fd)

            with gzip.open(path, "wt") as fh:
                json.dump(documents, fh)

            os.rename(path, path + ".json.gz")

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
            "protein_structure_locations": None,
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

        return separator.join(items)


class SupermatchConsumer(Process):
    def __init__(self, my_ippro: str, queue_in: Queue,
                 **kwargs):
        super().__init__()
        self.my_ippro = my_ippro
        self.queue_in = queue_in
        self.ora_ippro = kwargs.get("ora_ippro")
        self.threshold = kwargs.get("threshold", 0.75)
        self.chunk_size = kwargs.get("chunk_size", 10000)
        self.types = ("homologous_superfamily", "domain", "family", "repeat")

    def run(self):
        logging.info("{} ({}) started".format(self.name, os.getpid()))

        if self.ora_ippro:
            con, cur = dbms.connect(self.ora_ippro)

            try:
                cur.execute(
                    """
                    DROP TABLE INTERPRO.SUPERMATCH2
                    CASCADE CONSTRAINTS 
                    """
                )
            except:
                pass

            cur.execute(
                """
                CREATE TABLE INTERPRO.SUPERMATCH2
                (
                    PROTEIN_AC VARCHAR2(15) NOT NULL,
                    DBCODE CHAR(1) NOT NULL, 
                    ENTRY_AC VARCHAR2(9) NOT NULL,
                    POS_FROM NUMBER(5) NOT NULL,
                    POS_TO NUMBER(5) NOT NULL
                ) NOLOGGING            
                """
            )
        else:
            con = cur = None

        insert_query = """
            INSERT /*+APPEND*/ INTO INTERPRO.SUPERMATCH2 (
              PROTEIN_AC, DBCODE, ENTRY_AC, POS_FROM, POS_TO
            )
            VALUES (:1, :2, :3, :4, :5)
        """

        cnt = 0
        rows = []
        sets = {}
        overlaps = {}
        while True:
            task = self.queue_in.get()
            if task is None:
                break

            # Data for one single protein
            accession, dbcode, supermatches = task

            matches = {}
            for entry_ac, start, end in supermatches:
                cnt += 1

                if con is not None:
                    # Insert supermatch in Oracle
                    rows.append((
                        accession,
                        dbcode,
                        entry_ac,
                        start,
                        end
                    ))

                    if len(rows) == self.chunk_size:
                        cur.executemany(insert_query, rows)
                        con.commit()
                        rows = []

                # Current implementation: leftmost match only
                if entry_ac not in matches:
                    matches[entry_ac] = [(start, end)]

            self.intersect(matches, sets, overlaps)

        if con is not None:
            if rows:
                cur.executemany(insert_query, rows)
                con.commit()
                rows = []

            # Constraints
            cur.execute(
                """
                ALTER TABLE INTERPRO.SUPERMATCH2
                ADD CONSTRAINT PK_SUPERMATCH2
                PRIMARY KEY (PROTEIN_AC, ENTRY_AC, POS_FROM, POS_TO, DBCODE)
                """
            )

            try:
                cur.execute(
                    """
                    ALTER TABLE INTERPRO.SUPERMATCH2
                    ADD CONSTRAINT FK_SUPERMATCH2$PROTEIN_AC
                    FOREIGN KEY (PROTEIN_AC) 
                    REFERENCES INTERPRO.PROTEIN (PROTEIN_AC)
                    ON DELETE CASCADE
                    """
                )
            except:
                pass

            try:
                cur.execute(
                    """
                    ALTER TABLE INTERPRO.SUPERMATCH2
                    ADD CONSTRAINT FK_SUPERMATCH2$ENTRY_AC
                    FOREIGN KEY (ENTRY_AC) 
                    REFERENCES INTERPRO.ENTRY (ENTRY_AC)
                    ON DELETE CASCADE
                    """
                )
            except:
                pass

            # Indexes
            cur.execute(
                """
                CREATE INDEX I_SUPERMATCH2$PROTEIN 
                ON INTERPRO.SUPERMATCH2 (PROTEIN_AC) 
                NOLOGGING
                """
            )
            cur.execute(
                """
                CREATE INDEX I_SUPERMATCH2$ENTRY 
                ON INTERPRO.SUPERMATCH2 (ENTRY_AC) 
                NOLOGGING
                """
            )
            cur.execute(
                """
                CREATE INDEX I_SUPERMATCH2$DBCODE$ENTRY 
                ON INTERPRO.SUPERMATCH2 (DBCODE, ENTRY_AC) 
                NOLOGGING
                """
            )

            # Statistics
            cur.execute(
                """
                    BEGIN
                        DBMS_STATS.GATHER_TABLE_STATS(:1, :2, cascade => TRUE);
                    END;
                """,
                ("INTERPRO", "SUPERMATCH2")
            )

            # Privileges
            cur.execute(
                """
                GRANT SELECT 
                ON INTERPRO.SUPERMATCH2 
                TO INTERPRO_SELECT
                """
            )

            cur.close()
            con.close()

        entries = mysql.get_entries(self.my_ippro)
        con, cur = dbms.connect(self.my_ippro)

        # Compute Jaccard coefficients
        overlapping = {}

        for acc1 in overlaps:
            s1 = sets[acc1]

            for acc2 in overlaps[acc1]:
                s2 = sets[acc2]
                o1, o2 = overlaps[acc1][acc2]

                # Independent coefficients
                coef1 = o1 / (s1 + s2 - o1)
                coef2 = o2 / (s1 + s2 - o2)

                # Final coefficient: average of independent coefficients
                coef = (coef1 + coef2) * 0.5

                # Containment indices
                c1 = o1 / s1
                c2 = o2 / s2

                if any(map(lambda x: x >= self.threshold, (coef, c1, c2))):
                    t1 = entries[acc1.lower()]["type"]
                    t2 = entries[acc2.lower()]["type"]

                    if t1 == "homologous_superfamily":
                        if t2 not in self.types:
                            continue
                    elif t2 == "homologous_superfamily":
                        if t1 not in self.types:
                            continue

                    if acc1 in overlapping:
                        overlapping[acc1].append(acc2)
                    else:
                        overlapping[acc1] = [acc2]

                    if acc2 in overlapping:
                        overlapping[acc2].append(acc1)
                    else:
                        overlapping[acc2] = [acc1]

        for acc in overlapping:
            cur.execute(
                """
                UPDATE webfront_entry
                SET overlaps_with = %s
                WHERE accession = %s
                """,
                (json.dumps(overlapping[acc]), acc)
            )

        cur.close()
        con.commit()
        con.close()

        logging.info("{} ({}) terminated ({} supermatches)".format(
            self.name, os.getpid(), cnt)
        )

    @staticmethod
    def intersect(matches: dict, sets: dict, intersections: dict):
        for acc1 in matches:
            if acc1 in sets:
                sets[acc1] += 1
            else:
                sets[acc1] = 1

            for acc2 in matches:
                if acc1 >= acc2:
                    continue
                elif acc1 not in intersections:
                    intersections[acc1] = {acc2: [0, 0]}
                elif acc2 not in intersections[acc1]:
                    intersections[acc1][acc2] = [0, 0]

                m1 = matches[acc1][0]
                m2 = matches[acc2][0]
                o = min(m1[1], m2[1]) - max(m1[0], m2[0]) + 1

                l1 = m1[1] - m1[0] + 1
                l2 = m2[1] - m2[0] + 1

                if o > l1 * 0.5:
                    # acc1 is in acc2 (because it overlaps acc2 at least 50%)
                    intersections[acc1][acc2][0] += 1

                if o > l2 * 0.5:
                    # acc2 is in acc1
                    intersections[acc1][acc2][1] += 1


def create_documents(ora_ippro, my_ippro, proteins_f, descriptions_f,
                     comments_f, proteomes_f, prot_matches_f, outdir,
                     **kwargs):
    n_producers = kwargs.get("producers", 1)
    threshold = kwargs.get("threshold", 0.75)
    chunk_size = kwargs.get("chunk_size", 100000)
    limit = kwargs.get("limit", 0)
    populate_oracle = kwargs.get("populate_oracle", True)

    doc_queue = Queue(n_producers)
    supermatch_queue = Queue()

    consumer = SupermatchConsumer(
        my_ippro, supermatch_queue,
        threshold=threshold,
        ora_ippro=ora_ippro if populate_oracle else None
    )
    consumer.start()

    producers = [
        DocumentProducer(ora_ippro, my_ippro, doc_queue,
                         supermatch_queue, outdir)
        for _ in range(n_producers)
    ]

    for p in producers:
        p.start()

    # MySQL data
    logging.info("loading data from MySQL")
    taxa = mysql.get_taxa(my_ippro, slim=False)
    entries = set(mysql.get_entries(my_ippro))

    # Open stores
    proteins = disk.Store(proteins_f)
    descriptions = disk.Store(descriptions_f)
    comments = disk.Store(comments_f)
    proteomes = disk.Store(proteomes_f)
    prot_matches = disk.Store(prot_matches_f)

    set_taxa = set(taxa.keys())

    logging.info("starting")
    cnt = 0
    total = 0
    chunk = []
    entries_with_matches = set()
    ts = time.time()
    for acc, protein in proteins.iter():
        tax_id = protein["taxon"]
        taxon = taxa[tax_id]

        try:
            set_taxa.remove(tax_id)
        except KeyError:
            pass

        name, other_names = descriptions.get(acc, (None, None))
        matches = prot_matches.get(acc, [])

        # Enqueue protein
        chunk.append((
            acc,
            protein["identifier"],
            name,
            "reviewed" if protein["isReviewed"] else "unreviewed",
            protein["length"],
            comments.get(acc, []),
            matches,
            proteomes.get(acc, []),
            taxon
        ))

        if len(chunk) == chunk_size:
            doc_queue.put(("protein", chunk))
            chunk = []

        # Remove entries with protein matches
        for m in matches:
            entries_with_matches.add(m["method_ac"])
            if m["entry_ac"]:
                entries_with_matches.add(m["entry_ac"])

        total += 1
        cnt += 1
        if total == limit:
            break
        elif not total % 1000000:
            logging.info("{:>12} ({:.0f} proteins/sec)".format(
                total, cnt // (time.time() - ts)
            ))
            cnt = 0
            ts = time.time()

    if chunk:
        doc_queue.put(("protein", chunk))

    logging.info("{:>12} ({:.0f} proteins/sec)".format(
        total, cnt // (time.time() - ts)
    ))

    # Add entries without matches
    chunk = [(entry_ac,) for entry_ac in entries - entries_with_matches]
    for i in range(0, len(chunk), chunk_size):
        doc_queue.put(("entry", chunk[i:i+chunk_size]))

    # Add taxa without proteins
    chunk = [(taxa[tax_id],) for tax_id in set_taxa]
    for i in range(0, len(chunk), chunk_size):
        doc_queue.put(("taxonomy", chunk[i:i+chunk_size]))

    # Poison pill
    for _ in producers:
        doc_queue.put(None)

    # Close stores to free memory
    for store in (proteins, descriptions, comments, proteomes, prot_matches):
        store.close()

    # Wait for producers to finish
    for p in producers:
        p.join()

    # Once producers are done: poisons consumer
    supermatch_queue.put(None)
    consumer.join()

    # Delete loading file so Loaders know that all files are generated
    os.unlink(os.path.join(outdir, LOADING_FILE))

    logging.info("complete")


class DocumentLoader(Process):
    def __init__(self, host, doc_type, queue_in, queue_out, **kwargs):
        super().__init__()
        self.host = host
        self.type = doc_type
        self.queue_in = queue_in
        self.queue_out = queue_out
        self.max_bytes = kwargs.get("max_bytes", 100 * 1024 * 1024)
        self.suffix = kwargs.get("suffix", "")
        self.threads = kwargs.get("threads", 4)

    def run(self):
        logging.info("{} ({}) started".format(self.name, os.getpid()))
        es = Elasticsearch([self.host])

        # Disable Elastic logger
        tracer = logging.getLogger("elasticsearch")
        tracer.setLevel(logging.CRITICAL + 1)

        failed_files = []
        total_documents = 0
        total_errors = 0
        while True:
            filepath = self.queue_in.get()
            if filepath is None:
                break

            with gzip.open(filepath, "rt") as fh:
                documents = json.load(fh)

            actions = []
            for doc in documents:
                # Define which index to use
                if doc["entry_db"]:
                    index = doc["entry_db"] + self.suffix
                else:
                    index = EXTRA_INDEX + self.suffix

                actions.append({
                    '_op_type': 'index',
                    '_index': index,
                    '_type': self.type,
                    '_id': doc["id"],
                    '_source': doc
                })

            gen = helpers.parallel_bulk(
                es, actions,
                thread_count=self.threads,
                queue_size=self.threads,
                # disable chunk_size (num of docs)
                # to only rely on max_chunk_bytes (bytes)
                chunk_size=-1,
                max_chunk_bytes=self.max_bytes,
                raise_on_exception=False,
                raise_on_error=False
            )

            total_documents += len(actions)
            n_errors = len([item for status, item in gen if not status])
            total_errors += n_errors

            if n_errors:
                failed_files.append(filepath)

        self.queue_out.put(failed_files)
        logging.info(
            "{} ({}) terminated "
            "({} documents, {} errors)".format(
                self.name, os.getpid(), total_documents, total_errors
            )
        )


def index_documents(my_ippro: str, host: str, doc_type: str,
                    properties: str, src: str, **kwargs):
    indices = kwargs.get("indices")
    n_loaders = kwargs.get("loaders", 1)
    shards = kwargs.get("shards", 5)
    suffix = kwargs.get("suffix", "").lower()

    # Parse Elastic host (str -> dict)
    host = parse_host(host)

    logging.info("indexing documents to {}".format(host["host"]))

    # Load property mapping
    with open(properties, "rt") as fh:
        properties = json.load(fh)

    if indices:
        # Custom number of shards
        with open(indices, "rt") as fh:
            indices = json.load(fh)

        # Force keys to be in lower case
        indices = {k.lower(): v for k, v in indices.items()}
    else:
        indices = {}

    # Load databases (base of indices)
    databases = list(mysql.get_entry_databases(my_ippro).keys())

    # Establish connection
    es = Elasticsearch([host])

    # Disable Elastic logger
    tracer = logging.getLogger("elasticsearch")
    tracer.setLevel(logging.CRITICAL + 1)

    # Create indices
    for index in databases + [EXTRA_INDEX]:
        try:
            n_shards = indices[index]
        except KeyError:
            n_shards = shards

        index += suffix

        # https://www.elastic.co/guide/en/elasticsearch/reference/current/indices-update-settings.html
        # https://www.elastic.co/guide/en/elasticsearch/guide/current/indexing-performance.html
        while True:
            try:
                res = es.indices.delete(index)
            except exceptions.ConnectionTimeout:
                pass
            except exceptions.NotFoundError:
                break
            else:
                break

        body = {
            "mappings": {
                doc_type: {
                    "properties": properties
                }
            },
            "settings": {
                "number_of_shards": n_shards,
                "number_of_replicas": 0,    # default: 1
                "refresh_interval": -1      # default: 1s
            }
        }

        es.indices.create(index, body=body)

    queue_in = Queue()
    queue_out = Queue()
    loaders = [
        DocumentLoader(host, doc_type, queue_in, queue_out, suffix=suffix)
        for _ in range(n_loaders)
    ]
    for l in loaders:
        l.start()

    pathname = os.path.join(src, "**", "*.json.gz")
    files = set()
    stop = False
    while True:
        os.path.join(src, LOADING_FILE)
        for filepath in glob.iglob(pathname, recursive=True):
            if filepath not in files:
                files.add(filepath)
                queue_in.put(filepath)

        if stop:
            logging.info("{} files to load".format(len(files)))
            break
        elif not os.path.isfile(os.path.join(src, LOADING_FILE)):
            # All files ready, but loop one last time
            stop = True
        else:
            time.sleep(60)

    # At this point, all files are in the queue
    for _ in loaders:
        queue_in.put(None)

    for l in loaders:
        l.join()

    # Get files that failed to load
    files = []
    for _ in loaders:
        files += queue_out.get()

    # Repeat until all files are loaded
    while files:
        logging.info("{} files to load".format(len(files)))
        queue_in = Queue()
        queue_out = Queue()
        loaders = [
            DocumentLoader(host, doc_type, queue_in, queue_out, suffix=suffix)
            for _ in range(min(n_loaders, len(files)))
        ]
        for l in loaders:
            l.start()

        for filepath in files:
            queue_in.put(filepath)

        for _ in loaders:
            queue_in.put(None)

        for l in loaders:
            l.join()

        files = []
        for _ in loaders:
            files += queue_out.get()

    logging.info("complete")


def update_alias(my_ippro: str, hosts: list, alias: str, **kwargs):
    suffix = kwargs.get("suffix", "").lower()
    delete = kwargs.get("delete", True)

    databases = list(mysql.get_entry_databases(my_ippro).keys())
    new_indices = set()
    for index in databases + [EXTRA_INDEX]:
        new_indices.add(index + suffix)

    for host in hosts:
        es = Elasticsearch([parse_host(host)])

        # Disable Elastic logger
        tracer = logging.getLogger("elasticsearch")
        tracer.setLevel(logging.CRITICAL + 1)

        exists = es.indices.exists_alias(name=alias)
        if exists:
            # Alias already exists: update it

            # Indices currently using the alias
            indices = set(es.indices.get_alias(name=alias))

            actions = []
            for index in new_indices:
                try:
                    indices.remove(index)
                except KeyError:
                    # Index does not yet have this alias: add it
                    actions.append({
                        "add": {
                            "index": index,
                            "alias": alias
                        }
                    })

            for index in indices:
                # Remove outdated indices that have this alias
                actions.append({
                    "remove": {
                        "index": index,
                        "alias": alias
                    }
                })

            if actions:
                # Atomic operation
                # (alias removed from the old indices at the same time it's added to the new ones)
                es.indices.update_aliases(body={"actions": actions})

            if delete:
                for index in indices:
                    while True:
                        try:
                            res = es.indices.delete(index)
                        except exceptions.ConnectionTimeout:
                            pass
                        except exceptions.NotFoundError:
                            break
                        else:
                            break
        else:
            # Create alias
            es.indices.put_alias(index=','.join(new_indices), name=alias)

        # Update index settings
        for index in new_indices:
            es.indices.put_settings({
                # 'number_of_replicas': 1,
                'refresh_interval': None  # default (1s)
            }, index)
