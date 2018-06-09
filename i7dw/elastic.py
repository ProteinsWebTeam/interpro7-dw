#!/usr/bin/env python
# -*- coding: utf-8 -*-

import gzip
import hashlib
import json
import logging
import multiprocessing as mp
import os
import pickle
import shutil
import tempfile
import time

from elasticsearch import Elasticsearch, helpers, exceptions

from . import dbms, disk, mysql
from .ebi import interpro, pdbe

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s: %(levelname)s: %(message)s',
    datefmt='%y-%m-%d %H:%M:%S'
)


def disable_es_logger():
    # since lever is now greater than CRITICAL, no message will be processed
    tracer = logging.getLogger('elasticsearch')
    tracer.setLevel(logging.CRITICAL + 1)


def create_indices(databases, hosts, doc_type, properties_json=None, indices_json=None, default_shards=5, suffix=''):
    # Establish connection
    es = Elasticsearch(hosts)

    # Disable Elastic logger
    disable_es_logger()

    if properties_json:
        with open(properties_json, 'rt') as fh:
            properties = json.load(fh)
    else:
        properties = {}

    if indices_json:
        with open(indices_json, 'rt') as fh:
            custom_shards = json.load(fh)

        # Force keys to be in lower cases
        custom_shards = {k.lower(): custom_shards[k] for k in custom_shards}
    else:
        custom_shards = {}

    for index in databases + ['others']:
        if index in custom_shards:
            shards = custom_shards[index]
        else:
            shards = default_shards

        if suffix:
            index += suffix

        logging.info('creating ES index: {}'.format(index))

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
            'settings': {
                'number_of_shards': shards,
                'number_of_replicas': 0,    # default: 1
                'refresh_interval': -1      # default: 1s
            }
        }

        if properties:
            body['mappings'] = {
                doc_type: {
                    'properties': properties
                }
            }

        es.indices.create(index, body=body)


class ElasticDocProducer(mp.Process):
    def __init__(self, ora_uri, my_uri, task_queue, file_queue, sm_queue, outdir, **kwargs):
        super().__init__()
        self.ora_uri = ora_uri
        self.my_uri = my_uri
        self.task_queue = task_queue
        self.file_queue = file_queue
        self.sm_queue = sm_queue
        self.outdir = outdir
        self.min_overlap = kwargs.get('min_overlap', 20)
        self.max_size = kwargs.get('max_size', 1000000)
        self.chunk_size = kwargs.get('chunk_size', 10000)

        self.entries = None
        self.sets = None
        self.proteomes = None
        self.pfam = set()
        self.structures = None

    def run(self):
        logging.info('{} ({}) started'.format(self.name, os.getpid()))

        # Get PDBe structures
        self.structures = pdbe.get_structures(self.ora_uri, citations=True, fragments=True, by_protein=True)

        # Get entries/sets/proteomes
        self.entries = mysql.get_entries(self.my_uri)
        self.sets = mysql.get_sets(self.my_uri)
        self.proteomes = mysql.get_proteomes(self.my_uri)

        for e in self.entries.values():
            if e['database'] == 'pfam':
                self.pfam.add(e['accession'])

        docs = []
        cnt = 0
        ts = time.time()

        # Get proteins from main queue
        while True:
            task = self.task_queue.get()
            if task is None:
                break

            _type, chunk = task
            if _type == 'protein':
                fn = self.process_protein
            elif _type == 'entry':
                fn = self.process_entry
            elif _type == 'taxonomy':
                fn = self.process_taxonomy
            else:
                continue

            for args in chunk:
                docs += fn(*args)

                if len(docs) >= self.max_size:
                    cnt += len(docs)
                    for filepath in self._dump(docs, self.outdir, self.chunk_size):
                        if self.file_queue:
                            self.file_queue.put(filepath)
                    docs = []

        if docs:
            cnt += len(docs)
            for filepath in self._dump(docs, self.outdir, self.chunk_size):
                if self.file_queue:
                    self.file_queue.put(filepath)
            docs = []

        # Free memory (if process terminates but not garbage collected)
        self.entries = None
        self.sets = None
        self.proteomes = None
        self.pfam = set()
        self.structures = None

        logging.info('{} ({}) terminated ({} documents, {:.0f} documents/sec)'.format(
            self.name, os.getpid(), cnt, cnt // (time.time() - ts))
        )

    @staticmethod
    def init_doc():
        return {
            'id': None,

            # Protein
            'protein_acc': None,
            'protein_length': None,
            'protein_size': None,
            'protein_db': None,
            'text_protein': None,

            # Taxonomy
            'tax_id': None,
            'tax_name': None,
            'lineage': None,
            'rank': None,
            'text_taxonomy': None,

            # Proteome
            'proteome_acc': None,
            'proteome_name': None,
            'is_reference': None,
            'text_proteome': None,

            # Entry
            'entry_acc': None,
            'entry_db': None,
            'entry_type': None,
            'entry_date': None,
            'entry_protein_locations': None,
            'entry_go_terms': None,
            'integrated': None,
            'text_entry': None,

            # Set
            'set_acc': None,
            'set_db': None,
            'set_integrated': None,
            'text_set': None,

            # Structure
            'structure_acc': None,
            'structure_resolution': None,
            'structure_date': None,
            'structure_evidence': None,

            # Chain
            'chain': None,
            'protein_structure_locations': None,
            'structure_chain': None,
            'text_structure': None,

            # Domain architecture
            'ida_id': None,
            'ida': None
        }

    @staticmethod
    def _join(*args, separator=' '):
        items = []
        for item in args:
            if item is None:
                continue
            elif isinstance(item, (int, float)):
                item = str(item)
            elif isinstance(item, (list, tuple, set)):
                item = separator.join(map(str, item))
            elif isinstance(item, dict):
                item = separator.join(map(str, item.values()))
            elif not isinstance(item, str):
                continue

            items.append(item)

        return separator.join(items)

    def process_entry(self, accession):
        entry = self.entries[accession]
        doc = self.init_doc()

        go_terms = [t['identifier'] for t in entry['go_terms']]
        doc.update({
            'entry_acc': entry['accession'],
            'entry_db': entry['database'],
            'entry_type': entry['type'],
            'entry_date': entry['date'].strftime('%Y-%m-%d'),
            'integrated': entry['integrated'],
            'text_entry': self._join(
                entry['accession'], entry['name'], entry['type'], entry['descriptions'], *go_terms
            ),
            'entry_protein_locations': [],  # this entry does not match any protein
            'entry_go_terms': go_terms
        })

        sets = self.sets.get(accession)

        if sets:
            docs = []
            for set_ac, set_db in sets.items():
                _doc = doc.copy()
                _doc.update({
                    'set_acc': set_ac,
                    'set_db': set_db,
                    'set_integrated': [],  # todo: implement set integration (e.g. pathways)
                    'text_set': self._join(set_ac, set_db),
                    'id': self._join(entry['accession'], set_ac, separator='-')
                })
                docs.append(_doc)

            return docs
        else:
            doc['id'] = entry['accession']
            return [doc]

    def process_taxonomy(self, taxon):
        doc = self.init_doc()

        doc.update({
            'tax_id': taxon['taxId'],
            'tax_name': taxon['scientificName'],
            'lineage': taxon['lineage'].strip().split(),
            'rank': taxon['rank'],
            'text_taxonomy': self._join(taxon['taxId'], taxon['fullName'], taxon['rank']),
            'id': taxon['taxId']
        })

        return [doc]

    def process_protein(self, accession, identifier, name, database, length, comments, matches, proteomes, taxon):
        # PDBe structures for this protein
        structures = self.structures.get(accession, {})

        # Prepare matches/supermatches
        entry_matches = {}
        supermatches = []
        dom_arch = []
        dom_entries = set()
        for m in matches:
            entry_ac = m['entry_ac']
            method_ac = m['method_ac']

            if method_ac not in entry_matches:
                entry_matches[method_ac] = []

            entry_matches[method_ac].append((m['start'], m['end']))

            if method_ac in self.pfam:
                dom_entries.add(method_ac)
                if entry_ac:
                    dom_entries.add(entry_ac)
                    dom_arch.append('{}:{}'.format(method_ac, entry_ac))
                else:
                    dom_arch.append('{}:'.format(method_ac))

            if entry_ac:
                try:
                    entry = self.entries[entry_ac]
                except KeyError:
                    # TODO: remove after next refresh
                    continue

                supermatches.append(
                    interpro.Supermatch(entry_ac, entry['root'], m['start'], m['end'])
                )

        if dom_arch:
            dom_arch = '-'.join(dom_arch)
            dom_arch_id = hashlib.sha1(dom_arch.encode('utf-8')).hexdigest()
        else:
            dom_arch = dom_arch_id = None

        # Merge overlapping supermatches
        records = []
        sets = interpro.merge_supermatches(supermatches, min_overlap=self.min_overlap)
        for s in sets:
            for sm in s.supermatches:
                for entry_ac in sm.get_entries():
                    # Supermatch records for Oracle
                    records.append((accession, entry_ac, sm.start, sm.end, 'S' if database == 'reviewed' else 'T'))

                    # Add supermatches to matches for Elastic
                    if entry_ac not in entry_matches:
                        entry_matches[entry_ac] = []

                    entry_matches[entry_ac].append((sm.start, sm.end))

        if self.sm_queue and records:
            self.sm_queue.put(records)

        if length <= 100:
            size = 'small'
        elif length <= 1000:
            size = 'medium'
        else:
            size = 'large'

        doc = self.init_doc()

        # Accession in lower cases
        accession = accession.lower()

        # Add protein info
        doc.update({
            'protein_acc': accession,
            'protein_length': length,
            'protein_size': size,
            'protein_db': database,
            'text_protein': self._join(accession, identifier, name, database, comments),

            'tax_id': taxon['taxId'],
            'tax_name': taxon['scientificName'],
            'lineage': taxon['lineage'].strip().split(),
            'rank': taxon['rank'],
            'text_taxonomy': self._join(taxon['taxId'], taxon['fullName'], taxon['rank']),
        })

        docs = [doc]
        _docs = []
        for upid in proteomes:
            if upid in self.proteomes:
                p = self.proteomes[upid]
            else:
                logging.warning('invalid proteome {} for protein {}'.format(upid, accession))
                continue

            for doc in docs:
                _doc = doc.copy()
                _doc.update({
                    'proteome_acc': upid,
                    'proteome_name': p['name'],
                    'is_reference': p['is_reference'],
                    'text_proteome': self._join(upid, *list(p.values())),
                })
                _docs.append(_doc)

        if _docs:
            docs = _docs
            _docs = []

        for entry_ac in entry_matches:
            entry = self.entries[entry_ac]
            go_terms = [t['identifier'] for t in entry['go_terms']]

            for doc in docs:
                _doc = doc.copy()
                _doc.update({
                    'entry_acc': entry['accession'],
                    'entry_db': entry['database'],
                    'entry_type': entry['type'],
                    'entry_date': entry['date'].strftime('%Y-%m-%d'),
                    'integrated': entry['integrated'],
                    'text_entry': self._join(
                        entry['accession'], entry['name'], entry['type'], entry['descriptions'], *go_terms
                    ),
                    'entry_protein_locations': [
                        # todo: implement multi part fragments
                        {'fragments': [{'start': start, 'end': end}]} for start, end in entry_matches[entry_ac]
                    ],
                    'entry_go_terms': go_terms
                })

                if entry['accession'] in dom_entries:
                    # IDA makes only sense for InterPro/Pfam entries
                    _doc.update({
                        'ida_id': dom_arch_id,
                        'ida': dom_arch
                    })

                sets = self.sets.get(entry_ac)
                if sets:
                    for set_ac, set_db in sets.items():
                        _doc_set = _doc.copy()
                        _doc_set.update({
                            'set_acc': set_ac,
                            'set_db': set_db,
                            'set_integrated': [],  # todo: implement set integration (e.g. pathways)
                            'text_set': self._join(set_ac, set_db)
                        })
                        _docs.append(_doc_set)
                else:
                    _docs.append(_doc)

        if _docs:
            docs = _docs
            _docs = []

        for structure in structures.values():
            text = self._join(
                structure['accession'],
                structure['evidence'],
                structure['name'],
                ' '.join([pub['title'] for pub in structure['citations'].values() if pub.get('title')])
            )

            for doc in docs:
                _doc = doc.copy()
                _doc.update({
                    'structure_acc': structure['accession'],
                    'structure_resolution': structure['resolution'],
                    'structure_date': structure['date'].strftime('%Y-%m-%d'),
                    'structure_evidence': structure['evidence']
                })

                for chain, fragments in structure['chains'].items():
                    _doc_chain = _doc.copy()
                    _doc_chain.update({
                        'chain': chain,
                        'protein_structure_locations': [
                            {'fragments': [{'start': m['start'], 'end': m['end']}]} for m in fragments
                        ],
                        'structure_chain': '{} - {}'.format(structure['accession'], chain),
                        'text_structure': '{} {}'.format(chain, text)
                    })

                    _docs.append(_doc_chain)

        if _docs:
            docs = _docs
            _docs = []

        for doc in docs:
            doc['id'] = self._join(
                doc['protein_acc'], doc['proteome_acc'], doc['entry_acc'],
                doc['set_acc'], doc['structure_acc'], doc['chain'],
                separator='-'
            )

        return docs

    @staticmethod
    def _dump(docs, outdir, chunk_size):
        if len(docs) > chunk_size:
            # Too many documents for one single file: create a directory and write files inside
            outdir = tempfile.mkdtemp(dir=outdir)

            for i in range(0, len(docs), chunk_size):
                fd, filepath = tempfile.mkstemp(suffix='.json', dir=outdir)
                os.close(fd)

                with open(filepath, 'wt') as fh:
                    json.dump(docs[i:i+chunk_size], fh)

                yield filepath
        else:
            # All documents in one single file
            fd, filepath = tempfile.mkstemp(suffix='.json', dir=outdir)
            os.close(fd)

            with open(filepath, 'wt') as fh:
                json.dump(docs, fh)

            yield filepath


class OracleLoader(mp.Process):
    def __init__(self, ora_uri, my_uri, task_queue, chunk_size=100000):
        super().__init__()
        self.ora_uri = ora_uri
        self.my_uri = my_uri
        self.task_queue = task_queue
        self.chunk_size = chunk_size

    def run(self):
        logging.info('{} ({}) started'.format(self.name, os.getpid()))

        con, cur = dbms.connect(self.ora_uri)
        try:
            cur.execute('DROP TABLE INTERPRO.SUPERMATCH2 CASCADE CONSTRAINTS')
        except:
            pass

        cur.execute(
            """
            CREATE TABLE INTERPRO.SUPERMATCH2
            (
                PROTEIN_AC VARCHAR2(15) NOT NULL,
                ENTRY_AC VARCHAR2(9) NOT NULL,
                POS_FROM NUMBER(5) NOT NULL,
                POS_TO NUMBER(5) NOT NULL,
                DBCODE CHAR(1) NOT NULL 
            ) NOLOGGING
            """
        )

        cnt = 0         # num of supermatches
        data = []       # chunk of supermatches to be inserted
        sets = {}       # entry_ac -> num of proteins matched
        overlaps = {}   # entry_ac1 -> entry_ac2 -> [0, 0] (num of proteins in which on entry overlaps with the other)
        ts = time.time()
        while True:
            records = self.task_queue.get()

            if records is None:
                break

            matches = {}  # entry_ac -> [(start1, end1), (start2, end2), ...]

            for protein_ac, entry_ac, start, end, dbcode in records:
                cnt += 1
                data.append((protein_ac.upper(), entry_ac.upper(), start, end, dbcode))

                if len(data) == self.chunk_size:
                    cur.executemany(
                        """
                        INSERT /*+APPEND*/ INTO INTERPRO.SUPERMATCH2 (PROTEIN_AC, ENTRY_AC, POS_FROM, POS_TO, DBCODE)
                        VALUES (:1, :2, :3, :4, :5)
                        """,
                        data
                    )
                    con.commit()
                    data = []

                # Current implementation only considers the leftmost match
                if entry_ac not in matches:
                    matches[entry_ac] = [(start, end)]

            interpro.intersect_matches(matches, sets, overlaps)

        if data:
            cur.executemany(
                """
                INSERT /*+APPEND*/ INTO INTERPRO.SUPERMATCH2 (PROTEIN_AC, ENTRY_AC, POS_FROM, POS_TO, DBCODE)
                VALUES (:1, :2, :3, :4, :5)
                """,
                data
            )
            con.commit()

        logging.info('{} ({}): {} supermatches inserted ({:.0f} supermatches/sec)'.format(
            self.name, os.getpid(), cnt, cnt // (time.time() - ts))
        )

        # Constraints
        logging.info('{} ({}): adding constraints'.format(self.name, os.getpid()))
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
                FOREIGN KEY (PROTEIN_AC) REFERENCES INTERPRO.PROTEIN (PROTEIN_AC)
                ON DELETE CASCADE
                """
            )
        except:
            logging.error('{} ({}): could not create PROTEIN_AC constraint on INTERPRO.SUPERMATCH2'.format(
                self.name, os.getpid()
            ))

        try:
            cur.execute(
                """
                ALTER TABLE INTERPRO.SUPERMATCH2
                ADD CONSTRAINT FK_SUPERMATCH2$ENTRY_AC
                FOREIGN KEY (ENTRY_AC) REFERENCES INTERPRO.ENTRY (ENTRY_AC)
                ON DELETE CASCADE
                """
            )
        except:
            logging.error('{} ({}): could not create ENTRY_AC constraint on INTERPRO.SUPERMATCH2'.format(
                self.name, os.getpid()
            ))

        # Indexes
        logging.info('{} ({}): creating indexes'.format(self.name, os.getpid()))
        cur.execute('CREATE INDEX I_SUPERMATCH2$PROTEIN ON INTERPRO.SUPERMATCH2 (PROTEIN_AC) NOLOGGING')
        cur.execute('CREATE INDEX I_SUPERMATCH2$ENTRY ON INTERPRO.SUPERMATCH2 (ENTRY_AC) NOLOGGING')
        cur.execute('CREATE INDEX I_SUPERMATCH2$DBCODE$ENTRY ON INTERPRO.SUPERMATCH2 (DBCODE, ENTRY_AC) NOLOGGING')

        # Statistics
        logging.info('{} ({}): gathering statistics'.format(self.name, os.getpid()))
        cur.execute(
            """
                BEGIN
                    DBMS_STATS.GATHER_TABLE_STATS(:1, :2, cascade => TRUE);
                END;
            """,
            ('INTERPRO', 'SUPERMATCH2')
        )

        # Privileges
        cur.execute('GRANT SELECT ON INTERPRO.SUPERMATCH2 TO INTERPRO_SELECT')

        cur.close()
        con.close()

        # Compute Jaccard coefficients
        overlapping = {}

        for acc1 in overlaps:
            s1 = sets[acc1]

            for acc2 in overlaps[acc1]:
                s2 = sets[acc2]
                o1, o2 = overlaps[acc1][acc2]

                # Independent coefficients (does not have much sense, but hey, I did not come with this)
                coef1 = o1 / (s1 + s2 - o1)
                coef2 = o2 / (s1 + s2 - o2)

                # Final coefficient: average of independent coefficients
                coef = (coef1 + coef2) * 0.5

                # Containment indices
                c1 = o1 / s1
                c2 = o2 / s2

                # TODO: do not hardcode threshold
                if coef >= 0.75 or c1 >= 0.51 or c2 >= 0.51:
                    if acc1 in overlapping:
                        overlapping[acc1].append(acc2)
                    else:
                        overlapping[acc1] = [acc2]

                    if acc2 in overlapping:
                        overlapping[acc2].append(acc1)
                    else:
                        overlapping[acc2] = [acc1]

        con, cur = dbms.connect(self.my_uri)

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

        logging.info('{} ({}) terminated'.format(self.name, os.getpid()))


class _ElasticLoader(mp.Process):
    def __init__(self, hosts, doc_type, inqueue, outqueue, **kwargs):
        super().__init__()
        self.hosts = hosts
        self.type = doc_type
        self.inqueue = inqueue
        self.outqueue = outqueue
        self.gzip = kwargs.get('gzip', False)
        self.max_bytes = kwargs.get('max_bytes', 100 * 1024 * 1024)
        self.suffix = kwargs.get('suffix')

    def run(self):
        es = Elasticsearch(self.hosts)

        # Disable logger (if not already disabled by parent process)
        disable_es_logger()

        while True:
            filepath = self.inqueue.get()

            if filepath is None:
                break
            elif filepath.lower().endswith('.gz'):
                _open = gzip.open
                compress = False  # file is already compressed
            else:
                _open = open
                compress = self.gzip

            with _open(filepath, 'rt') as fh:
                docs = json.load(fh)

            actions = []
            for doc in docs:
                if doc['entry_db']:
                    index = doc['entry_db']
                else:
                    index = 'others'

                if self.suffix:
                    index += self.suffix

                actions.append({
                    '_op_type': 'index',
                    '_index': index,
                    '_type': self.type,
                    '_id': doc['id'],
                    '_source': doc
                })

            ts = time.time()
            n_successfull, errors = helpers.bulk(
                es, actions,
                chunk_size=-1,  # disable chunk_size (num of docs) to only rely on max_chunk_bytes (bytes)
                max_chunk_bytes=self.max_bytes,
                raise_on_exception=False,
                raise_on_error=False
            )
            index_time = time.time() - ts

            if n_successfull == len(actions):
                if compress:
                    with gzip.open(filepath + '.gz', 'wt') as fh:
                        json.dump(docs, fh)

                    os.unlink(filepath)

                if os.path.isfile(filepath + '.p'):
                    os.unlink(filepath + '.p')

                self.outqueue.put((filepath, True, index_time))
            else:
                with open(filepath + '.p', 'wb') as fh:
                    pickle.dump(errors, fh)

                self.outqueue.put((filepath, False, index_time))


class ElasticLoaderPool(mp.Process):
    def __init__(self, hosts, doc_type, inqueue, **kwargs):
        super().__init__()
        self.hosts = hosts
        self.type = doc_type
        self.inqueue = inqueue

        self.gzip = kwargs.get('gzip', False)
        self.processes = kwargs.get('processes', 3)
        self.outqueue = kwargs.get('outqueue')
        self.suffix= kwargs.get('suffix')

    def run(self):
        logging.info('{} ({}) started'.format(self.name, os.getpid()))

        """
        Queues for loaders
        in: filepath
        out: tuple (filepath, status, index_time)
        """
        inqueue = mp.Queue()
        outqueue = mp.Queue()

        # Starts with a chunk size of 100MB
        max_bytes = 100 * 1024 * 1024

        # Starts with one loader only
        processes = 1
        loaders = [
            _ElasticLoader(self.hosts, self.type, inqueue, outqueue,
                           gzip=self.gzip, max_bytes=max_bytes, suffix=self.suffix)
        ]

        loaders[0].start()

        cnt = 0
        n_indexed = 0     # number of completely indexed files
        _n_indexed = 0
        avg_time = 0    # average time to (try to) index a file
        _avg_time = 0
        failed_files = []
        while True:
            filepath = self.inqueue.get()

            if filepath is None:
                break
            elif self.hosts is None:
                continue

            inqueue.put(filepath)
            cnt += 1

            if cnt == 10:
                for _ in loaders:
                    inqueue.put(None)

                for l in loaders:
                    l.join()

                for _ in range(cnt):
                    filepath, status, index_time = outqueue.get()
                    avg_time += index_time

                    if status:
                        n_indexed += 1
                    else:
                        failed_files.append(filepath)

                avg_time /= cnt
                logging.info(
                    '{}/{} files indexed (avg: {:.1f} secs); workers: {}/{}'.format(
                        n_indexed, cnt, avg_time, processes, self.processes
                    )
                )

                if n_indexed == cnt and processes < self.processes:
                    # All documents were indexed, and we are not using all workers: add one
                    processes += 1
                elif n_indexed < 0.75 * cnt and processes > 1:
                    # More than 25% of files failed to load
                    processes -= 1

                _n_indexed = n_indexed
                _avg_time = avg_time
                n_indexed = 0
                cnt = 0

                loaders = [
                    _ElasticLoader(self.hosts, self.type, inqueue, outqueue,
                                   gzip=self.gzip, max_bytes=max_bytes, suffix=self.suffix)
                    for _ in range(processes)
                ]

                for l in loaders:
                    l.start()

        if cnt:
            for _ in loaders:
                inqueue.put(None)
            inqueue.close()

            for l in loaders:
                l.join()

            for _ in range(cnt):
                filepath, status, index_time = outqueue.get()
                avg_time += index_time

                if status:
                    n_indexed += 1
                else:
                    failed_files.append(filepath)

            avg_time /= cnt
            logging.info(
                '{}/{} files indexed (avg: {:.1f} secs); workers: {}/{}'.format(
                    n_indexed, cnt, avg_time, processes, self.processes
                )
            )

        if self.outqueue:
            # Pass the file that could not be completely indexed to the out queue
            self.outqueue.put(failed_files)


class ElasticLoader(mp.Process):
    def __init__(self, hosts, doc_type, version, task_queue, **kwargs):
        super().__init__()
        self.hosts = hosts
        self.type = doc_type
        self.version = version
        self.task_queue = task_queue

        self.chunk_size = kwargs.get('chunk_size', 500)
        self.retries = kwargs.get('retries', 0)
        self.compress = kwargs.get('gzip', False)
        # queue for files containing at least one document that could not be indexed
        self.errors = kwargs.get('errors')

    def run(self):
        logging.info('{} ({}) started'.format(self.name, os.getpid()))

        es = Elasticsearch(self.hosts)

        # Disable logger (if not already disabled by parent process)
        disable_es_logger()

        n_indexed = 0
        n_errors = 0
        failed_files = []
        ts = time.time()
        while True:
            filepath = self.task_queue.get()

            if filepath is None:
                break
            elif self.hosts is None:
                continue
            elif filepath.lower().endswith('.gz'):
                _open = gzip.open
                compress = False    # file is already compressed
            else:
                _open = open
                compress = self.compress

            with _open(filepath, 'rt') as fh:
                docs = json.load(fh)

            actions = []
            for doc in docs:
                if doc['entry_db']:
                    index = doc['entry_db'] + self.version
                else:
                    index = 'others' + self.version

                actions.append({
                    '_op_type': 'index',
                    '_index': index,
                    '_type': self.type,
                    '_id': doc['id'],
                    '_source': doc
                })

            n_successfull, errors = helpers.bulk(
                es, actions,
                chunk_size=self.chunk_size,
                max_chunk_bytes=100 * 1024 * 1024,
                raise_on_exception=False,
                raise_on_error=False,
                max_retries=self.retries
            )

            logging.info('{} ({}): {}: {} indexed, {} errors'.format(
                self.name, os.getpid(), os.path.basename(filepath), n_successfull, len(errors))
            )

            n_indexed += n_successfull
            n_errors += len(errors)

            if n_successfull == len(actions):
                if compress:
                    with gzip.open(filepath + '.gz', 'wt') as fh:
                        json.dump(docs, fh)

                    os.unlink(filepath)

                if os.path.isfile(filepath + '.p'):
                    os.unlink(filepath + '.p')
            else:
                with open(filepath + '.p', 'wb') as fh:
                    pickle.dump(errors, fh)

                failed_files.append(filepath)

        if self.errors:
            self.errors.put(failed_files)

        logging.info('{} ({}) terminated ({} documents, {:.0f} documents/sec, {} errors in {} files)'.format(
            self.name, os.getpid(), n_indexed, n_indexed // (time.time() - ts), n_errors, len(failed_files))
        )


def create_relationships(ora_uri, my_uri, proteins_f, descriptions_f, comments_f, proteomes_f, prot_matches_f,
                         outdir, **kwargs):

    # Workers (# of total proc: +1 for parent, +1 for OracleLoader, +1 for ElasticLoaderPool)
    n_producers = kwargs.get('producers', 3)
    n_loaders = kwargs.get('loaders', 3)

    chunk_size = kwargs.get('chunk_size', 100000)
    limit = kwargs.get('limit', 0)
    supermatches = kwargs.get('supermatches', False)

    # Indices-related stuff
    hosts = kwargs.get('hosts')
    doc_type = kwargs.get('doc_type')
    suffix = kwargs.get('suffix')
    indices_json = kwargs.get('indices_json')
    properties_json = kwargs.get('properties_json')
    default_shards = kwargs.get('shards', 5)

    # Protein queue: passed to child process that write JSON of relationships
    protein_queue = mp.Queue(n_producers)

    if supermatches:
        # Supermatches, then overlapping entries (Jaccard) will be calculated
        sm_queue = mp.Queue()
        ora_loader = OracleLoader(ora_uri, my_uri, sm_queue)
        ora_loader.start()
    else:
        sm_queue = None
        ora_loader = None

    if hosts and n_loaders:
        # Generated JSON files will be indexed into Elasticsearch
        file_queue = mp.Queue()
        es_loader = ElasticLoaderPool(hosts, doc_type, file_queue, gzip=True, processes=n_loaders, suffix=suffix)
        es_loader.start()
    else:
        file_queue = None
        es_loader = None

    producers = [
        ElasticDocProducer(ora_uri, my_uri, protein_queue, file_queue, sm_queue, outdir)
        for _ in range(n_producers)
    ]

    for w in producers:
        w.start()

    # MySQL data
    logging.info('loading data from MySQL')
    taxa = mysql.get_taxa(my_uri, slim=False)
    entries = mysql.get_entries(my_uri, minimal=True)

    # Create output directory for files
    logging.info('preparing output directory')
    if os.path.isdir(outdir):
        shutil.rmtree(outdir)
    os.makedirs(outdir)

    if es_loader:
        # Create indices (delete them first if they exist)
        databases = list(mysql.get_entry_databases(my_uri).keys())
        create_indices(databases, hosts, doc_type, properties_json, indices_json, default_shards, suffix)

    # Open stores
    proteins = disk.Store(proteins_f)
    descriptions = disk.Store(descriptions_f)
    comments = disk.Store(comments_f)
    proteomes = disk.Store(proteomes_f)
    prot_matches = disk.Store(prot_matches_f)

    set_taxa = set(taxa.keys())

    logging.info('starting')
    cnt = 0
    chunk = []
    ts = time.time()
    for acc, protein in proteins.iter():
        tax_id = protein['taxon']

        try:
            taxon = taxa[tax_id]
        except KeyError:
            logging.warning('invalid taxon ({}) for protein {}'.format(protein['taxon'], acc))
            continue

        if tax_id in set_taxa:
            set_taxa.remove(tax_id)

        name, other_names = descriptions.get(acc, (None, None))
        matches = prot_matches.get(acc, [])

        # Enqueue data for creating Elastic documents
        chunk.append((
            acc,
            protein['identifier'],
            name,
            'reviewed' if protein['isReviewed'] else 'unreviewed',
            protein['length'],
            comments.get(acc, []),
            matches,
            proteomes.get(acc, []),
            taxon
        ))

        if len(chunk) == chunk_size:
            protein_queue.put(('protein', chunk))
            chunk = []

        # Remove entries with protein matches
        _entries = []
        for m in matches:
            _entries.append(m['method_ac'])
            if m['entry_ac']:
                _entries.append(m['entry_ac'])

        entries -= set(_entries)

        cnt += 1
        if not cnt % 1000000:
            logging.info('{:>12} ({:.0f} proteins/sec)'.format(cnt, cnt // (time.time() - ts)))

        if cnt == limit:
            break

    logging.info('{:>12} ({:.0f} proteins/sec)'.format(cnt, cnt // (time.time() - ts)))

    if chunk:
        protein_queue.put(('protein', chunk))

    # Add entries without matches
    chunk = [(entry_ac,) for entry_ac in entries]
    for i in range(0, len(chunk), chunk_size):
        protein_queue.put(('entry', chunk[i:i+chunk_size]))

    # Add taxa without proteins
    chunk = [(taxa[tax_id],) for tax_id in set_taxa]
    for i in range(0, len(chunk), chunk_size):
        protein_queue.put(('taxonomy', chunk[i:i+chunk_size]))

    # All proteins have been enqueued: poison pill
    for _ in producers:
        protein_queue.put(None)
    protein_queue.close()

    # Close stores to free memory
    for store in (proteins, descriptions, comments, proteomes, prot_matches):
        store.close()

    # Wait for ElasticDocProducers to finish
    for w in producers:
        w.join()

    # When ElasticDocProducers are done, no more file can be passed to the file queue: poison pill
    if file_queue:
        file_queue.put(None)
        file_queue.close()

    # And same for supermatches (this will trigger constraints/indexes creation)
    if sm_queue:
        sm_queue.put(None)
        sm_queue.close()

    if es_loader:
        # Wait for ElasticLoader to terminate (some files still could have failed to be indexed)
        es_loader.join()

        # Then try to index files that failed earlier
        index_documents(hosts, doc_type, outdir, processes=n_loaders, extension='.json', gzip=True, suffix=suffix)

    # index_documents() may take a while (if being run): in the meantime let's wait for supermatches to be inserted
    if ora_loader:
        ora_loader.join()

    logging.info('complete')


def find_files(src, extension, limit=0):
    paths = []
    for root, dirs, files in os.walk(src):
        for f in files:
            if f.endswith(extension):
                paths.append(os.path.join(root, f))

                if len(paths) == limit:
                    return paths

    return paths


def index_documents(hosts, doc_type, src, **kwargs):
    extension = kwargs.get('extension', '.json')
    files = kwargs.get('files')
    limit = kwargs.get('limit', 0)

    # Any additional keyword arguments will be passed to ElasticLoaderPool constructor (e.g. processes, gzip, suffix)

    if not files:
        # Get files from source directory
        files = find_files(src, extension, limit)

        if not files:
            logging.warning('0 files to load')
            return

    logging.info('{} files to load'.format(len(files)))

    inqueue = mp.Queue()    # files to index
    outqueue = mp.Queue()   # files that failed to be indexed

    loader = ElasticLoaderPool(hosts, doc_type, inqueue, outqueue=outqueue, **kwargs)
    loader.start()

    for filepath in files:
        inqueue.put(filepath)

    inqueue.put(None)
    inqueue.close()
    loader.join()
    files = outqueue.get()

    return files


def collect(uri, hosts, doc_type, src, **kwargs):
    n_iter = kwargs.get('iterations', 3)
    seconds = kwargs.get('seconds', 60)
    suffix = kwargs.get('suffix')

    # To create indices
    _create_indices = kwargs.get('create_indices', False)
    indices_json = kwargs.get('indices_json')
    properties_json = kwargs.get('properties_json')
    default_shards = kwargs.get('shards', 5)

    # Define indices using the entry databases and the version
    databases = list(mysql.get_entry_databases(uri).keys())

    if _create_indices:
        create_indices(databases, hosts, doc_type, properties_json, indices_json, default_shards, suffix)

    # Add the "other" index for documents without an entry (hence without DB)
    indices = databases + ['others']

    # First iteration (files that could not be completely loaded are returned)
    files = index_documents(hosts, doc_type, src, **kwargs)

    i_iter = 1  # start at 1 to perform (n_iter - 1) iterations (as we already performed one above)
    while files and i_iter < n_iter:
        time.sleep(seconds)  # let Elastic some time to complete tasks
        i_iter += 1

        # Force properties to be None, to avoid indices to be deleted and re-created
        kwargs['properties'] = None
        kwargs['files'] = files  # files to be loaded
        files = index_documents(hosts, doc_type, src, **kwargs)

    if files:
        raise RuntimeError('could not load all files')

    # All files loaded at this point: it's time to update settings and create an alias
    es = Elasticsearch(hosts)

    # Disable Elastic logger
    disable_es_logger()

    # Update indices settings
    for index in indices:
        if suffix:
            index += suffix

        es.indices.put_settings({
            # 'number_of_replicas': 1,
            'refresh_interval': None  # default (1s)
        }, index)

        # # Perform a force merge to reduce the number of segments (disabled as painfully slow)
        # es.indices.forcemerge(index)

    logging.info('complete')


def update_alias(hosts, version, alias, indices=None, uri=None):
    if not indices:
        databases = mysql.get_entry_databases(uri)
        indices = list(databases.keys()) + ['others']

    es = Elasticsearch(hosts)
    disable_es_logger()

    exists = es.indices.exists_alias(name=alias)
    if exists:
        actions = []

        for index in indices:
            actions.append({'add': {'index': index + version, 'alias': alias}})

        # Indices currently using the alias
        indices = es.indices.get_alias(name=alias)
        for index in indices:
            actions.append({'remove': {'index': index, 'alias': alias}})

        # Atomic operation (alias removed from the old indices at the same time it's added to the new ones)
        es.indices.update_aliases(body={'actions': actions})
    else:
        # Create alias
        es.indices.put_alias(index=','.join([index + version for index in indices]), name=alias)
