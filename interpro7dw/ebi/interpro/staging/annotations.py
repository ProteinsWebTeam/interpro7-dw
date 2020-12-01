import gzip
import json
import os
from io import StringIO
from multiprocessing import Process, Queue

import MySQLdb

from interpro7dw import hmmer, logger
from interpro7dw.ebi import pfam
from interpro7dw.ebi.interpro import production as ippro
from interpro7dw.utils import DumpFile, DirectoryTree, Store
from interpro7dw.utils import url2dict


def _export_hmms(p_uniprot2matches: str, pro_url: str, dt: DirectoryTree,
                 buffer_size: int = 1000):
    logger.info("counting hits per model")
    signatures = {}
    with Store(p_uniprot2matches) as u2matches:
        cnt = 0
        for entries in u2matches.values():
            for entry_acc, locations in entries.items():
                for loc in locations:
                    if loc["model"] is None:
                        continue  # InterPro entries

                    try:
                        models = signatures[entry_acc]
                    except KeyError:
                        models = signatures[entry_acc] = {}

                    try:
                        models[loc["model"]] += 1
                    except KeyError:
                        models[loc["model"]] = 1

            cnt += 1
            if not cnt % 10e6:
                logger.info(f"{cnt:>12,}")

        logger.info(f"{cnt:>12,}")

    for entry_acc, models in signatures.items():
        # Select the model with the most hits
        model_acc = sorted(models, key=lambda k: (-models[k], k))[0]
        signatures[entry_acc] = model_acc

    df = DumpFile(dt.mktemp(), compress=True)
    cnt = 0

    iterator = ippro.get_hmms(pro_url, multi_models=True)
    for entry_acc, model_acc, hmm_bytes in iterator:
        if model_acc and signatures[entry_acc] != model_acc:
            continue

        hmm_str = gzip.decompress(hmm_bytes).decode("utf-8")
        df.dump((
            entry_acc,
            "hmm",
            hmm_bytes,
            "application/gzip",
            None
        ))

        with StringIO(hmm_str) as stream:
            hmm = hmmer.HMMFile(stream)

        df.dump((
            entry_acc,
            "logo",
            json.dumps(hmm.logo("info_content_all", "hmm")),
            "application/json",
            None
        ))

        cnt += 2
        if cnt >= buffer_size:
            df.close()
            yield df.path
            df = DumpFile(dt.mktemp(), compress=True)
            cnt = 0

    df.close()
    yield df.path


def _export_alns(pfam_url: str, dt: DirectoryTree, buffer_size: int = 1000):
    df = DumpFile(dt.mktemp(), compress=True)
    cnt = 0

    iterator = pfam.get_alignments(pfam_url)
    for entry_acc, aln_type, aln_bytes, count in iterator:
        df.dump((
            entry_acc,
            f"alignment:{aln_type}",
            aln_bytes,
            "application/gzip",
            count
        ))

        cnt += 1
        if cnt == buffer_size:
            df.close()
            yield df.path
            df = DumpFile(dt.mktemp(), compress=True)
            cnt = 0

    df.close()
    yield df.path


def _insert(url: str, queue: Queue):
    for path in iter(queue.get, None):
        with DumpFile(path) as df:
            con = MySQLdb.connect(**url2dict(url))
            cur = con.cursor()

            for acc, anntype, value, mime, count in df:
                cur.execute(
                    """
                        INSERT INTO webfront_entryannotation (
                          accession, type, value, mime_type, num_sequences
                        )
                        VALUES (%s, %s, %s, %s, %s)
                    """,
                    (acc, anntype, value, mime, count)
                )

            con.commit()
            cur.close()
            con.close()

        os.remove(path)


def insert_annotations(pro_url: str, p_uniprot2matches: str, pfam_url: str,
                       stg_url: str, **kwargs):
    tmpdir = kwargs.get("tmpdir")

    con = MySQLdb.connect(**url2dict(stg_url))
    cur = con.cursor()
    cur.execute("DROP TABLE IF EXISTS webfront_entryannotation")
    cur.execute(
        """
        CREATE TABLE webfront_entryannotation
        (
            annotation_id INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
            accession VARCHAR(25) NOT NULL,
            type VARCHAR(20) NOT NULL,
            value LONGBLOB NOT NULL,
            mime_type VARCHAR(32) NOT NULL,
            num_sequences INT
        ) CHARSET=utf8mb4 DEFAULT COLLATE=utf8mb4_unicode_ci
        """
    )
    cur.close()
    con.close()

    queue = Queue()
    consumer = Process(target=_insert,
                       args=(stg_url, queue))
    consumer.start()

    dt = DirectoryTree(root=tmpdir)

    # Get HMMs from InterPro Oracle database
    for path in _export_hmms(p_uniprot2matches, pro_url, dt):
        queue.put(path)

    # Get alignments from Pfam MySQL database
    for path in _export_alns(pfam_url, dt):
        queue.put(path)

    queue.put(None)
    consumer.join()
    dt.remove()

    con = MySQLdb.connect(**url2dict(stg_url))
    cur = con.cursor()
    cur.execute("CREATE INDEX i_entryannotation "
                "ON webfront_entryannotation (accession)")
    cur.close()
    con.close()
