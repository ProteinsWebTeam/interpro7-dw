import gzip
import json
import os
from io import StringIO
from multiprocessing import Process, Queue

import MySQLdb

from interpro7dw import hmmer
from interpro7dw.ebi import pfam
from interpro7dw.ebi.interpro import production as ippro
from interpro7dw.utils import DumpFile, DirectoryTree
from interpro7dw.utils import url2dict


def _export_annotations(pro_url: str, pfam_url: str, dt: DirectoryTree,
                        queue: Queue, buffer_size: int = 1000):
    cnt = 0
    df = DumpFile(dt.mktemp(), compress=True)

    iterator = ippro.get_hmms(pro_url, multi_models=False)
    for signature_acc, model_acc, hmm_bytes in iterator:
        hmm_str = gzip.decompress(hmm_bytes).decode("utf-8")
        df.dump((
            signature_acc,
            "hmm",
            hmm_bytes,
            "application/gzip",
            None
        ))

        with StringIO(hmm_str) as stream:
            hmm = hmmer.HMMFile(stream)

        df.dump((
            signature_acc,
            "logo",
            json.dumps(hmm.logo("info_content_all", "hmm")),
            "application/json",
            None
        ))

        cnt += 2
        if cnt >= buffer_size:
            df.close()
            queue.put(df.path)
            df = DumpFile(dt.mktemp(), compress=True)
            cnt = 0

    iterator = pfam.get_alignments(pfam_url)
    for signature_acc, aln_type, aln_bytes, count in iterator:
        df.dump((
            signature_acc,
            f"alignment:{aln_type}",
            aln_bytes,
            "application/gzip",
            count
        ))

        cnt += 1
        if cnt >= buffer_size:
            df.close()
            queue.put(df.path)
            df = DumpFile(dt.mktemp(), compress=True)
            cnt = 0

    df.close()
    queue.put(df.path)
    queue.put(None)


def insert_annotations(pro_url: str, pfam_url: str, stg_url: str, **kwargs):
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

    dt = DirectoryTree(root=tmpdir)
    queue = Queue()
    producer = Process(target=_export_annotations,
                       args=(pro_url, pfam_url, dt, queue))
    producer.start()

    for path in iter(queue.get, None):
        with DumpFile(path) as df:
            con = MySQLdb.connect(**url2dict(stg_url))
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

    producer.join()
    dt.remove()

    con = MySQLdb.connect(**url2dict(stg_url))
    cur = con.cursor()
    cur.execute("CREATE INDEX i_entryannotation "
                "ON webfront_entryannotation (accession)")
    cur.close()
    con.close()
