# -*- coding: utf-8 -*-

import multiprocessing as mp
from multiprocessing import Process, Queue
import os
import tarfile
from tempfile import mkstemp
from typing import List, Optional, Sequence, Tuple
from xml.dom.minidom import getDOMImplementation

import cx_Oracle

from interpro7dw import logger
from interpro7dw.ebi.interpro import production as ippro
from interpro7dw.utils import KVdb, Store


def dump_proteins(proteins_file: str, matches_file: str, signatures: dict,
                   inqueue: mp.Queue, outqueue: mp.Queue):
    doc = getDOMImplementation().createDocument(None, None, None)
    with KVdb(proteins_file) as kvdb, Store(matches_file) as store:
        for from_upi, to_upi, filepath in iter(inqueue.get, None):
            with open(filepath, "wt") as fh:
                fh.write('<?xml version="1.0" encoding="UTF-8"?>\n')
                for upi, matches in store.range(from_upi, to_upi):
                    try:
                        length, crc64 = kvdb[upi]
                    except KeyError:
                        """
                        This may happen because UNIPARC.PROTEIN is refreshed 
                        using IPREAD while match data come from ISPRO, 
                        which uses UAPRO  (more up-to-date than UAREAD)
                        """
                        continue

                    protein = doc.createElement("protein")
                    protein.setAttribute("id", upi)
                    protein.setAttribute("length", str(length))
                    protein.setAttribute("crc64", crc64)

                    for signature_acc, model, locations in matches:
                        signature = signatures[signature_acc]

                        match = doc.createElement("match")
                        match.setAttribute("id", signature_acc)
                        match.setAttribute("name", signature["name"])
                        match.setAttribute("dbname", signature["database"])
                        match.setAttribute("status", 'T')
                        match.setAttribute("evd", signature["evidence"])
                        match.setAttribute("model", model)

                        if signature["interpro"]:
                            ipr = doc.createElement("ipr")
                            for attname, value in signature["interpro"]:
                                if value:
                                    ipr.setAttribute(attname, value)

                            match.appendChild(ipr)

                        for start, end, score, aln, frags in locations:
                            lcn = doc.createElement("lcn")
                            lcn.setAttribute("start", str(start))
                            lcn.setAttribute("end", str(end))

                            if frags:
                                lcn.setAttribute("fragments", frags)

                            if aln:
                                lcn.setAttribute("alignment", aln)

                            lcn.setAttribute("score", str(score))
                            match.appendChild(lcn)

                        protein.appendChild(match)

                    protein.writexml(fh, addindent="  ", newl="\n")

            outqueue.put(filepath)


def get_proteins(cur: cx_Oracle.Cursor, from_upi: str, to_upi: Optional[str]):
    if to_upi:
        sql = "WHERE M.UPI >= :1 AND M.UPI < :2"
        params = (from_upi, to_upi)
    else:
        sql = "WHERE M.UPI >= :1"
        params = (from_upi,)

    cur.execute(
        f"""
        SELECT M.UPI, P.LEN, P.CRC64, M.METHOD_AC, M.MODEL_AC, M.SEQ_START, 
               M.SEQ_END, M.SCORE, M.SEQ_FEATURE, M.FRAGMENTS
        FROM IPRSCAN.MV_IPRSCAN M
        INNER JOIN UNIPARC.PROTEIN P ON M.UPI = P.UPI
        {sql}
        ORDER BY UPI
        """, params
    )
    _upi = _length = _crc64 = None
    matches = []
    for upi, length, crc64, acc, model, start, stop, score, aln, frags in cur:
        if upi != _upi:
            if _upi:
                yield _upi, _length, _crc64, merge_matches(matches)

            _upi = upi
            _length = length,
            _crc64 = crc64
            matches = []

    if _upi:
        yield _upi, _length, _crc64, merge_matches(matches)


def dump_proteins2(url: str, inqueue: mp.Queue, outqueue: mp.Queue):
    doc = getDOMImplementation().createDocument(None, None, None)

    con = cx_Oracle.connect(url)
    cur = con.cursor()

    signatures = ippro.get_signatures(cur)

    for from_upi, to_upi, filepath in iter(inqueue.get, None):
        with open(filepath, "wt") as fh:
            proteins = get_proteins(cur, from_upi, to_upi)
            for upi, length, crc64, matches in proteins:
                protein = doc.createElement("protein")
                protein.setAttribute("id", upi)
                protein.setAttribute("length", str(length))
                protein.setAttribute("crc64", crc64)

                for signature_acc, model, locations in matches:
                    signature = signatures[signature_acc]

                    match = doc.createElement("match")
                    match.setAttribute("id", signature_acc)
                    match.setAttribute("name", signature["name"])
                    match.setAttribute("dbname", signature["database"])
                    match.setAttribute("status", 'T')
                    match.setAttribute("evd", signature["evidence"])
                    match.setAttribute("model", model)

                    if signature["interpro"]:
                        ipr = doc.createElement("ipr")
                        for attname, value in signature["interpro"]:
                            if value:
                                ipr.setAttribute(attname, value)

                        match.appendChild(ipr)

                    for start, end, score, aln, frags in locations:
                        lcn = doc.createElement("lcn")
                        lcn.setAttribute("start", str(start))
                        lcn.setAttribute("end", str(end))

                        if frags:
                            lcn.setAttribute("fragments", frags)

                        if aln:
                            lcn.setAttribute("alignment", aln)

                        lcn.setAttribute("score", str(score))
                        match.appendChild(lcn)

                    protein.appendChild(match)

                protein.writexml(fh, addindent="  ", newl="\n")

        outqueue.put(filepath)

    cur.close()
    con.close()


def merge_matches(matches: Sequence[dict]) -> List[Tuple[str, str, List]]:
    signatures = {}
    for acc, model, start, stop, score, aln, frags in matches:
        try:
            s = signatures[acc]
        except KeyError:
            s = signatures[acc] = {
                "model": model or acc,
                "locations": []
            }
        finally:
            s["locations"].append((start, stop, score, aln, frags))

    result = []
    for acc in sorted(signatures):
        s = signatures[acc]
        s["locations"].sort()
        result.append((acc, s["model"], s["locations"]))

    return result


def export_matches2(url: str, outdir: str, processes: int = 8,
                    proteins_per_file: int = 1000000):
    logger.info("spawning tasks")
    inqueue = Queue()
    outqueue = Queue()
    workers = []
    for _ in range(max(1, processes - 1)):
        p = Process(target=dump_proteins2,
                    args=(url, inqueue, outqueue))
        p.start()
        workers.append(p)

    num_files = 0
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    from_upi = None
    cur.execute("SELECT UPI FROM UNIPARC.PROTEIN ORDER BY UPI")
    for i, (upi, ) in enumerate(cur):
        if not i % proteins_per_file:
            if from_upi:
                num_files += 1
                filename = f"uniparc_match_{num_files}.dump"
                filepath = os.path.join(outdir, filename)
                inqueue.put((from_upi, upi, filepath))

            from_upi = upi

    num_files += 1
    filename = f"uniparc_match_{num_files}.dump"
    filepath = os.path.join(outdir, filename)
    inqueue.put((from_upi, None, filepath))

    cur.close()
    con.close()

    logger.info("creating XML archive")
    output = os.path.join(outdir, "uniparc_match.tar.gz")
    with tarfile.open(output, "w:gz") as fh:
        for i in range(num_files):
            filepath = outqueue.get()
            fh.add(filepath, arcname=os.path.basename(filepath))
            os.remove(filepath)
            logger.info(f"{i+1:>6}/{num_files}")

    for p in workers:
        p.join()

    logger.info("complete")


def export_matches(url: str, outdir: str, dir: Optional[str] = None,
                   processes: int = 8, proteins_per_file: int = 1000000):
    fd, proteins_file = mkstemp(dir=dir)
    os.close(fd)
    os.remove(proteins_file)

    logger.info("exporting UniParc proteins")
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    keys = []
    with KVdb(proteins_file, writeback=True) as kvdb:
        cur.execute(
            """
            SELECT UPI, LEN, CRC64
            FROM UNIPARC.PROTEIN
            ORDER BY UPI
            """
        )
        for i, (upi, length, crc64) in enumerate(cur):
            kvdb[upi] = (length, crc64)
            if not i % 1000000:
                kvdb.sync()

            if not i % 10000:
                keys.append(upi)

        kvdb.sync()

    logger.info("exporting UniParc matches")
    fd, matches_file = mkstemp(dir=dir)
    os.close(fd)
    with Store(matches_file, keys, dir) as store:
        cur.execute(
            """
            SELECT MA.UPI, MA.METHOD_AC, MA.MODEL_AC,
                   MA.SEQ_START, MA.SEQ_END, MA.SCORE, MA.SEQ_FEATURE,
                   MA.FRAGMENTS
            FROM IPRSCAN.MV_IPRSCAN MA
            INNER JOIN INTERPRO.METHOD ME
              ON MA.METHOD_AC = ME.METHOD_AC
            """
        )

        i = 0
        for row in cur:
            store.append(row[0], row[1:])

            i += 1
            if not i % 1000000:
                store.sync()

                if not i % 100000000:
                    logger.info(f"{i:>15,}")

        logger.info(f"{i:>15,}")
        size = store.merge(fn=merge_matches, processes=processes)

    logger.info("loading signatures")
    signatures = ippro.get_signatures(cur)
    cur.close()
    con.close()

    logger.info("spawning processes")
    ctx = mp.get_context(method="spawn")
    inqueue = ctx.Queue()
    outqueue = ctx.Queue()
    workers = []
    for _ in range(max(1, processes - 1)):
        p = ctx.Process(target=dump_proteins,
                        args=(proteins_file, matches_file, signatures,
                              inqueue, outqueue))
        p.start()
        workers.append(p)

    with Store(matches_file) as store:
        num_files = 0

        i = 0
        from_upi = None
        for upi in store:
            i += 1
            if not i % 10000000:
                logger.info(f"{i:>13,}")

            if i % proteins_per_file == 1:
                if from_upi:
                    num_files += 1
                    filename = f"uniparc_match_{num_files}.dump"
                    filepath = os.path.join(outdir, filename)
                    inqueue.put((from_upi, upi, filepath))

                from_upi = upi

        num_files += 1
        filename = f"uniparc_match_{num_files}.dump"
        filepath = os.path.join(outdir, filename)
        inqueue.put((from_upi, None, filepath))
        logger.info(f"{i:>13,}")

    for _ in workers:
        inqueue.put(None)

    logger.info("creating XML archive")
    output = os.path.join(outdir, "uniparc_match.tar.gz")
    with tarfile.open(output, "w:gz") as fh:
        for i in range(num_files):
            filepath = outqueue.get()
            fh.add(filepath, arcname=os.path.basename(filepath))
            os.remove(filepath)
            logger.info(f"{i+1:>6}/{num_files}")

    for p in workers:
        p.join()

    size += os.path.getsize(proteins_file)
    os.remove(proteins_file)
    size += os.path.getsize(matches_file)
    os.remove(matches_file)
    logger.info(f"temporary files: {size/1024/1024:.0f} MB")
    logger.info("complete")
