import math
import multiprocessing as mp
import os
import tarfile
from xml.dom.minidom import getDOMImplementation

from interpro7dw.utils import logger
from interpro7dw.utils.store import KVStore


_ARCHIVE = "uniparc_match.tar.gz"


def chunk(file: str, chunksize: int):
    with KVStore(file) as store:
        keys = store.get_keys()
        store_chunksize = 0
        for _ in store.range(keys[0], keys[1]):
            store_chunksize += 1

        if chunksize % store_chunksize == 0 and chunksize >= store_chunksize:
            start = None
            i = 0
            for key in keys:
                if i % chunksize == 0:
                    if start:
                        yield start, key

                    start = key

                i += store_chunksize

            yield start, None
        else:
            start = None
            for i, key in enumerate(store):
                if i % chunksize == 0:
                    if start:
                        yield start, key

                    start = key

            yield start, None


def write_xml(proteins_file: str, matches_file: str, inqeue: mp.Queue,
              outqueue: mp.Queue):
    with KVStore(proteins_file) as s1, KVStore(matches_file) as s2:
        for start, stop, output in iter(inqeue.get, None):
            count = 0
            with open(output, "wt") as fh:
                fh.write('<?xml version="1.0" encoding="UTF-8"?>\n')
                doc = getDOMImplementation().createDocument(None, None, None)
                for upi, matches in s2.range(start, stop):
                    try:
                        length, crc64, md5 = s1[upi]
                    except KeyError:
                        """
                        This may happen if matches are calculated 
                        against sequences in UAPRO instead of UAREAD
                        """
                        continue

                    protein = doc.createElement("protein")
                    protein.setAttribute("id", upi)
                    protein.setAttribute("length", str(length))
                    protein.setAttribute("crc64", crc64)

                    for signature_acc in sorted(matches):
                        signature = matches[signature_acc]

                        match = doc.createElement("match")
                        match.setAttribute("id", signature_acc)
                        match.setAttribute("name", signature["name"])
                        match.setAttribute("dbname", signature["database"])
                        match.setAttribute("status", 'T')
                        match.setAttribute("evd", signature["evidence"])
                        match.setAttribute("model", signature["model"])

                        if signature["entry"]:
                            entry = signature["entry"]

                            ipr = doc.createElement("ipr")
                            ipr.setAttribute("id", entry["accession"])
                            ipr.setAttribute("name", entry["name"])
                            ipr.setAttribute("type", entry["type"])

                            if entry["parent"]:
                                ipr.setAttribute("parent_id", entry["parent"])

                            match.appendChild(ipr)

                        for loc in signature["locations"]:
                            pos_start, pos_end, score, aln, frags = loc

                            lcn = doc.createElement("lcn")
                            lcn.setAttribute("start", str(pos_start))
                            lcn.setAttribute("end", str(pos_end))
                            lcn.setAttribute("score", str(score))

                            if frags:
                                lcn.setAttribute("fragments", frags)

                            if aln:
                                lcn.setAttribute("alignment", aln)

                            match.appendChild(lcn)

                        protein.appendChild(match)

                    protein.writexml(fh, addindent="  ", newl="\n")
                    count += 1

            outqueue.put((output, count))


def archive_matches(proteins_file: str, matches_file: str, outdir: str,
                    processes: int = 8, proteins_per_file: int = 1000000):
    logger.info("Writing XML files")
    os.makedirs(outdir, exist_ok=True)

    inqueue = mp.Queue()
    outqueue = mp.Queue()
    workers = []
    for _ in range(max(1, processes - 1)):
        p = mp.Process(target=write_xml,
                       args=(proteins_file, matches_file, inqueue, outqueue))
        p.start()
        workers.append(p)

    num_files = 0
    for start, stop in chunk(matches_file, proteins_per_file):
        num_files += 1
        filename = f"uniparc_match_{num_files}.dump"
        filepath = os.path.join(outdir, filename)
        inqueue.put((start, stop, filepath))

    for _ in workers:
        inqueue.put(None)

    logger.info("Archiving XML files")
    with tarfile.open(os.path.join(outdir, _ARCHIVE), "w:gz") as fh:
        progress = 0
        milestone = step = math.ceil(0.1 * num_files)
        for _ in range(num_files):
            filepath, num_proteins = outqueue.get()
            if num_proteins > 0:
                fh.add(filepath, arcname=os.path.basename(filepath))

            os.unlink(filepath)
            progress += 1
            if progress == milestone:
                logger.info(f"{progress:>15,.0f} / {num_files:,}")
                milestone += step

        logger.info(f"{progress:>15,.0f} / {num_files:,}")

    for p in workers:
        p.join()

    logger.info("done")
