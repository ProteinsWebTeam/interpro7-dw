import os
import tarfile
from xml.dom.minidom import getDOMImplementation

from interpro7dw.utils import logger
from interpro7dw.utils.store import SimpleStore


_ARCHIVE = "uniparc_match.tar.gz"


def archive_uniparc_matches(matches_file: str, outdir: str,
                            proteins_per_file: int = 1000000):
    logger.info("writing proteins to XML files")

    os.makedirs(outdir, exist_ok=True)
    files = []
    filename = filepath = fh = None

    doc = getDOMImplementation().createDocument(None, None, None)
    with SimpleStore(matches_file) as store:
        for i, (upi, length, crc64, matches) in enumerate(store):
            if i % proteins_per_file == 0:
                filename = f"uniparc_match_{len(files) + 1}.dump"
                filepath = os.path.join(outdir, filename)
                files.append(filepath)
                fh = open(filepath, "wt")
                fh.write('<?xml version="1.0" encoding="UTF-8"?>\n')

            for signature in matches:
                accession = signature[0]
                name = signature[1]
                database = signature[2]
                evidence = signature[3]
                model = signature[4]
                entry = signature[5]
                locations = signature[6]

                protein = doc.createElement("protein")
                protein.setAttribute("id", upi)
                protein.setAttribute("length", str(length))
                protein.setAttribute("crc64", crc64)

                match = doc.createElement("match")
                match.setAttribute("id", accession)
                match.setAttribute("name", name)
                match.setAttribute("dbname", database)
                match.setAttribute("status", 'T')
                match.setAttribute("evd", evidence)
                match.setAttribute("model", model)

                if entry:
                    ipr = doc.createElement("ipr")
                    ipr.setAttribute("id", entry["id"])
                    ipr.setAttribute("name", entry["name"])
                    ipr.setAttribute("type", entry["type"])

                    if entry["parent_id"]:
                        ipr.setAttribute("parent_id", entry["parent_id"])

                    match.appendChild(ipr)

                for start, end, score, seq_feature, fragments in locations:
                    lcn = doc.createElement("lcn")
                    lcn.setAttribute("start", str(start))
                    lcn.setAttribute("end", str(end))

                    if fragments:
                        lcn.setAttribute("fragments", fragments)

                    if seq_feature:
                        lcn.setAttribute("alignment", seq_feature)

                    lcn.setAttribute("score", str(score))
                    match.appendChild(lcn)

                protein.appendChild(match)

            protein.writexml(fh, addindent="  ", newl="\n")

            if (i + 1) % 100e6 == 0:
                logger.info(f"{i + 1:>15,}")

        logger.info(f"{i + 1:>15,}")

    logger.info("creating XML archive")
    with tarfile.open(os.path.join(outdir, _ARCHIVE), "w:gz") as fh:
        for i, filepath in enumerate(files):
            logger.info(f"{i + 1:>6} / {len(files)}")
            fh.add(filepath, arcname=os.path.basename(filepath))
            os.remove(filepath)

    logger.info("done")
