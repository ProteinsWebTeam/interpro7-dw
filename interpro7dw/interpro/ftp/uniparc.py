import os
import tarfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from xml.dom.minidom import getDOMImplementation

from interpro7dw.utils import logger
from interpro7dw.utils.store import BasicStore


_ARCHIVE = "uniparc_match.tar.gz"


def write_xml(bsfile: str, xmlfile: str):
    with BasicStore(bsfile) as bs, open(xmlfile, "wt") as fh:
        fh.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        doc = getDOMImplementation().createDocument(None, None, None)
        for proteins in bs:
            for upi, protein in proteins.items():
                protein_elem = doc.createElement("protein")
                protein_elem.setAttribute("id", upi)
                protein_elem.setAttribute("length", str(protein["length"]))
                protein_elem.setAttribute("crc64", protein["crc64"])

                for match in protein["matches"]:
                    signature = match["signature"]

                    match_elem = doc.createElement("match")
                    match_elem.setAttribute("id", signature["accession"])
                    match_elem.setAttribute("name", signature["name"])
                    match_elem.setAttribute("dbname", signature["signatureLibraryRelease"]["name"])
                    match_elem.setAttribute("status", "T")
                    match_elem.setAttribute("evd", match["extra"]["evidence"])
                    match_elem.setAttribute("model", match["model-ac"])

                    if signature["entry"]:
                        entry = signature["entry"]

                        ipr_elem = doc.createElement("ipr")
                        ipr_elem.setAttribute("id", entry["accession"])
                        ipr_elem.setAttribute("name", entry["name"])
                        ipr_elem.setAttribute("type", entry["type"])

                        if entry["parent"]:
                            ipr_elem.setAttribute("parent_id", entry["parent"])

                        match_elem.appendChild(ipr_elem)

                    for loc in match["locations"]:
                        lcn_elem = doc.createElement("lcn")
                        lcn_elem.setAttribute("start", str(loc["start"]))
                        lcn_elem.setAttribute("end", str(loc["end"]))
                        lcn_elem.setAttribute("score", str(loc["score"]))

                        if loc["extra"]["fragments"]:
                            lcn_elem.setAttribute("fragments", loc["extra"]["fragments"])

                        # TODO: add this for HAMAP and PROSITE
                        # if loc["sequence-feature"]:
                        #     lcn_elem.setAttribute("alignment", loc["sequence-feature"])

                        match_elem.appendChild(lcn_elem)

                    protein_elem.appendChild(match)

                protein_elem.writexml(fh, addindent="  ", newl="\n")


def archive_matches(indir: str, outdir: str, processes: int = 8):
    logger.info("Writing XML files")
    os.makedirs(outdir, exist_ok=True)

    errors = 0
    with (ProcessPoolExecutor(max_workers=max(1, processes - 1)) as executor,
          tarfile.open(os.path.join(outdir, _ARCHIVE), "w:gz") as fh):
        fs = {}
        for i, filename in enumerate(sorted(os.listdir(indir))):
            bsfile = os.path.join(indir, filename)
            xmlfile = os.path.join(outdir, f"{str(i+1).zfill(5)}.xml")
            f = executor.submit(write_xml, bsfile, xmlfile)
            fs[f] = xmlfile

        done = 0
        milestone = step = 5
        for f in as_completed(fs):
            try:
                f.result()
            except Exception as exc:
                logger.error(exc)
                errors += 1
            else:
                xmlfile = fs[f]
                fh.add(xmlfile, arcname=os.path.basename(xmlfile))

                os.unlink(xmlfile)

                done += 1
                progress = done * 100 / len(fs)
                if progress >= milestone:
                    logger.info(f"{progress:.0f%}")
                    milestone += step

    if errors:
        raise RuntimeError(f"{errors} errors")

    logger.info("done")
