import glob
import os
import tarfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from xml.dom.minidom import getDOMImplementation

from interpro7dw.utils import logger
from interpro7dw.utils.store import BasicStore


_ARCHIVE = "uniparc_match.tar.gz"
_WITH_ALIGNMENT = {
    "HAMAP",
    "PROSITE patterns",
    "PROSITE profiles",
}


def archive_matches(indir: str, outdir: str, processes: int = 8):
    logger.info("Writing XML files")
    os.makedirs(outdir, exist_ok=True)

    errors = 0
    with (ProcessPoolExecutor(max_workers=max(1, processes - 1)) as executor,
          tarfile.open(os.path.join(outdir, _ARCHIVE), "w:gz") as fh):
        files = glob.glob(os.path.join(indir, "*.dat"))

        fs = {}
        for i, bspath in enumerate(sorted(files)):
            xmlpath = os.path.join(outdir, f"{str(i+1).zfill(6)}.xml")
            f = executor.submit(write_xml, bspath, xmlpath)
            fs[f] = xmlpath

        done = 0
        milestone = step = 5
        for f in as_completed(fs):
            try:
                f.result()
            except Exception as exc:
                logger.error(exc)
                errors += 1
            else:
                xmlpath = fs[f]
                fh.add(xmlpath, arcname=os.path.basename(xmlpath))

                os.unlink(xmlpath)

                done += 1
                progress = done * 100 / len(fs)
                if progress >= milestone:
                    logger.info(f"{progress:.0f}%")
                    milestone += step

    if errors:
        raise RuntimeError(f"{errors} errors")

    logger.info("done")


def write_xml(bspath: str, xmlpath: str):
    with BasicStore(bspath) as bs, open(xmlpath, "wt") as fh:
        fh.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        doc = getDOMImplementation().createDocument(None, None, None)
        for proteins in bs:
            for upi, protein in proteins.items():
                protein_elem = doc.createElement("protein")
                protein_elem.setAttribute("id", upi)
                protein_elem.setAttribute("length", str(protein["length"]))
                protein_elem.setAttribute("crc64", protein["crc64"])

                for match in protein["matches"]:
                    """PIRSR matches do not represent 'real' family/domain hits, so
                    these matches should not be exported alongside the other member db matches
                    """
                    if match["signature"]["signatureLibraryRelease"]["name"] == "PIRSR":
                        continue

                    signature = match["signature"]
                    database = signature["signatureLibraryRelease"]

                    match_elem = doc.createElement("match")
                    match_elem.setAttribute("id", signature["accession"])
                    match_elem.setAttribute("name", signature["name"])
                    match_elem.setAttribute("dbname", database["library"])
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

                        feature = loc["sequence-feature"]
                        if feature:
                            if database["library"] in _WITH_ALIGNMENT:
                                attname = "alignment"
                            else:
                                attname = "sequence-feature"

                            lcn_elem.setAttribute(attname, feature)

                        match_elem.appendChild(lcn_elem)

                    protein_elem.appendChild(match_elem)

                protein_elem.writexml(fh, addindent="  ", newl="\n")
