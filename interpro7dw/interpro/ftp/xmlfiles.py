import gzip
import math
import multiprocessing as mp
import os
import pickle
import re
import shutil
from typing import Optional
from xml.dom.minidom import getDOMImplementation, parseString
from xml.parsers.expat import ExpatError

from interpro7dw.interpro.oracle.matches import DC_STATUSES
from interpro7dw.utils import logger
from interpro7dw.utils.store import BasicStore, KVStore


_FEATURES_DTD = "extra.dtd"
_FEATURES_XML = "extra.xml.gz"
_INTERPRO_DTD = "interpro.dtd"
_INTERPRO_XML = "interpro.xml.gz"
_MATCHES_DTD = "match_complete.dtd"
_MATCHES_XML = "match_complete.xml.gz"
_STRUCTURES_DTD = "feature.dtd"
_STRUCTURES_XML = "feature.xml.gz"
_DC_STATUSES = {value: key for key, value in DC_STATUSES.items()}
_KEY_SPECIES = {
    "3702",  # Arabidopsis thaliana
    "6239",  # Caenorhabditis elegans
    "7955",  # Danio rerio
    "7227",  # Drosophila melanogaster
    "9606",  # Homo sapiens
    "10090",  # Mus musculus
    "367110",  # Neurospora crassa
    "10116",  # Rattus norvegicus
    "559292",  # Saccharomyces cerevisiae
    "284812",  # Schizosaccharomyces pombe
    "4577",  # Zea mays
}
_SUPERKINGDOMS = {
    "Archaea",
    "Bacteria",
    "Eukaryota",
    "Viruses"
}
_TAGS = {
    "cazy": "CAZY",
    "cog": "COG",
    "genprop": "GENPROP",
    "ec": "EC",
    "intenz": "EC",
    "interpro": "INTERPRO",
    "pfam": "PFAM",
    "pdbe": "PDBE",
    "pirsf": "PIRSF",
    "prosite": "PROSITE",
    "prositedoc": "PROSITEDOC",
    "superfamily": "SSF",
    "swissprot": "SWISSPROT",
    "tigrfams": "TIGRFAMs"
}


def _restore_tags(match: re.Match) -> str:
    tag, key = match.groups()
    tag = tag.lower()
    if tag == "cite":
        return f'<cite idref="{key}"/>'
    elif tag in _TAGS:
        return f'<db_xref db="{_TAGS[tag]}" dbkey="{key}"/>'
    elif tag not in ["mim", "pmid", "pubmed"]:
        logger.warning(match.group(0))


def _restore_abstract(data: str) -> str:
    return re.sub(pattern=r"\[([a-z]+):([a-z0-9_.:]+)\]",
                  repl=_restore_tags,
                  string=data,
                  flags=re.I)


def export_interpro(entries_file: str, entry2xrefs_file: str,
                    databases_file: str, taxa_file: str, outdir: str):
    os.makedirs(outdir, exist_ok=True)
    shutil.copy(os.path.join(os.path.dirname(__file__), _INTERPRO_DTD),
                outdir)

    logger.info("loading entries")
    with open(entries_file, "rb") as fh:
        entries = pickle.load(fh)

    public_entries = set()
    deleted_entries = set()
    entry2children = {}
    entry2signatures = {}
    for entry in entries.values():
        if entry.database.lower() != "interpro":
            if entry.integrated_in:
                if entry.integrated_in in entry2signatures:
                    entry2signatures[entry.integrated_in].append(entry)
                else:
                    entry2signatures[entry.integrated_in] = [entry]
        elif entry.deletion_date is None:
            public_entries.add(entry.accession)

            if entry.parent:
                if entry.parent in entry2children:
                    entry2children[entry.parent].append(entry.accession)
                else:
                    entry2children[entry.parent] = [entry.accession]
        else:
            deleted_entries.add(entry.accession)

    logger.info("loading taxa")
    with open(taxa_file, "rb") as fh:
        taxa = pickle.load(fh)

    entry2ancestors = {}
    entry2proteins = {}
    entry2species = {}
    entry2structures = {}

    logger.info("reading cross-references")
    with BasicStore(entry2xrefs_file, mode="r") as store:
        for entry_acc, entry_xrefs in store:
            entry2proteins[entry_acc] = len(entry_xrefs["proteins"])

            if entry_acc not in public_entries:
                continue

            entry = entries[entry_acc]
            if entry_xrefs["enzymes"]:
                entry.cross_references["EC"] = sorted(entry_xrefs["enzymes"])

            if entry_xrefs["metacyc"]:
                # item: tuple (pathway ID, pathway name)
                pathways = [item[0] for item in entry_xrefs["metacyc"]]
                entry.cross_references["METACYC"] = sorted(pathways)

            if entry_xrefs["reactome"]:
                # item: tuple (pathway ID, pathway name)
                pathways = [item[0] for item in entry_xrefs["reactome"]]
                entry.cross_references["REACTOME"] = sorted(pathways)

            if entry_xrefs["structures"]:
                entry2structures[entry_acc] = sorted(entry_xrefs["structures"])

            superkingdoms = {}
            entry_taxa = entry_xrefs["taxa"]["all"]
            for taxon_id, num_proteins in entry_taxa.items():
                lineage = taxa[taxon_id]["lineage"]
                superkingdom_id = lineage[0]

                try:
                    other_lineage = superkingdoms[superkingdom_id]
                except KeyError:
                    superkingdoms[superkingdom_id] = lineage
                else:
                    # Compare lineages and find lowest common ancestor
                    i = 0
                    while i < len(lineage) and i < len(other_lineage):
                        if lineage[i] != other_lineage[i]:
                            break
                        i += 1

                    # Path to the lowest common ancestor
                    superkingdoms[superkingdom_id] = lineage[:i]

            # Get lowest common ancestor for each represented superkingdom
            entry2ancestors[entry_acc] = []
            for lineage in superkingdoms.values():
                # Lowest common ancestor
                taxon_id = lineage[-1]
                taxon = taxa[taxon_id]
                num_proteins = entry_taxa[taxon_id]
                entry2ancestors[entry_acc].append((taxon["sci_name"],
                                                   num_proteins))

            entry2species[entry_acc] = []
            for taxon_id in _KEY_SPECIES:
                try:
                    num_proteins = entry_taxa[taxon_id]
                except KeyError:
                    continue
                else:
                    taxon = taxa[taxon_id]
                    entry2species[entry_acc].append((taxon["sci_name"],
                                                     num_proteins))

    file = os.path.join(outdir, _INTERPRO_XML)
    with gzip.open(file, "wt", encoding="utf-8") as fh:
        fh.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        fh.write('<!DOCTYPE interprodb SYSTEM "interpro.dtd">\n')
        fh.write("<interprodb>\n")

        doc = getDOMImplementation().createDocument(None, None, None)
        # writing <release> section (do not log progress, < 1 sec)
        elem = doc.createElement("release")
        with open(databases_file, "rb") as fh2:
            for key, info in pickle.load(fh2).items():
                if info["type"] in ("entry", "protein"):
                    release = info["release"]
                    version = release["version"]
                    date = release["date"].strftime("%d-%b-%y").upper()

                    dbinfo = doc.createElement("dbinfo")
                    dbinfo.setAttribute("version", version)
                    dbinfo.setAttribute("dbname", key)
                    dbinfo.setAttribute("entry_count", str(info["entries"]))
                    dbinfo.setAttribute("file_date", date)
                    elem.appendChild(dbinfo)

        elem.writexml(fh, addindent="  ", newl="\n")

        logger.info("writing entries")
        for entry_acc in sorted(public_entries):
            entry = entries[entry_acc]
            elem = doc.createElement("interpro")
            elem.setAttribute("id", entry.accession)
            elem.setAttribute("protein_count", str(entry2proteins[entry_acc]))
            elem.setAttribute("short_name", entry.short_name)
            elem.setAttribute("type", entry.type)

            name = doc.createElement("name")
            name.appendChild(doc.createTextNode(entry.name))
            elem.appendChild(name)

            text = _restore_abstract('\n'.join(entry.descriptions))
            try:
                _doc = parseString(f"<abstract>{text}</abstract>")
            except ExpatError as exc:
                # TODO: use CDATA section for all entries
                logger.warning(f"{entry_acc}: {exc}")
                # abstract = doc.createElement("abstract")
                # abstract.appendChild(doc.createCDATASection(text))
            else:
                abstract = _doc.documentElement
                elem.appendChild(abstract)

            if entry.go_terms:
                go_list = doc.createElement("class_list")

                for term in entry.go_terms:
                    go_elem = doc.createElement("classification")
                    go_elem.setAttribute("id", term["identifier"])
                    go_elem.setAttribute("class_type", "GO")

                    _elem = doc.createElement("category")
                    _elem.appendChild(
                        doc.createTextNode(term["category"]["name"])
                    )
                    go_elem.appendChild(_elem)

                    _elem = doc.createElement("description")
                    _elem.appendChild(
                        doc.createTextNode(term["name"])
                    )
                    go_elem.appendChild(_elem)

                    go_list.appendChild(go_elem)

                elem.appendChild(go_list)

            if entry.literature:
                pub_list = doc.createElement("pub_list")
                for pub_id in sorted(entry.literature):
                    pub = entry.literature[pub_id]

                    pub_elem = doc.createElement("publication")
                    pub_elem.setAttribute("id", pub_id)

                    _elem = doc.createElement("author_list")
                    if pub["authors"]:
                        _elem.appendChild(
                            doc.createTextNode(", ".join(pub['authors']))
                        )
                    else:
                        _elem.appendChild(doc.createTextNode("Unknown"))
                    pub_elem.appendChild(_elem)

                    if pub["title"]:
                        _elem = doc.createElement("title")
                        _elem.appendChild(
                            doc.createTextNode(pub["title"])
                        )
                        pub_elem.appendChild(_elem)

                    if pub["URL"]:
                        _elem = doc.createElement("url")
                        _elem.appendChild(doc.createTextNode(pub["URL"]))
                        pub_elem.appendChild(_elem)

                    _elem = doc.createElement("db_xref")
                    if pub["PMID"]:
                        _elem.setAttribute("db", "PUBMED")
                        _elem.setAttribute("dbkey", str(pub["PMID"]))
                    else:
                        _elem.setAttribute("db", "MEDLINE")
                        _elem.setAttribute("dbkey", "MEDLINE")
                    pub_elem.appendChild(_elem)

                    if pub["ISO_journal"]:
                        _elem = doc.createElement("journal")
                        _elem.appendChild(
                            doc.createTextNode(pub["ISO_journal"])
                        )
                        pub_elem.appendChild(_elem)

                    if pub["ISBN"]:
                        _elem = doc.createElement("book_title")
                        isbn = f"ISBN:{pub['ISBN']}"
                        _elem.appendChild(doc.createTextNode(isbn))
                        pub_elem.appendChild(_elem)

                    if pub["volume"] or pub["issue"] or pub["raw_pages"]:
                        _elem = doc.createElement("location")
                        if pub["volume"]:
                            _elem.setAttribute("volume", pub["volume"])

                        if pub["issue"]:
                            _elem.setAttribute("issue", pub["issue"])

                        if pub["raw_pages"]:
                            _elem.setAttribute("pages", pub["raw_pages"])

                        pub_elem.appendChild(_elem)

                    if pub["year"]:
                        _elem = doc.createElement("year")
                        _elem.appendChild(
                            doc.createTextNode(str(pub["year"]))
                        )
                        pub_elem.appendChild(_elem)

                    pub_list.appendChild(pub_elem)

                elem.appendChild(pub_list)

            if entry.parent:
                par_elem = doc.createElement("parent_list")
                _elem = doc.createElement("rel_ref")
                _elem.setAttribute("ipr_ref", entry.parent)
                par_elem.appendChild(_elem)
                elem.appendChild(par_elem)

            if entry_acc in entry2children:
                child_list = doc.createElement("child_list")
                for child_acc in sorted(entry2children[entry_acc]):
                    _elem = doc.createElement("rel_ref")
                    _elem.setAttribute("ipr_ref", child_acc)
                    child_list.appendChild(_elem)

                elem.appendChild(child_list)

            members = entry2signatures.get(entry_acc, [])
            mem_list = doc.createElement("member_list")
            for member in sorted(members, key=lambda x: x.accession):
                _elem = doc.createElement("db_xref")
                _elem.setAttribute("protein_count",
                                   str(num_proteins[member.accession]))
                _elem.setAttribute("db", member.database)
                _elem.setAttribute("dbkey", member.accession)
                _elem.setAttribute("name", member.short_name)
                mem_list.appendChild(_elem)

            elem.appendChild(mem_list)

            if entry.cross_references:
                xref_list = doc.createElement("external_doc_list")
                for ref_db in sorted(entry.cross_references):
                    for ref_id in sorted(entry.cross_references[ref_db]):
                        _elem = doc.createElement("db_xref")
                        _elem.setAttribute("db", ref_db)
                        _elem.setAttribute("dbkey", ref_id)
                        xref_list.appendChild(_elem)
                elem.appendChild(xref_list)

            if entry_acc in entry2structures:
                xref_list = doc.createElement("structure_db_links")
                for pdb_id in entry2structures[entry_acc]:
                    _elem = doc.createElement("db_xref")
                    _elem.setAttribute("db", "PDB")
                    _elem.setAttribute("dbkey", pdb_id)
                    xref_list.appendChild(_elem)
                elem.appendChild(xref_list)

            # Write taxonomic distribution
            tax_dist = doc.createElement("taxonomy_distribution")
            for name in sorted(entry2ancestors):
                num_proteins = entry2ancestors[entry_acc]
                _elem = doc.createElement("taxon_data")
                _elem.setAttribute("name", name)
                _elem.setAttribute("proteins_count", str(num_proteins))
                tax_dist.appendChild(_elem)
            elem.appendChild(tax_dist)

            if entry2species[entry_acc]:
                # Write key species
                key_spec = doc.createElement("key_species")
                for name, num_proteins in sorted(entry2species[entry_acc]):
                    _elem = doc.createElement("taxon_data")
                    _elem.setAttribute("name", name)
                    _elem.setAttribute("proteins_count", str(num_proteins))
                    key_spec.appendChild(_elem)
                elem.appendChild(key_spec)

            elem.writexml(fh, addindent="  ", newl="\n")

        if deleted_entries:
            block = doc.createElement("deleted_entries")
            for entry_acc in sorted(deleted_entries):
                elem = doc.createElement("del_ref")
                elem.setAttribute("id", entry_acc)
                block.appendChild(elem)

            block.writexml(fh, addindent="  ", newl="\n")

        fh.write("</interprodb>\n")

    logger.info("complete")


def _export_matches(proteins_file: str, matches_file: str,
                    protein2isoforms: dict, start: str, stop: Optional[str],
                    output: str):
    with open(output, "wt") as fh:
        with KVStore(proteins_file) as st1, KVStore(matches_file) as st2:
            doc = getDOMImplementation().createDocument(None, None, None)

            for protein_acc, protein in st1.range(start, stop):
                elem = doc.createElement("protein")
                elem.setAttribute("id", protein_acc)
                elem.setAttribute("name", protein["identifier"])
                elem.setAttribute("length", str(protein["length"]))
                elem.setAttribute("crc64", protein["crc64"])

                signatures, entries = st2.get(protein_acc, ({}, {}))
                for signature_acc in sorted(signatures):
                    signature = signatures[signature_acc]
                    entry_acc = signature["entry"]
                    entry = entries[entry_acc] if entry_acc else None
                    elem.appendChild(create_match(doc, signature, entry))

                elem.writexml(fh, addindent="  ", newl="\n")

                isoforms = protein2isoforms.get(protein_acc, [])
                for variant_acc, length, crc64, matches in isoforms:
                    elem = doc.createElement("protein")
                    elem.setAttribute("id", variant_acc)
                    elem.setAttribute("name", variant_acc)
                    elem.setAttribute("length", str(length))
                    elem.setAttribute("crc64", crc64)

                    signatures, entries = matches
                    for signature_acc in sorted(signatures):
                        signature = signatures[signature_acc]
                        entry_acc = signature["entry"]
                        entry = entries[entry_acc] if entry_acc else None
                        elem.appendChild(create_match(doc, signature, entry))

                    elem.writexml(fh, addindent="  ", newl="\n")


def export_matches(databases_file: str, isoforms_file: str,
                   proteins_file: str, matches_file: str, outdir: str,
                   processes: int = 8):
    logger.info("starting")
    os.makedirs(outdir, exist_ok=True)
    shutil.copy(os.path.join(os.path.dirname(__file__), _MATCHES_DTD),
                outdir)

    logger.info("loading isoforms")
    protein2isoforms = {}
    with BasicStore(isoforms_file, mode="r") as store:
        for isoform in store:
            protein_acc = isoform["protein"]
            try:
                isoforms = protein2isoforms[protein_acc]
            except KeyError:
                isoforms = protein2isoforms[protein_acc] = []
            finally:
                isoforms.append((
                    isoform["accession"],
                    isoform["length"],
                    isoform["crc64"],
                    isoform["matches"]
                ))

    # Sorting isoforms by accession (so XXXX-1 comes before XXXX-2)
    for isoforms in protein2isoforms.values():
        isoforms.sort(key=lambda x: x[0])

    logger.info("writing XML files")
    with KVStore(matches_file) as store:
        keys = store.get_keys()

    processes = max(1, processes - 1)
    chunksize = math.ceil(len(keys) / processes)
    output = os.path.join(outdir, _MATCHES_XML)
    workers = []
    for i in range(processes):
        start = keys[i * chunksize]
        try:
            stop = keys[(i + 1) * chunksize]
        except IndexError:
            stop = None

        tempfile = f"{output}.{i+1}"
        p = mp.Process(target=_export_matches,
                       args=(proteins_file, matches_file, protein2isoforms,
                             start, stop, tempfile))
        p.start()
        workers.append((p, tempfile))

    with gzip.open(output, "wt", encoding="utf-8") as fh:
        fh.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        fh.write('<!DOCTYPE interpromatch SYSTEM "match_complete.dtd">\n')
        fh.write('<interpromatch>\n')

        doc = getDOMImplementation().createDocument(None, None, None)
        elem = doc.createElement("release")
        with open(databases_file, "rb") as fh2:
            for key, info in pickle.load(fh2).items():
                if info["type"] == "entry":
                    release = info["release"]
                    version = release["version"]
                    date = release["date"].strftime("%d-%b-%y").upper()

                    dbinfo = doc.createElement("dbinfo")
                    dbinfo.setAttribute("dbname", key)
                    dbinfo.setAttribute("version", version)
                    dbinfo.setAttribute("entry_count", str(info["entries"]))
                    dbinfo.setAttribute("file_date", date)
                    elem.appendChild(dbinfo)

        elem.writexml(fh, addindent="  ", newl="\n")

        for i, (p, tempfile) in enumerate(workers):
            p.join()

            with open(tempfile, "rt", encoding="utf-8") as fh2:
                while (block := fh2.read(1024)) != '':
                    fh.write(block)

            os.unlink(tempfile)
            logger.info(f"{i + 1:>6} / {len(workers)}")

        fh.write('</interpromatch>\n')

    logger.info("done")


def create_match(doc, match: dict, entry: Optional[dict]):
    elem = doc.createElement("match")
    elem.setAttribute("id", match["accession"])
    elem.setAttribute("name", match["name"])
    elem.setAttribute("dbname", match["database"])
    elem.setAttribute("status", 'T')
    """
    The model is stored in locations, so we get the model
    from the first location for the match's 'model' attribute
    """
    locations = match["locations"]
    elem.setAttribute("model", locations[0]["model"])
    elem.setAttribute("evd", match["evidence"])

    if entry:
        ipr = doc.createElement("ipr")
        ipr.setAttribute("id", match["entry"])
        ipr.setAttribute("name", entry["name"])
        ipr.setAttribute("type", entry["type"])

        if entry["parent"]:
            ipr.setAttribute("parent_id", entry["parent"])

        elem.appendChild(ipr)

    for loc in locations:
        elem.appendChild(create_lcn(doc, loc))

    return match


def create_lcn(doc, location: dict):
    fragments = location["fragments"]

    """
    We do not have to original start/end match positions,
    so we use the leftmost/rightmost fragment positions.

    We also reconstruct the fragment string (START-END-STATUS)
    """
    fragments_obj = []
    start = fragments[0]["start"]
    end = 0

    for frag in fragments:
        if frag["end"] > end:
            end = frag["end"]

        status = _DC_STATUSES[frag["dc-status"]]
        fragments_obj.append(f"{frag['start']}-{frag['end']}-{status}")

    lcn = doc.createElement("lcn")
    lcn.setAttribute("start", str(start))
    lcn.setAttribute("end", str(end))
    lcn.setAttribute("fragments", ','.join(fragments_obj))
    lcn.setAttribute("score", str(location["score"]))
    return lcn


def export_feature_matches(databases_file: str, proteins_file: str,
                           features_file: str, outdir: str):
    os.makedirs(outdir, exist_ok=True)
    shutil.copy(os.path.join(os.path.dirname(__file__), _FEATURES_DTD),
                outdir)

    logger.info("starting")
    file = os.path.join(outdir, _FEATURES_XML)
    with gzip.open(file, "wt", encoding="utf-8") as fh:
        fh.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        fh.write('<!DOCTYPE interproextra SYSTEM "extra.dtd">\n')
        fh.write('<interproextra>\n')

        doc = getDOMImplementation().createDocument(None, None, None)
        elem = doc.createElement("release")
        with open(databases_file, "rb") as fh2:
            for key, info in pickle.load(fh2).items():
                if info["type"] in "feature":
                    dbinfo = doc.createElement("dbinfo")
                    dbinfo.setAttribute("dbname", key)
                    version = info["release"]["version"]
                    if version:
                        dbinfo.setAttribute("version", version)

                    elem.appendChild(dbinfo)

        elem.writexml(fh, addindent="  ", newl="\n")

        with KVStore(proteins_file) as st1, KVStore(features_file) as st2:
            for i, (protein_acc, features) in enumerate(st2.items()):
                protein = st1[protein_acc]
                elem = doc.createElement("protein")
                elem.setAttribute("id", protein_acc)
                elem.setAttribute("name", protein["identifier"])
                elem.setAttribute("length", str(protein["length"]))
                elem.setAttribute("crc64", protein["crc64"])

                for feature_acc in sorted(features):
                    feature = features[feature_acc]

                    match = doc.createElement("match")
                    match.setAttribute("id", feature_acc)
                    match.setAttribute("name", feature["name"])
                    match.setAttribute("dbname", feature["database"])
                    match.setAttribute("status", 'T')
                    match.setAttribute("model", feature_acc)
                    match.setAttribute("evd", feature["evidence"])

                    for loc in feature["locations"]:
                        pos_start, pos_end, seq_feature = loc

                        lcn = doc.createElement("lcn")
                        lcn.setAttribute("start", str(pos_start))
                        lcn.setAttribute("end", str(pos_end))

                        if seq_feature:
                            lcn.setAttribute("sequence-feature", seq_feature)

                        match.appendChild(lcn)

                    elem.appendChild(match)

                elem.writexml(fh, addindent="  ", newl="\n")

                if (i + 1) % 1e7 == 0:
                    logger.info(f"{i + 1:>15,}")

            logger.info(f"{i + 1:>15,}")

        fh.write('</interproextra>\n')

    logger.info("done")


def export_structure_matches(structures_file: str, proteins_file: str,
                             protein2structures_file: str, outdir: str):
    os.makedirs(outdir, exist_ok=True)
    shutil.copy(os.path.join(os.path.dirname(__file__), "feature.dtd"),
                outdir)

    logger.info("loading PDBe data")
    with open(structures_file, "rb") as fh:
        data = pickle.load(fh)

    protein2cath = data["cath"]
    protein2scop = data["scop"]
    del data

    with open(protein2structures_file, "rb") as fh:
        protein2structures = pickle.load(fh)

    logger.info("writing file")
    output = os.path.join(outdir, "feature.xml.gz")
    with gzip.open(output, "wt", encoding="utf-8") as fh:
        fh.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        fh.write('<!DOCTYPE interprofeature SYSTEM "feature.dtd">\n')
        fh.write('<interprofeature>\n')

        with KVStore(proteins_file) as proteins:
            doc = getDOMImplementation().createDocument(None, None, None)

            for protein_acc, protein in proteins.items():
                pdbe_entries = protein2structures.get(protein_acc, {})
                cath_entries = protein2cath.get(protein_acc, {})
                scop_entries = protein2scop.get(protein_acc, {})

                if pdbe_entries or cath_entries or scop_entries:
                    elem = doc.createElement("protein")
                    elem.setAttribute("id", protein_acc)
                    elem.setAttribute("name", protein["identifier"])
                    elem.setAttribute("length", str(protein["length"]))
                    elem.setAttribute("crc64", protein["crc64"])

                    for pdbe_id in sorted(pdbe_entries):
                        chains = pdbe_entries[pdbe_id]
                        for chain_id in sorted(chains):
                            domain = doc.createElement("domain")
                            domain.setAttribute("id", f"{pdbe_id}{chain_id}")
                            domain.setAttribute("dbname", "PDB")

                            for loc in chains[chain_id]:
                                start = loc["protein_start"]
                                end = loc["protein_end"]

                                coord = doc.createElement("coord")
                                coord.setAttribute("pdb", pdbe_id)
                                coord.setAttribute("chain", chain_id)
                                coord.setAttribute("start", str(start))
                                coord.setAttribute("end", str(end))
                                domain.appendChild(coord)

                            elem.appendChild(domain)

                    for domain_id in sorted(cath_entries):
                        entry = cath_entries[domain_id]

                        domain = doc.createElement("domain")
                        domain.setAttribute("id", domain_id)
                        domain.setAttribute("cfn", entry["superfamily"]["id"])
                        domain.setAttribute("dbname", "CATH")

                        for loc in entry["locations"]:
                            coord = doc.createElement("coord")
                            coord.setAttribute("pdb", entry["pdb_id"])
                            coord.setAttribute("chain", entry["chain"])
                            coord.setAttribute("start", str(loc["start"]))
                            coord.setAttribute("end", str(loc["end"]))
                            domain.appendChild(coord)

                        elem.appendChild(domain)

                    for domain_id in sorted(scop_entries):
                        entry = scop_entries[domain_id]

                        domain = doc.createElement("domain")
                        domain.setAttribute("id", domain_id)
                        domain.setAttribute("cfn", entry["superfamily"]["id"])
                        domain.setAttribute("dbname", "SCOP")

                        for loc in entry["locations"]:
                            coord = doc.createElement("coord")
                            coord.setAttribute("pdb", entry["pdb_id"])
                            coord.setAttribute("chain", entry["chain"])
                            coord.setAttribute("start", str(loc["start"]))
                            coord.setAttribute("end", str(loc["end"]))
                            domain.appendChild(coord)

                        elem.appendChild(domain)

                    elem.writexml(fh, addindent="  ", newl="\n")

        fh.write('</interprofeature>\n')

    logger.info("complete")
