import gzip
import math
import multiprocessing as mp
import os
import pickle
import re
import shutil

from xml.dom.minidom import getDOMImplementation, parseString
from xml.parsers.expat import ExpatError

from interpro7dw.interpro.oracle.matches import DC_STATUSES
from interpro7dw.utils import logger
from interpro7dw.utils.store import BasicStore, KVStore

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
    "83333",  # Escherichia coli
    "9606",  # Homo sapiens
    "10090",  # Mus musculus
    "367110",  # Neurospora crassa
    "39947",  # Oryza sativa subsp. japonica
    "10116",  # Rattus norvegicus
    "559292",  # Saccharomyces cerevisiae
    "284812",  # Schizosaccharomyces pombe
    "4577",  # Zea mays
}
_TAGS = {
    "cazy": "CAZY",
    "cog": "COG",
    "genprop": "GENPROP",
    "ec": "EC",
    "intenz": "EC",
    "interpro": "INTERPRO",
    "ncbifam": "NCBIFAM",
    "pfam": "PFAM",
    "pdbe": "PDBE",
    "pirsf": "PIRSF",
    "prosite": "PROSITE",
    "prositedoc": "PROSITEDOC",
    "superfamily": "SSF",
    "swissprot": "SWISSPROT",
}


def _restore_tags(match: re.Match) -> str:
    tag, key = match.groups()
    tag = tag.lower()
    if tag == "cite":
        return f'<cite idref="{key}"/>'
    elif tag in _TAGS:
        return f'<db_xref db="{_TAGS[tag]}" dbkey="{key}"/>'
    elif tag not in ["omim", "pmid", "pubmed"]:
        logger.warning(match.group(0))


def _restore_abstract(data: str) -> str:
    return re.sub(pattern=r"\[([a-z]+):([a-z0-9_.:]+)\]",
                  repl=_restore_tags,
                  string=data,
                  flags=re.I)


def export_interpro(
    entries_file: str,
    entry2xrefs_file: str,
    databases_file: str,
    taxa_file: str,
    outdir: str,
):
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
        if entry.database.lower() == "interpro":
            if entry.deletion_date is None:
                public_entries.add(entry.accession)

                if entry.parent:
                    if entry.parent in entry2children:
                        entry2children[entry.parent].append(entry.accession)
                    else:
                        entry2children[entry.parent] = [entry.accession]
            else:
                deleted_entries.add(entry.accession)
        elif entry.integrated_in:
            if entry.integrated_in in entry2signatures:
                entry2signatures[entry.integrated_in].append(entry)
            else:
                entry2signatures[entry.integrated_in] = [entry]

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
                entry2structures[entry_acc] = sorted(
                    [pdb_id for pdb_id, _ in entry_xrefs["structures"]]
                )

            superkingdoms = {}
            entry_taxa = entry_xrefs["taxa"]
            for taxon_id, num_proteins in entry_taxa["hit"].items():
                lineage = taxa[taxon_id]["lineage"]

                i = 0
                superkingdom_id = None
                for i, node_id in enumerate(lineage):
                    if node_id not in ("1", "131567"):
                        """
                        Skip root (1) and meta-superkingdom (131567) containing:
                            * Bacteria (2)
                            * Archaea (2157)
                            * Eukaryota (2759)
                        """
                        superkingdom_id = node_id
                        break

                if superkingdom_id is None:
                    raise ValueError(f"no superkingdom found "
                                     f"(entry {entry_acc}, taxon {taxon_id})")

                lineage = lineage[i:]

                try:
                    other_lineage = superkingdoms[superkingdom_id]
                except KeyError:
                    superkingdoms[superkingdom_id] = lineage
                else:
                    # Compare lineages and find the lowest common ancestor
                    i = 0
                    while i < len(lineage) and i < len(other_lineage):
                        if lineage[i] != other_lineage[i]:
                            break
                        i += 1

                    # Path to the lowest common ancestor
                    superkingdoms[superkingdom_id] = lineage[:i]

            # Get the lowest common ancestor for each represented superkingdom
            entry2ancestors[entry_acc] = []
            for lineage in superkingdoms.values():
                # Lowest common ancestor
                taxon_id = lineage[-1]
                taxon = taxa[taxon_id]
                num_proteins = entry_taxa["all"][taxon_id]
                entry2ancestors[entry_acc].append((taxon["sci_name"],
                                                   num_proteins))

            entry2species[entry_acc] = []
            for taxon_id in _KEY_SPECIES:
                try:
                    num_proteins = entry_taxa["hit"][taxon_id]
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
            elem.setAttribute("protein_count",
                              str(entry2proteins.get(entry_acc, 0)))
            elem.setAttribute("short_name", entry.short_name)
            elem.setAttribute("type", entry.type)

            elem.setAttribute("is-llm", "true" if entry.llm else "false")
            elem.setAttribute("is-llm-reviewed", "true" if entry.llm_reviewed else "false")

            name = doc.createElement("name")
            name.appendChild(doc.createTextNode(entry.name))
            elem.appendChild(name)

            # label abstract (ab) is ai-generated, and if reviewed
            llm_descrs = reviewed_llm_descrs = 0
            blocks = []
            for descr in entry.descriptions:
                blocks.append(descr["text"])
                if descr["llm"]:
                    llm_descrs += 1
                    if descr["checked"]:
                        reviewed_llm_descrs += 1

            text = _restore_abstract("\n".join(blocks))
            ab_is_llm = ab_is_reviewed_llm = "false"
            if llm_descrs > 0:
                # At least one AI-generated description
                ab_is_llm = "true"
                if reviewed_llm_descrs == llm_descrs:
                    # Considered reviewed if all AI descriptions are reviewed
                    ab_is_reviewed_llm = "true"

            try:
                _doc = parseString(f"<abstract>{text}</abstract>")
            except ExpatError as exc:
                # TODO: use CDATA section for all entries
                logger.warning(f"{entry_acc}: {exc} -- -*-")
                # abstract = doc.createElement("abstract")
                # abstract.appendChild(doc.createCDATASection(text))
            else:
                abstract = _doc.documentElement
                abstract.setAttribute("is-llm", ab_is_llm)
                abstract.setAttribute("is-llm-reviewed", ab_is_reviewed_llm)
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
            for mem in sorted(members, key=lambda x: x.accession):
                _elem = doc.createElement("db_xref")
                _elem.setAttribute("protein_count",
                                   str(entry2proteins.get(mem.accession, 0)))
                _elem.setAttribute("db", mem.database)
                _elem.setAttribute("dbkey", mem.accession)
                _elem.setAttribute("name", mem.short_name)
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
            ancestors = entry2ancestors.get(entry_acc, [])
            for name, num_proteins in sorted(ancestors):
                _elem = doc.createElement("taxon_data")
                _elem.setAttribute("name", name)
                _elem.setAttribute("proteins_count", str(num_proteins))
                tax_dist.appendChild(_elem)
            elem.appendChild(tax_dist)

            key_species = entry2species.get(entry_acc, [])
            if key_species:
                # Write key species
                key_spec = doc.createElement("key_species")
                for name, num_proteins in sorted(key_species):
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


def _export_matches(proteins_file: str, matches_file: str, features_file: str, 
                    protein2isoforms: dict, start: str, stop: str | None,
                    output: str):
    
    with open(output, "wt") as fh:

        with KVStore(proteins_file) as st1, \
            KVStore(matches_file) as st2, \
            KVStore(features_file) as ff:

            doc = getDOMImplementation().createDocument(None, None, None)
            for protein_acc, protein in st1.range(start, stop):
                elem = doc.createElement("protein")
                elem.setAttribute("id", protein_acc)
                elem.setAttribute("name", protein["identifier"])
                elem.setAttribute("length", str(protein["length"]))
                elem.setAttribute("crc64", protein["crc64"])
                elem.setAttribute("taxid", protein["taxid"])
                
                status = "reviewed" if protein["reviewed"] else "unreviewed"
                elem.setAttribute("status", status)

                signatures, entries = st2.get(protein_acc, ({}, {}))
                for signature_acc in sorted(signatures):
                    signature = signatures[signature_acc]

                    if signature["database"].lower() == "antifam":
                        # Ignore AntiFam families
                        continue

                    entry_acc = signature["entry"]
                    entry = entries[entry_acc] if entry_acc else None
                    for match in create_matches(doc, signature_acc, signature,
                                                entry):
                        elem.appendChild(match)

                # Add extra features
                features = ff.get(protein_acc, [])
                for feature in features:
                    extra_match = create_extra_match(doc, feature)
                    elem.appendChild(extra_match)
                    
                elem.writexml(fh, addindent="  ", newl="\n")

                isoforms = protein2isoforms.get(protein_acc, [])
                for variant_acc, length, crc64, matches in isoforms:
                    elem = doc.createElement("protein")
                    elem.setAttribute("id", variant_acc)
                    elem.setAttribute("name", variant_acc)
                    elem.setAttribute("length", str(length))
                    elem.setAttribute("crc64", crc64)
                    elem.setAttribute("status", status)
                    elem.setAttribute("taxid", protein["taxid"])

                    signatures, entries = matches
                    for signature_acc in sorted(signatures):
                        signature = signatures[signature_acc]
                        entry_acc = signature["entry"]
                        entry = entries[entry_acc] if entry_acc else None
                        for match in create_matches(doc, signature_acc,
                                                    signature, entry):
                            elem.appendChild(match)

                    elem.writexml(fh, addindent="  ", newl="\n")


def export_matches(databases_file: str, isoforms_file: str,
                   proteins_file: str, features_file: str, 
                   matches_file: str, outdir: str,
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
                       args=(proteins_file, matches_file, 
                             features_file, protein2isoforms,
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


def create_matches(doc, match_acc: str, match: dict, entry: dict | None):
    models = {}

    for location in match["locations"]:
        model_acc = location["model"]
        try:
            models[model_acc].append(location)
        except KeyError:
            models[model_acc] = [location]

    for model_acc, locations in models.items():
        elem = doc.createElement("match")
        elem.setAttribute("id", match_acc)
        elem.setAttribute("name", match["name"])
        elem.setAttribute("dbname", match["database"])
        elem.setAttribute("status", 'T')
        elem.setAttribute("model", model_acc)
        elem.setAttribute("evd", match["evidence"])
        elem.setAttribute("type", match["type"])

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

        yield elem


def create_extra_match(doc, feature):

    match = doc.createElement("match")
    match.setAttribute("id", feature["accession"])
    match.setAttribute("name", feature["name"])
    match.setAttribute("dbname", feature["database"])
    match.setAttribute("status", 'T')
    match.setAttribute("model", feature["accession"])
    match.setAttribute("evd", feature["evidence"])

    for loc in feature["locations"]:
        pos_start, pos_end, seq_feature = loc

        lcn = doc.createElement("lcn")
        lcn.setAttribute("start", str(pos_start))
        lcn.setAttribute("end", str(pos_end))

        if seq_feature:
            lcn.setAttribute("sequence-feature", seq_feature)
            
        match.appendChild(lcn)

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
    if location.get("representative"):
        lcn.setAttribute("representative", "true")
    else:
        lcn.setAttribute("representative", "false")

    return lcn

# def export_structure_matches(structures_file: str, proteins_file: str,
#                              protein2structures_file: str, outdir: str):
#     os.makedirs(outdir, exist_ok=True)
#     shutil.copy(os.path.join(os.path.dirname(__file__), _STRUCTURES_DTD),
#                 outdir)
#
#     logger.info("loading PDBe data")
#     with open(structures_file, "rb") as fh:
#         data = pickle.load(fh)
#
#     protein2cath = data["cath"]
#     protein2scop = data["scop"]
#     del data
#
#     with open(protein2structures_file, "rb") as fh:
#         protein2structures = pickle.load(fh)
#
#     logger.info("writing file")
#     output = os.path.join(outdir, _STRUCTURES_XML)
#     with gzip.open(output, "wt", encoding="utf-8") as fh:
#         fh.write('<?xml version="1.0" encoding="UTF-8"?>\n')
#         fh.write('<!DOCTYPE interprofeature SYSTEM "feature.dtd">\n')
#         fh.write('<interprofeature>\n')
#
#         with KVStore(proteins_file) as proteins:
#             doc = getDOMImplementation().createDocument(None, None, None)
#
#             for protein_acc, protein in proteins.items():
#                 pdbe_entries = protein2structures.get(protein_acc, {})
#                 cath_entries = protein2cath.get(protein_acc, {})
#                 scop_entries = protein2scop.get(protein_acc, {})
#
#                 if pdbe_entries or cath_entries or scop_entries:
#                     elem = doc.createElement("protein")
#                     elem.setAttribute("id", protein_acc)
#                     elem.setAttribute("name", protein["identifier"])
#                     elem.setAttribute("length", str(protein["length"]))
#                     elem.setAttribute("crc64", protein["crc64"])
#
#                     for pdbe_id in sorted(pdbe_entries):
#                         chains = pdbe_entries[pdbe_id]
#                         for chain_id in sorted(chains):
#                             domain = doc.createElement("domain")
#                             domain.setAttribute("id", f"{pdbe_id}{chain_id}")
#                             domain.setAttribute("dbname", "PDB")
#
#                             for loc in chains[chain_id]:
#                                 start = loc["protein_start"]
#                                 end = loc["protein_end"]
#
#                                 coord = doc.createElement("coord")
#                                 coord.setAttribute("pdb", pdbe_id)
#                                 coord.setAttribute("chain", chain_id)
#                                 coord.setAttribute("start", str(start))
#                                 coord.setAttribute("end", str(end))
#                                 domain.appendChild(coord)
#
#                             elem.appendChild(domain)
#
#                     for domain_id in sorted(cath_entries):
#                         entry = cath_entries[domain_id]
#
#                         domain = doc.createElement("domain")
#                         domain.setAttribute("id", domain_id)
#                         domain.setAttribute("cfn", entry["superfamily"]["id"])
#                         domain.setAttribute("dbname", "CATH")
#
#                         for loc in entry["locations"]:
#                             coord = doc.createElement("coord")
#                             coord.setAttribute("pdb", entry["pdb_id"])
#                             coord.setAttribute("chain", entry["chain"])
#                             coord.setAttribute("start", str(loc["start"]))
#                             coord.setAttribute("end", str(loc["end"]))
#                             domain.appendChild(coord)
#
#                         elem.appendChild(domain)
#
#                     for domain_id in sorted(scop_entries):
#                         entry = scop_entries[domain_id]
#
#                         domain = doc.createElement("domain")
#                         domain.setAttribute("id", domain_id)
#                         domain.setAttribute("cfn", entry["superfamily"]["id"])
#                         domain.setAttribute("dbname", "SCOP")
#
#                         for loc in entry["locations"]:
#                             coord = doc.createElement("coord")
#                             coord.setAttribute("pdb", entry["pdb_id"])
#                             coord.setAttribute("chain", entry["chain"])
#                             coord.setAttribute("start", str(loc["start"]))
#                             coord.setAttribute("end", str(loc["end"]))
#                             domain.appendChild(coord)
#
#                         elem.appendChild(domain)
#
#                     elem.writexml(fh, addindent="  ", newl="\n")
#
#         fh.write('</interprofeature>\n')
#
#     logger.info("complete")
