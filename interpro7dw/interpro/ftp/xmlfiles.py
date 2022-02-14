import gzip
import os
import re
import shutil
from typing import Callable, List, Sequence
from xml.dom.minidom import getDOMImplementation, parseString
from xml.parsers.expat import ExpatError

from interpro7dw.interpro.oracle.proteins import DC_STATUSES
from interpro7dw.utils import logger
from interpro7dw.utils.store import SimpleStore, Store, loadobj


_FEATURES_DTD = "extra.dtd"
_FEATURES_XML = "extra.xml.gz"
_INTERPRO_DTD = "interpro.dtd"
_INTERPRO_XML = "interpro.xml.gz"
_MATCHES_DTD = "match_complete.dtd"
_MATCHES_XML = "match_complete.xml.gz"
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


def _find_nodes(node: dict, fn: Callable) -> List[dict]:
    nodes = []

    if fn(node):
        nodes.append(node)

    for child in node["children"]:
        nodes += _find_nodes(child, fn)

    return nodes


def _get_lineage(node: dict, node_id: str) -> List[str]:
    if node["id"] == node_id:
        return [node["id"]]

    for child in node["children"]:
        lineage = _get_lineage(child, node_id)
        if lineage:
            # Searched taxon was found in child
            if node["name"] != "root":
                lineage.append(node["id"])
            return lineage

    return []


def export_interpro(entries_file: str, entry2xrefs_file: str,
                    databases_file: str, outdir: str):
    logger.info("starting")
    shutil.copy(os.path.join(os.path.dirname(__file__), _INTERPRO_DTD),
                outdir)

    entries = loadobj(entries_file)
    public_entries = set()
    deleted_entries = set()
    for entry in entries.values():
        if entry.database != "interpro":
            continue
        elif entry.is_public:
            public_entries.add(entry.accession)
        else:
            deleted_entries.add(entry.accession)

    entry2species = {}
    entry2ancestors = {}

    logger.info("reading cross-references")
    with SimpleStore(entry2xrefs_file) as store:
        for entry_acc, entry_xrefs in store:
            if entry_acc not in public_entries:
                continue

            superkingdoms = {}
            tree = entry_xrefs["taxa"]["tree"]
            for taxon_id in entry_xrefs["taxa"]["all"]:
                lineage = _get_lineage(tree, taxon_id)[::-1]

                try:
                    superkingdom_id = lineage[0]
                except IndexError:
                    continue

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

                # Node of the LCA
                node = _find_nodes(tree, lambda x: x["id"] == taxon_id)[0]
                entry2ancestors[entry_acc].append((node["name"],
                                                   node["proteins"]))

            entry2species[entry_acc] = []
            for n in _find_nodes(tree, lambda x: x["id"] in _KEY_SPECIES):
                entry2species[entry_acc].append((n["name"], n["proteins"]))

    file = os.path.join(outdir, _INTERPRO_XML)
    with gzip.open(file, "wt", encoding="utf-8") as fh:
        fh.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        fh.write('<!DOCTYPE interprodb SYSTEM "interpro.dtd">\n')
        fh.write("<interprodb>\n")

        doc = getDOMImplementation().createDocument(None, None, None)
        # writing <release> section (do not log progress, < 1 sec)
        elem = doc.createElement("release")
        databases = {}
        for db in loadobj(databases_file):
            db_name = db[0]
            db_altname = db[1]
            db_type = db[4]
            db_num_entries = db[5]
            db_version = db[6]
            db_date = db[7]

            databases[db_name] = db_altname  # e.g. pfam -> PFAM

            if db_type in ("entry", "protein"):
                dbinfo = doc.createElement("dbinfo")
                dbinfo.setAttribute("version", db_version)
                dbinfo.setAttribute("dbname", db_altname)
                dbinfo.setAttribute("entry_count", str(db_num_entries))
                dbinfo.setAttribute("file_date",
                                    db_date.strftime("%d-%b-%y").upper())
                elem.appendChild(dbinfo)

        elem.writexml(fh, addindent="  ", newl="\n")

        logger.info("writing entries")
        for entry_acc in sorted(public_entries):
            entry = entries[entry_acc]
            elem = doc.createElement("interpro")
            elem.setAttribute("id", entry.accession)
            elem.setAttribute("protein_count", entry.counts["proteins"])
            elem.setAttribute("short_name", entry.short_name)
            elem.setAttribute("type", entry.type)

            name = doc.createElement("name")
            name.appendChild(doc.createTextNode(entry.name))
            elem.appendChild(name)

            text = _restore_abstract('\n'.join(entry.description))
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

            parent, children = entry.relations
            if parent:
                par_elem = doc.createElement("parent_list")
                _elem = doc.createElement("rel_ref")
                _elem.setAttribute("ipr_ref", parent)
                par_elem.appendChild(_elem)
                elem.appendChild(par_elem)

            if children:
                child_list = doc.createElement("child_list")
                for child in children:
                    _elem = doc.createElement("rel_ref")
                    _elem.setAttribute("ipr_ref", child)
                    child_list.appendChild(_elem)

                elem.appendChild(child_list)

            members = []
            for database, signatures in entry.integrates.items():
                for signature_acc in signatures:
                    members.append((
                        signature_acc,
                        entries[signature_acc].short_name,
                        database,
                        entries[signature_acc].counts["proteins"]
                    ))

            mem_list = doc.createElement("member_list")
            for member in sorted(members):
                _elem = doc.createElement("db_xref")
                _elem.setAttribute("protein_count", member[3])
                _elem.setAttribute("db", databases[member[2]])
                _elem.setAttribute("dbkey", member[0])
                _elem.setAttribute("name", member[1])
                mem_list.appendChild(_elem)
            elem.appendChild(mem_list)

            # Merge cross-references and pathways
            cross_refs = {}
            for key, values in entry.cross_references.items():
                cross_refs[databases[key]] = values

            for key, values in entry.pathways.items():
                cross_refs[databases[key]] = [val["id"] for val in values]

            if cross_refs:
                xref_list = doc.createElement("external_doc_list")
                for ref_db in sorted(cross_refs):
                    for ref_id in sorted(cross_refs[ref_db]):
                        _elem = doc.createElement("db_xref")
                        _elem.setAttribute("db", ref_db)
                        _elem.setAttribute("dbkey", ref_id)
                        xref_list.appendChild(_elem)
                elem.appendChild(xref_list)

            if entry_xrefs["structures"]:
                xref_list = doc.createElement("structure_db_links")
                for pdbe_id in sorted(entry_xrefs["structures"]):
                    _elem = doc.createElement("db_xref")
                    _elem.setAttribute("db", "PDB")
                    _elem.setAttribute("dbkey", pdbe_id)
                    xref_list.appendChild(_elem)
                elem.appendChild(xref_list)

            # Write taxonomic distribution
            tax_dist = doc.createElement("taxonomy_distribution")
            for name, num_proteins in sorted(entry2ancestors[entry_acc]):
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


def export_matches(databases_file: str, entries_file: str, isoforms_file: str,
                   proteins_file: str, matches_file: str, outdir: str):
    shutil.copy(os.path.join(os.path.dirname(__file__), _MATCHES_DTD),
                outdir)

    logger.info("loading isoforms")
    protein2isoforms = {}
    with SimpleStore(isoforms_file) as store:
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

    logger.info("loading entries")
    entries = loadobj(entries_file)
    signatures = {}
    for entry in entries.values():
        if entry.database == "interpro":
            continue
        elif entry.integrated_in:
            interpro_entry = {
                "id": entry.integrated_in,
                "name": entries[entry.integrated_in].name,
                "type": entries[entry.integrated_in].type,
                "parent": entries[entry.integrated_in].relations[0]
            }
        else:
            interpro_entry = None

        signatures[entry.accession] = {
            "accession": entry.accession,
            "name": entry.name or entry.accession,
            "database": entry.source_database,
            "evidence": entry.evidence,
            "interpro": interpro_entry
        }

    logger.info("starting")
    file = os.path.join(outdir, _MATCHES_XML)
    with gzip.open(file, "wt", encoding="utf-8") as fh:
        fh.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        fh.write('<!DOCTYPE interpromatch SYSTEM "match_complete.dtd">\n')
        fh.write('<interpromatch>\n')

        doc = getDOMImplementation().createDocument(None, None, None)

        elem = doc.createElement("release")
        for db in loadobj(databases_file):
            db_altname = db[1]
            db_type = db[4]
            db_num_entries = db[5]
            db_version = db[6]
            db_date = db[7]

            if db_type == "entry":
                dbinfo = doc.createElement("dbinfo")
                dbinfo.setAttribute("dbname", db_altname)
                if db_version:
                    dbinfo.setAttribute("version", db_version)
                if db_num_entries:
                    dbinfo.setAttribute("entry_count", str(db_num_entries))
                if db_date:
                    dbinfo.setAttribute("file_date",
                                        db_date.strftime("%d-%b-%y").upper())
                elem.appendChild(dbinfo)

        elem.writexml(fh, addindent="  ", newl="\n")

        with Store(proteins_file) as proteins, Store(matches_file) as matches:
            for i, (protein_acc, protein) in enumerate(proteins.items()):
                elem = doc.createElement("protein")
                elem.setAttribute("id", protein_acc)
                elem.setAttribute("name", protein["identifier"])
                elem.setAttribute("length", str(protein["length"]))
                elem.setAttribute("crc64", protein["crc64"])

                protein_entries = matches.get(protein_acc, {})
                for signature_acc in sorted(protein_entries):
                    try:
                        signature = signatures[signature_acc]
                    except KeyError:
                        # InterPro entry
                        continue

                    locations = protein_entries[signature_acc]
                    elem.appendChild(_create_match(doc, signature, locations))

                elem.writexml(fh, addindent="  ", newl="\n")

                isoforms = protein2isoforms.get(protein_acc, [])
                for acc, length, crc64, _matches in isoforms:
                    elem = doc.createElement("protein")
                    elem.setAttribute("id", acc)
                    elem.setAttribute("name", acc)
                    elem.setAttribute("length", str(length))
                    elem.setAttribute("crc64", crc64)

                    for signature_acc in sorted(_matches):
                        try:
                            signature = signatures[signature_acc]
                        except KeyError:
                            # InterPro entry
                            continue

                        locations = _matches[signature_acc]
                        elem.appendChild(_create_match(doc, signature,
                                                       locations))

                    elem.writexml(fh, addindent="  ", newl="\n")

                if (i + 1) % 1e7 == 0:
                    logger.info(f"{i + 1:>15,}")

            logger.info(f"{i + 1:>15,}")

        fh.write('</interpromatch>\n')


def _create_match(doc, signature: dict, locations: Sequence[dict]):
    match = doc.createElement("match")
    match.setAttribute("id", signature["accession"])
    match.setAttribute("name", signature["name"])
    match.setAttribute("dbname", signature["database"])
    match.setAttribute("status", 'T')
    """
    The model is stored in locations, so we get the model
    from the first location for the match's 'model' attribute
    """
    match.setAttribute("model", locations[0]["model"])
    match.setAttribute("evd", signature["evidence"])

    if signature["interpro"]:
        ipr = doc.createElement("ipr")
        for attname, value in signature["interpro"]:
            if value:
                ipr.setAttribute(attname, value)

        match.appendChild(ipr)

    for loc in locations:
        match.appendChild(create_lcn(doc, loc))

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
        for db in sorted(loadobj(databases_file)):
            db_altname = db[1]
            db_type = db[4]
            db_version = db[6]
            if db_type == "feature":
                dbinfo = doc.createElement("dbinfo")
                dbinfo.setAttribute("dbname", db_altname)
                if db_version:
                    dbinfo.setAttribute("version", db_version)

                elem.appendChild(dbinfo)

        elem.writexml(fh, addindent="  ", newl="\n")

        with Store(proteins_file) as proteins, Store(features_file) as fts:
            for i, (protein_acc, protein_features) in enumerate(fts.items()):
                protein = proteins[protein_acc]
                elem = doc.createElement("protein")
                elem.setAttribute("id", protein_acc)
                elem.setAttribute("name", protein["identifier"])
                elem.setAttribute("length", str(protein["length"]))
                elem.setAttribute("crc64", protein["crc64"])

                for feature_acc in sorted(protein_features):
                    feature = protein_features[feature_acc]

                    match = doc.createElement("match")
                    match.setAttribute("id", feature_acc)
                    match.setAttribute("name", feature["name"])
                    match.setAttribute("dbname", feature["database"])
                    match.setAttribute("status", 'T')
                    match.setAttribute("model", feature_acc)
                    match.setAttribute("evd", feature["evidence"])

                    for loc in sorted(feature["locations"]):
                        # there is only one fragment per location
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

    logger.info("complete")


def export_structure_matches():
    pass
