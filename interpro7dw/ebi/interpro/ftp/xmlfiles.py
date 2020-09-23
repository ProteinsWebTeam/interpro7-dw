# -*- coding: utf-8 -*-

import json
import gzip
import math
import multiprocessing as mp
import os
import re
import shutil
from tempfile import mkstemp
from typing import Optional, Sequence
from xml.dom.minidom import getDOMImplementation, parseString
from xml.parsers.expat import ExpatError

import cx_Oracle
import MySQLdb
import MySQLdb.cursors

from interpro7dw import logger
from interpro7dw.ebi import pdbe
from interpro7dw.ebi.interpro import production as ippro, utils
from interpro7dw.utils import DumpFile, KVdb, Store, loadobj, url2dict

_DC_STATUSES = {v: k for k, v in utils.DC_STATUSES.items()}
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


def export_interpro(url: str, p_entries: str, p_entry2xrefs: str,
                    p_interpro2taxonomy: str, outdir: str,
                    tmpdir: Optional[str] = None):
    shutil.copy(os.path.join(os.path.dirname(__file__), "interpro.dtd"),
                outdir)

    logger.info("loading entries")
    entries = loadobj(p_entries)
    interpro_entries = []
    deleted_entries = []
    for e in entries.values():
        if e.database != "interpro":
            continue
        elif e.is_deleted:
            deleted_entries.append(e.accession)
        else:
            interpro_entries.append(e.accession)

    logger.info("creating entry-taxon database")
    fd, taxdb = mkstemp(dir=tmpdir)
    os.close(fd)
    os.remove(taxdb)
    with DumpFile(p_interpro2taxonomy) as interpro2taxonomy:
        with KVdb(taxdb, writeback=True) as kvdb:
            i = 0
            for entry_acc, taxon_id, counts in interpro2taxonomy:
                kvdb[f"{entry_acc}-{taxon_id}"] = str(counts)

                i += 1
                if not i % 1000000:
                    kvdb.sync()

    logger.info("loading protein counts")
    con = MySQLdb.connect(**url2dict(url))
    cur = MySQLdb.cursors.SSCursor(con)
    cur.execute(
        """
        SELECT accession, counts
        FROM webfront_entry
        """
    )
    num_proteins = {}
    for entry_acc, counts in cur:
        num_proteins[entry_acc] = str(json.loads(counts)["proteins"])

    output = os.path.join(outdir, "interpro.xml.gz")
    with gzip.open(output, "wt", encoding="utf-8") as fh:
        fh.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        fh.write('<!DOCTYPE interprodb SYSTEM "interpro.dtd">\n')
        fh.write("<interprodb>\n")

        doc = getDOMImplementation().createDocument(None, None, None)

        # writing <release> section (do not log progress, < 1 sec)
        elem = doc.createElement("release")
        databases = {}
        cur.execute(
            """
            SELECT name, name_alt, type, num_entries, version, release_date
            FROM webfront_database
            ORDER BY name_long
            """
        )

        for name, name_alt, db_type, entry_count, version, date in cur:
            databases[name] = name_alt
            if db_type == "entry":
                dbinfo = doc.createElement("dbinfo")
                dbinfo.setAttribute("version", version)
                dbinfo.setAttribute("dbname", name_alt)
                dbinfo.setAttribute("entry_count", str(entry_count))
                dbinfo.setAttribute("file_date",
                                    date.strftime("%d-%b-%y").upper())
                elem.appendChild(dbinfo)

        elem.writexml(fh, addindent="  ", newl="\n")

        logger.info("loading taxonomic data")
        key_species = {
            3702,    # Arabidopsis thaliana
            6239,    # Caenorhabditis elegans
            7955,    # Danio rerio
            7227,    # Drosophila melanogaster
            9606,    # Homo sapiens
            10090,   # Mus musculus
            367110,  # Neurospora crassa
            10116,   # Rattus norvegicus
            559292,  # Saccharomyces cerevisiae
            284812,  # Schizosaccharomyces pombe
            4577,    # Zea mays
        }
        superkingdoms = {
            "Archaea": None,
            "Bacteria": None,
            "Eukaryota": None,
            "Viruses": None
        }
        cur.execute(
            """
            SELECT accession, scientific_name, full_name, lineage
            FROM webfront_taxonomy
            """
        )
        taxa = {}
        for tax_id, sci_name, full_name, lineage in cur:
            """
            lineage stored as a string with heading/leading whitespaces,
            and a whitespace between taxa 
            """
            taxa[tax_id] = (full_name, lineage.strip().split())

            if sci_name in superkingdoms:
                superkingdoms[sci_name] = tax_id

        cur.close()
        con.close()

        # Raise if a superkingdom is not in the table
        for sci_name, tax_id in superkingdoms.items():
            if tax_id is None:
                raise ValueError(f"{sci_name}: missing taxon ID")

        superkingdoms = {tax_id for tax_id in superkingdoms.values()}

        logger.info("writing entries")
        with DumpFile(p_entry2xrefs) as entry2xrefs, KVdb(taxdb) as kvdb:
            for entry_acc, xrefs in entry2xrefs:
                entry = entries[entry_acc]
                if entry.database != "interpro" or entry.is_deleted:
                    continue

                elem = doc.createElement("interpro")
                elem.setAttribute("id", entry.accession)
                elem.setAttribute("protein_count", num_proteins[entry_acc])
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
                            num_proteins[signature_acc],
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

                if entry.cross_references:
                    xref_list = doc.createElement("external_doc_list")
                    for ref_db in sorted(entry.cross_references):
                        for ref_id in sorted(entry.cross_references[ref_db]):
                            _elem = doc.createElement("db_xref")
                            _elem.setAttribute("db", databases[ref_db])
                            _elem.setAttribute("dbkey", ref_id)
                            xref_list.appendChild(_elem)
                    elem.appendChild(xref_list)

                if xrefs["structures"]:
                    xref_list = doc.createElement("structure_db_links")
                    for pdb_id in sorted(xrefs["structures"]):
                        _elem = doc.createElement("db_xref")
                        _elem.setAttribute("db", "PDB")
                        _elem.setAttribute("dbkey", pdb_id)
                        xref_list.appendChild(_elem)
                    elem.appendChild(xref_list)

                # Find key species and taxonomic distribution
                entry_key_species = []
                entry_superkingdoms = {}
                for tax_id in xrefs["taxa"]:
                    full_name, lineage = taxa[tax_id]

                    if tax_id in key_species:
                        entry_key_species.append((full_name, tax_id))

                    # Find the superkingdom contain this taxon
                    for superkingdom_id in superkingdoms:
                        if superkingdom_id in lineage:
                            break
                    else:
                        continue

                    try:
                        other_lineage = entry_superkingdoms[superkingdom_id]
                    except KeyError:
                        entry_superkingdoms[superkingdom_id] = lineage
                    else:
                        # Compare lineages and find lowest common ancestor
                        i = 0
                        while i < len(lineage) and i < len(other_lineage):
                            if lineage[i] != other_lineage[i]:
                                break
                            i += 1

                        # Path to the lowest common ancestor
                        entry_superkingdoms[superkingdom_id] = lineage[:i]

                # Get lowest common ancestor for each represented superkingdom
                lowest_common_ancestors = []
                for lineage in entry_superkingdoms.values():
                    # Lowest common ancestor
                    tax_id = lineage[-1]
                    full_name, _ = taxa[tax_id]
                    lowest_common_ancestors.append((full_name, tax_id))

                # Write taxonomic distribution
                tax_dist = doc.createElement("taxonomy_distribution")
                for full_name, tax_id in sorted(lowest_common_ancestors):
                    _elem = doc.createElement("taxon_data")
                    _elem.setAttribute("name", full_name)
                    key = f"{entry_acc}-{tax_id}"
                    _elem.setAttribute("proteins_count", kvdb[key])
                    tax_dist.appendChild(_elem)
                elem.appendChild(tax_dist)

                if entry_key_species:
                    # Write key species
                    key_spec = doc.createElement("key_species")
                    for full_name, tax_id in sorted(entry_key_species):
                        _elem = doc.createElement("taxon_data")
                        _elem.setAttribute("name", full_name)
                        key = f"{entry_acc}-{tax_id}"
                        _elem.setAttribute("proteins_count", kvdb[key])
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

    logger.info(f"temporary file: {os.path.getsize(taxdb)/1024/1024:,.0f} MB")
    os.remove(taxdb)
    logger.info("complete")


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
    We do not have to orginal start/end match positions, 
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


def _write_match_tmp(signatures: dict, u2variants: dict, p_proteins: str,
                     p_uniprot2matches: str, start: str, stop: Optional[str],
                     output: str):
    proteins = Store(p_proteins)
    u2matches = Store(p_uniprot2matches)
    with open(output, "wt", encoding="utf-8") as fh:
        doc = getDOMImplementation().createDocument(None, None, None)

        for uniprot_acc, protein in proteins.range(start, stop):
            elem = doc.createElement("protein")
            elem.setAttribute("id", uniprot_acc)
            elem.setAttribute("name", protein["identifier"])
            elem.setAttribute("length", str(protein["length"]))
            elem.setAttribute("crc64", protein["crc64"])

            try:
                protein_entries = u2matches[uniprot_acc]
            except KeyError:
                pass
            else:
                for signature_acc in sorted(protein_entries):
                    try:
                        signature = signatures[signature_acc]
                    except KeyError:
                        # InterPro entry
                        continue

                    elem.appendChild(
                        _create_match(doc, signature,
                                      protein_entries[signature_acc])
                    )
            finally:
                elem.writexml(fh, addindent="  ", newl="\n")

            protein_variants = u2variants.get(uniprot_acc, [])
            for variant, length, crc64, matches in protein_variants:
                elem = doc.createElement("protein")
                elem.setAttribute("id", variant)
                elem.setAttribute("name", variant)
                elem.setAttribute("length", str(length))
                elem.setAttribute("crc64", crc64)

                for signature_acc in sorted(matches):
                    try:
                        signature = signatures[signature_acc]
                    except KeyError:
                        # InterPro entry
                        continue

                    elem.appendChild(
                        _create_match(doc, signature,
                                      matches[signature_acc])
                    )

                elem.writexml(fh, addindent="  ", newl="\n")


def export_matches(pro_url: str, stg_url: str, p_proteins: str,
                   p_uniprot2matches: str, outdir: str, processes: int = 8):
    shutil.copy(os.path.join(os.path.dirname(__file__), "match_complete.dtd"),
                outdir)

    logger.info("loading isoforms")
    u2variants = {}
    for accession, variant in ippro.get_isoforms(pro_url).items():
        protein_acc = variant["protein_acc"]
        try:
            variants = u2variants[protein_acc]
        except KeyError:
            variants = u2variants[protein_acc] = []
        finally:
            variants.append((
                accession,
                variant["length"],
                variant["crc64"],
                variant["matches"]
            ))

    logger.info("loading signatures")
    con = cx_Oracle.connect(pro_url)
    cur = con.cursor()
    signatures = ippro.get_signatures(cur)
    cur.close()
    con.close()

    logger.info("spawning processes")
    processes = max(1, processes - 1)
    ctx = mp.get_context(method="spawn")
    workers = []
    with Store(p_proteins) as proteins:
        proteins_per_file = math.ceil(len(proteins) / processes)
        start_acc = None
        for i, uniprot_acc in enumerate(proteins):
            if not i % proteins_per_file:
                if start_acc:
                    filename = f"match_{len(workers)+1}.xml"
                    filepath = os.path.join(outdir, filename)
                    p = ctx.Process(target=_write_match_tmp,
                                    args=(signatures, u2variants, p_proteins,
                                          p_uniprot2matches, start_acc,
                                          uniprot_acc, filepath))
                    p.start()
                    workers.append((p, filepath))

                start_acc = uniprot_acc

        filename = f"match_{len(workers) + 1}.xml"
        filepath = os.path.join(outdir, filename)
        p = ctx.Process(target=_write_match_tmp,
                        args=(signatures, u2variants, p_proteins,
                              p_uniprot2matches, start_acc, None, filepath))
        p.start()
        workers.append((p, filepath))

    logger.info("concatenating XML files")
    con = MySQLdb.connect(**url2dict(stg_url))
    cur = con.cursor()
    cur.execute(
        """
        SELECT name, name_alt, type, num_entries, version, release_date
        FROM webfront_database
        ORDER BY name_long
        """
    )

    doc = getDOMImplementation().createDocument(None, None, None)
    elem = doc.createElement("release")
    for name, name_alt, db_type, entry_count, version, date in cur:
        if db_type == "entry":
            dbinfo = doc.createElement("dbinfo")
            dbinfo.setAttribute("dbname", name_alt)
            if version:
                dbinfo.setAttribute("version", version)
            if entry_count:
                dbinfo.setAttribute("entry_count", str(entry_count))
            if date:
                dbinfo.setAttribute("file_date",
                                    date.strftime("%d-%b-%y").upper())
            elem.appendChild(dbinfo)
    cur.close()
    con.close()

    output = os.path.join(outdir, "match_complete.xml.gz")
    with gzip.open(output, "wt", encoding="utf-8") as fh:
        fh.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        fh.write('<!DOCTYPE interpromatch SYSTEM "match_complete.dtd">\n')
        fh.write('<interpromatch>\n')
        elem.writexml(fh, addindent="  ", newl="\n")

        for i, (p, filepath) in enumerate(workers):
            p.join()
            with open(filepath, "rt", encoding="utf-8") as tfh:
                for line in tfh:
                    fh.write(line)

            os.remove(filepath)
            logger.info(f"\t{i+1} / {len(workers)}")

        fh.write('</interpromatch>\n')

    logger.info("complete")


def _write_feature_tmp(features: dict, p_proteins: str,
                       p_uniprot2features: str, start: str,
                       stop: Optional[str], output: str):
    proteins = Store(p_proteins)
    u2features = Store(p_uniprot2features)

    with open(output, "wt", encoding="utf-8") as fh:
        doc = getDOMImplementation().createDocument(None, None, None)

        # for uniprot_acc, protein in proteins.range(start, stop):
        for uniprot_acc, protein_features in u2features.range(start, stop):
            protein = proteins[uniprot_acc]
            elem = doc.createElement("protein")
            elem.setAttribute("id", uniprot_acc)
            elem.setAttribute("name", protein["identifier"])
            elem.setAttribute("length", str(protein["length"]))
            elem.setAttribute("crc64", protein["crc64"])

            for feature_acc in sorted(protein_features):
                feature = features[feature_acc]
                feature_match = protein_features[feature_acc]

                match = doc.createElement("match")
                match.setAttribute("id", feature_acc)
                match.setAttribute("name", feature["name"])
                match.setAttribute("dbname", feature["database"])
                match.setAttribute("status", 'T')
                match.setAttribute("model", feature_acc)
                match.setAttribute("evd", feature["evidence"])

                for loc in feature_match["locations"]:
                    # there is only one fragment per location
                    frag = loc["fragments"][0]

                    lcn = doc.createElement("lcn")
                    lcn.setAttribute("start", str(frag["start"]))
                    lcn.setAttribute("end", str(frag["end"]))
                    lcn.setAttribute("fragments",
                                     f"{frag['start']}-"
                                     f"{frag['end']}-"
                                     f"{_DC_STATUSES['CONTINUOUS']}")

                    if frag["seq_feature"]:
                        lcn.setAttribute("sequence-feature",
                                         frag["seq_feature"])

                    match.appendChild(lcn)

                elem.appendChild(match)

            elem.writexml(fh, addindent="  ", newl="\n")


def export_features_matches(url: str, p_proteins: str, p_uniprot2features: str,
                            outdir: str, processes: int = 8):
    shutil.copy(os.path.join(os.path.dirname(__file__), "extra.dtd"),
                outdir)

    logger.info("loading features")
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    features = ippro.get_features(cur)
    cur.close()
    con.close()

    logger.info("spawning processes")
    processes = max(1, processes - 1)
    ctx = mp.get_context(method="spawn")
    workers = []
    with Store(p_uniprot2features) as proteins:
        proteins_per_file = math.ceil(len(proteins) / processes)
        start_acc = None
        for i, uniprot_acc in enumerate(proteins):
            if not i % proteins_per_file:
                if start_acc:
                    filename = f"extra_{len(workers) + 1}.xml"
                    filepath = os.path.join(outdir, filename)
                    p = ctx.Process(target=_write_feature_tmp,
                                    args=(features, p_proteins,
                                          p_uniprot2features, start_acc,
                                          uniprot_acc, filepath))
                    p.start()
                    workers.append((p, filepath))

                start_acc = uniprot_acc

        filename = f"extra_{len(workers) + 1}.xml"
        filepath = os.path.join(outdir, filename)
        p = ctx.Process(target=_write_feature_tmp,
                        args=(features, p_proteins, p_uniprot2features,
                              start_acc, None, filepath))
        p.start()
        workers.append((p, filepath))

    logger.info("concatenating XML files")
    output = os.path.join(outdir, "extra.xml.gz")
    with gzip.open(output, "wt", encoding="utf-8") as fh:
        fh.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        fh.write('<!DOCTYPE interproextra SYSTEM "extra.dtd">\n')
        fh.write('<interproextra>\n')

        doc = getDOMImplementation().createDocument(None, None, None)
        elem = doc.createElement("release")
        databases = {(f["database"], f["version"]) for f in features.values()}
        for name, version in sorted(databases):
            dbinfo = doc.createElement("dbinfo")
            dbinfo.setAttribute("dbname", name)

            if version:
                dbinfo.setAttribute("version", version)

            elem.appendChild(dbinfo)

        elem.writexml(fh, addindent="  ", newl="\n")

        for i, (p, filepath) in enumerate(workers):
            p.join()
            with open(filepath, "rt", encoding="utf-8") as tfh:
                for line in tfh:
                    fh.write(line)

            os.remove(filepath)
            logger.info(f"\t{i+1} / {len(workers)}")

        fh.write('</interproextra>\n')

    logger.info("complete")


def export_structure_matches(url: str, p_proteins: str, p_structures: str,
                             outdir:str):
    shutil.copy(os.path.join(os.path.dirname(__file__), "feature.dtd"),
                outdir)

    logger.info("loading structures")
    uniprot2pdbe = {}
    for pdb_id, entry in loadobj(p_structures).items():
        for uniprot_acc, chains in entry["proteins"].items():
            try:
                uniprot2pdbe[uniprot_acc][pdb_id] = chains
            except KeyError:
                uniprot2pdbe[uniprot_acc] = {pdb_id: chains}

    logger.info("loading CATH/SCOP domains")
    uni2prot2cath = pdbe.get_cath_domains(url)
    uni2prot2scop = pdbe.get_scop_domains(url)

    logger.info("writing file")
    output = os.path.join(outdir, "feature.xml.gz")
    with gzip.open(output, "wt", encoding="utf-8") as fh:
        fh.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        fh.write('<!DOCTYPE interprofeature SYSTEM "feature.dtd">\n')
        fh.write('<interprofeature>\n')

        with Store(p_proteins) as proteins:
            doc = getDOMImplementation().createDocument(None, None, None)

            for uniprot_acc, protein in proteins.items():
                pdb_entries = uniprot2pdbe.get(uniprot_acc, {})
                cath_entries = uni2prot2cath.get(uniprot_acc, {})
                scop_entries = uni2prot2scop.get(uniprot_acc, {})

                if pdb_entries or cath_entries or scop_entries:
                    elem = doc.createElement("protein")
                    elem.setAttribute("id", uniprot_acc)
                    elem.setAttribute("name", protein["identifier"])
                    elem.setAttribute("length", str(protein["length"]))
                    elem.setAttribute("crc64", protein["crc64"])

                    for pdb_id in sorted(pdb_entries):
                        chains = pdb_entries[pdb_id]
                        for chain_id in sorted(chains):
                            domain = doc.createElement("domain")
                            domain.setAttribute("id", f"{pdb_id}{chain_id}")
                            domain.setAttribute("dbname", "PDB")

                            for loc in chains[chain_id]:
                                start = loc["protein_start"]
                                end = loc["protein_end"]

                                coord = doc.createElement("coord")
                                coord.setAttribute("pdb", pdb_id)
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
