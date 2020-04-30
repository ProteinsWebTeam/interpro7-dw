# -*- coding: utf-8 -*-

import gzip
import os
from typing import Sequence
from xml.dom.minidom import getDOMImplementation

import cx_Oracle
import MySQLdb

from interpro7dw import logger
from interpro7dw.ebi.interpro import production as ippro, utils
from interpro7dw.utils import DumpFile, Store, loadobj, url2dict

_DC_STATUSES = {v: k for k, v in utils.DC_STATUSES.items()}


def _create_lcn_elem(doc, location):
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
    lcn.setAttribute("stop", str(end))
    lcn.setAttribute("fragments", ','.join(fragments_obj))
    lcn.setAttribute("score", str(location["score"]))
    return lcn


def _create_match_elem(doc, signature: dict, locations: Sequence[dict]):
    match = doc.createElement("match")
    match.setAttribute("id", signature["accession"])
    match.setAttribute("name", signature["name"])
    match.setAttribute("dbname", signature["database"])
    match.setAttribute("status", 'T')
    """
    The model is stored in locations,  so we get the model 
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
        match.appendChild(_create_lcn_elem(doc, loc))

    return match


def export_interpro(url: str, p_entries: str, p_entry2xrefs: str, outdir: str):
    logger.info("loading protein counts")
    entries = loadobj(p_entries)
    num_proteins = {}
    interpro2pdb = {}
    with DumpFile(p_entry2xrefs) as entry2xrefs:
        for entry_acc, xrefs in entry2xrefs:
            entry = entries[entry_acc]
            if entry.database != "interpro" or entry.is_deleted:
                continue

            num_proteins[entry_acc] = str(len(xrefs["proteins"]))
            interpro2pdb[entry_acc] = xrefs["structures"]

    logger.info("writing XML")
    with gzip.open(os.path.join(outdir, "interpro.xml.gz"), "wt") as fh:
        fh.write('<?xml version="1.0" encoding="ISO-8859-1"?>\n')
        fh.write('<!DOCTYPE interprodb SYSTEM "interpro.dtd">\n')
        fh.write("<interprodb>\n")

        doc = getDOMImplementation().createDocument(None, None, None)

        elem = doc.createElement("release")
        databases = {}
        con = MySQLdb.connect(**url2dict(url))
        cur = con.cursor()
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
                dbinfo.setAttribute("dbname", name_alt)
                dbinfo.setAttribute("entry_count", str(entry_count))
                dbinfo.setAttribute("file_date",
                                    date.strftime("%d-%b-%Y").upper())
                elem.appendChild(dbinfo)
        cur.close()
        con.close()
        elem.writexml(fh, addindent="  ", newl="\n")

        for entry_acc in sorted(entries):
            entry = entries[entry_acc]
            if entry.database != "interpro" or entry.is_deleted:
                continue

            elem = doc.createElement("interpro")
            elem.setAttribute("id", entry.accession)
            elem.setAttribute("protein_count", num_proteins[entry_acc])
            elem.setAttribute("short_name", entry.short_name)
            elem.setAttribute("type", entry.type)

            abstract = doc.createElement("abstract")
            node = doc.createTextNode("\n".join(entry.description))
            abstract.appendChild(node)
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

                    if pub["raw_pages"] or pub["volume"] or pub["issue"]:
                        _elem = doc.createElement("location")
                        if pub["raw_pages"]:
                            _elem.setAttribute("pages", pub["raw_pages"])

                        if pub["volume"]:
                            _elem.setAttribute("volume", pub["volume"])

                        if pub["issue"]:
                            _elem.setAttribute("issue", pub["issue"])

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

            if entry_acc in interpro2pdb:
                xref_list = doc.createElement("structure_db_links")
                for pdb_id in sorted(interpro2pdb[entry_acc]):
                    _elem = doc.createElement("db_xref")
                    _elem.setAttribute("db", "PDB")
                    _elem.setAttribute("dbkey", pdb_id)
                    xref_list.appendChild(_elem)
                elem.appendChild(xref_list)

            # TODO: <taxonomy_distribution>
            #  <taxon_data name="Metazoa" proteins_count="32"/>

            elem.writexml(fh, addindent="  ", newl="\n")

        fh.write("</interprodb>\n")


def export_matches(pro_url: str, stg_url: str, p_proteins: str,
                   p_uniprot2matches: str, outdir: str):
    logger.info("loading isoforms")
    u2variants = {}
    for accession, variant in ippro.get_isoforms(pro_url).items():
        protein_acc = variant["protein"]
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

    logger.info("writing match_complete.xml.gz")
    with gzip.open(os.path.join(outdir, "match_complete.xml.gz"), "wt") as fh:
        fh.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        fh.write('<!DOCTYPE interpromatch SYSTEM "match_complete.dtd">\n')
        fh.write('<interpromatch>\n')

        doc = getDOMImplementation().createDocument(None, None, None)

        elem = doc.createElement("release")
        con = MySQLdb.connect(**url2dict(stg_url))
        cur = con.cursor()
        cur.execute(
            """
            SELECT name, name_alt, type, num_entries, version, release_date
            FROM webfront_database
            ORDER BY name_long
            """
        )

        for name, name_alt, db_type, entry_count, version, date in cur:
            if db_type == "entry":
                dbinfo = doc.createElement("dbinfo")
                dbinfo.setAttribute("dbname", name_alt)
                dbinfo.setAttribute("entry_count", str(entry_count))
                dbinfo.setAttribute("file_date",
                                    date.strftime("%d-%b-%Y").upper())
                elem.appendChild(dbinfo)
        cur.close()
        con.close()
        elem.writexml(fh, addindent="  ", newl="\n")

        proteins = Store(p_proteins)
        u2matches = Store(p_uniprot2matches)

        i = 0
        for uniprot_acc, protein in proteins.items():
            protein_entries = u2matches.get(uniprot_acc, {})

            if protein_entries:
                elem = doc.createElement("protein")
                elem.setAttribute("id", uniprot_acc)
                elem.setAttribute("name", protein["identifier"])
                elem.setAttribute("length", protein["length"])
                elem.setAttribute("crc64", protein["crc64"])

                for signature_acc in sorted(protein_entries):
                    try:
                        signature = signatures[signature_acc]
                    except KeyError:
                        # InterPro entry
                        continue

                    locations = protein_entries[signature_acc]
                    elem.appendChild(_create_match_elem(doc, signature,
                                                        locations))

                elem.writexml(fh, addindent="  ", newl="\n")

            protein_variants = u2variants.get(uniprot_acc, [])
            for variant, length, crc64, matches in protein_variants:
                elem = doc.createElement("protein")
                elem.setAttribute("id", variant)
                elem.setAttribute("name", variant)
                elem.setAttribute("length", length)
                elem.setAttribute("crc64", crc64)

                for signature_acc in sorted(matches):
                    try:
                        signature = signatures[signature_acc]
                    except KeyError:
                        # InterPro entry
                        continue

                    locations = matches[signature_acc]
                    elem.appendChild(_create_match_elem(doc, signature,
                                                        locations))

                elem.writexml(fh, addindent="  ", newl="\n")

            i += 1
            if not i % 10000000:
                logger.info(f"{i:>13}")

        logger.info(f"{i:>13}")
        proteins.close()
        u2matches.close()
        fh.write('</interpromatch>\n')

    logger.info("complete")
