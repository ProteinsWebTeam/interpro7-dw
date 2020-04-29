# -*- coding: utf-8 -*-

import gzip
import multiprocessing as mp
import os
import tarfile
from datetime import datetime
from typing import List, Optional, Sequence, Tuple
from xml.dom.minidom import getDOMImplementation

import cx_Oracle
import MySQLdb

from interpro7dw import logger
from interpro7dw.ebi import pdbe
from interpro7dw.ebi.interpro import production as ippro, utils
from interpro7dw.utils import DirectoryTree, DumpFile, KVdb, Store
from interpro7dw.utils import loadobj, url2dict

_DC_STATUSES = {v: k for k, v in utils.DC_STATUSES.items()}


def _export_pdb2interpro2go2uniprot(cur: cx_Oracle.Cursor, output: str):
    # PDBe sequences from UniParc
    cur.execute(
        """
        SELECT UPI, AC
        FROM UNIPARC.XREF
        WHERE DBID = 21
        AND DELETED = 'N'
        """
    )
    sequences = {}
    for upi, pdb_acc in cur:
        if upi in sequences:
            sequences[upi]["structures"].add(pdb_acc)
        else:
            sequences[upi] = {
                "structures": {pdb_acc},
                "entries": set()
            }

    # Integrated signatures whose entry is checked and has GO terms
    cur.execute(
        """
        SELECT METHOD_AC, ENTRY_AC
        FROM INTERPRO.ENTRY2METHOD
        WHERE ENTRY_AC IN (
          SELECT ENTRY_AC
          FROM INTERPRO.ENTRY
          WHERE CHECKED = 'Y'
          INTERSECT
          SELECT DISTINCT ENTRY_AC
          FROM INTERPRO.INTERPRO2GO
        )
        """
    )
    signatures = dict(cur.fetchall())

    # GO terms in InterPro
    cur.execute(
        """
        SELECT ENTRY_AC, GO_ID
        FROM INTERPRO.INTERPRO2GO
        WHERE ENTRY_AC IN (
          SELECT ENTRY_AC
          FROM INTERPRO.ENTRY
          WHERE CHECKED = 'Y'
        )
        """
    )
    entries = {}
    for entry_acc, go_id in cur:
        if entry_acc in entries:
            entries[entry_acc].add(go_id)
        else:
            entries[entry_acc] = {go_id}

    # PDBe matches
    cur.execute(
        """
        SELECT DISTINCT UPI, METHOD_AC
        FROM IPRSCAN.MV_IPRSCAN
        WHERE UPI IN (
            SELECT UPI
            FROM UNIPARC.XREF
            WHERE DBID = 21
            AND DELETED = 'N'
        )
        """
    )
    for upi, signature_acc in cur:
        try:
            entry_acc = signatures[signature_acc]
        except KeyError:
            pass
        else:
            sequences[upi]["entries"].add(entry_acc)

    # PDBe taxonomy
    structures = pdbe.get_chain_taxonomy(cur)

    # UniProt accessions
    cur.execute(
        """
        SELECT DISTINCT A.AC, B.AC
        FROM UNIPARC.XREF A
        LEFT OUTER JOIN UNIPARC.XREF B ON A.UPI = B.UPI
        WHERE A.DBID = 21
        AND A.DELETED = 'N'
        AND B.DBID IN (2, 3)
        AND B.DELETED = 'N'
        """
    )
    pdb2uniprot = {}
    for pdb_acc, protein_acc in cur:
        if not protein_acc:
            continue
        elif pdb_acc in pdb2uniprot:
            pdb2uniprot[pdb_acc].add(protein_acc)
        else:
            pdb2uniprot[pdb_acc] = {protein_acc}

    tmp_path = f"{output}.tmp"
    with open(tmp_path, "wt") as fh:
        fh.write("#PDBe ID\tchain\tTaxon ID\t"
                 "InterPro accession\tGO ID\tUniProt accession\n")

        for seq in sequences.values():
            for pdb_acc in seq["structures"]:
                try:
                    s = structures[pdb_acc]
                except KeyError:
                    # Structure does not exist in PDBe database
                    continue

                pdb_id = s["id"]
                chain = s["chain"]
                proteins = pdb2uniprot.get(pdb_acc, {''})

                for tax_id in s["taxa"]:
                    for entry_acc in seq["entries"]:
                        for go_id in entries[entry_acc]:
                            for protein_acc in proteins:
                                fh.write(f"{pdb_id}\t{chain}\t{tax_id}\t"
                                         f"{entry_acc}\t{go_id}\t"
                                         f"{protein_acc}\n")

    try:
        os.remove(output)
    except FileNotFoundError:
        pass
    finally:
        os.rename(tmp_path, output)
        os.chmod(output, 0o775)


def _export_interpro2go2uniprot(cur: cx_Oracle.Cursor, output: str):
    cur.execute(
        """
        SELECT DISTINCT IG.ENTRY_AC, IG.GO_ID, M.PROTEIN_AC
        FROM INTERPRO.MATCH M
        INNER JOIN INTERPRO.ENTRY2METHOD EM
          ON M.METHOD_AC = EM.METHOD_AC
        INNER JOIN INTERPRO.ENTRY E
          ON EM.ENTRY_AC = E.ENTRY_AC
        INNER JOIN INTERPRO.INTERPRO2GO IG
          ON E.ENTRY_AC = IG.ENTRY_AC
        WHERE E.CHECKED = 'Y'
        """
    )

    tmp_path = f"{output}.tmp"
    with open(tmp_path, "wt") as fh:
        fh.write("#InterPro accession\tGO ID\tUniProt accession\n")

        for row in cur:
            fh.write('\t'.join(row) + '\n')

    try:
        os.remove(output)
    except FileNotFoundError:
        pass
    finally:
        os.rename(tmp_path, output)
        os.chmod(output, 0o775)


def export_goa_mapping(pro_url: str, stg_url: str, outdir: str):
    os.makedirs(outdir, exist_ok=True)

    con = cx_Oracle.connect(pro_url)
    cur = con.cursor()

    logger.info("exporting PDB-InterPro-GO-UniProt mapping")
    filepath = os.path.join(outdir, "pdb2interpro2go.tsv")
    _export_pdb2interpro2go2uniprot(cur, filepath)

    logger.info("exporting InterPro-GO-UniProt mapping")
    filepath = os.path.join(outdir, "interpro2go2uniprot.tsv")
    _export_interpro2go2uniprot(cur, filepath)
    cur.close()
    con.close()

    logger.info("exporting release info")
    con = MySQLdb.connect(**url2dict(stg_url))
    cur = con.cursor()
    cur.execute(
        """
        SELECT version, release_date 
        FROM webfront_database 
        WHERE name='interpro'
        """
    )
    version, date = cur.fetchone()
    cur.close()
    con.close()

    filepath = os.path.join(outdir, "release.txt")
    with open(filepath, "wt") as fh:
        fh.write(f"InterPro version:    {version}\n")
        fh.write(f"Release date:        {date:%A, %d %B %Y}\n")
        fh.write(f"Generated on:        {datetime.now():%Y-%m-%d %H:%M}\n")

    os.chmod(filepath, 0o775)
    logger.info("complete")


def _write_node(node, fh, level):
    fh.write(f"{'-'*2*level}{node['accession']}::{node['name']}\n")

    for child in node["children"]:
        _write_node(child, fh, level+1)


def export_flat_files(p_entries: str, p_uniprot2matches: str, outdir: str):
    logger.info("loading entries")
    entries = []
    integrated = {}
    for e in loadobj(p_entries).values():
        if e.database == "interpro" and not e.is_deleted:
            entries.append(e)

            for signatures in e.integrates.values():
                for signature_acc in signatures:
                    integrated[signature_acc] = (e.accession, e.name)

    logger.info("writing entry.list")
    with open(os.path.join(outdir, "entry.list"), "wt") as fh:
        fh.write("ENTRY_AC\tENTRY_TYPE\tENTRY_NAME\n")

        for e in sorted(entries, key=lambda e: (e.type, e.accession)):
            fh.write(f"{e.accession}\t{e.type}\t{e.name}\n")

    logger.info("writing names.dat")
    with open(os.path.join(outdir, "names.dat"), "wt") as fh:
        for e in sorted(entries, key=lambda e: e.accession):
            fh.write(f"{e.accession}\t{e.name}\n")

    logger.info("writing short_names.dat")
    with open(os.path.join(outdir, "short_names.dat"), "wt") as fh:
        for e in sorted(entries, key=lambda e: e.accession):
            fh.write(f"{e.accession}\t{e.short_name}\n")

    logger.info("writing interpro2go")
    with open(os.path.join(outdir, "interpro2go"), "wt") as fh:
        fh.write(f"!date: {datetime.now():%Y/%m/%d %H:%M:%S}\n")
        fh.write("!Mapping of InterPro entries to GO\n")
        fh.write("!\n")

        for e in sorted(entries, key=lambda e: e.accession):
            for term in e.go_terms:
                fh.write(f"InterPro:{e.accession} {e.name} > "
                         f"GO:{term['name']} ; {term['identifier']}\n")

    logger.info("writing ParentChildTreeFile.txt")
    with open(os.path.join(outdir, "ParentChildTreeFile.txt"), "wt") as fh:
        for e in sorted(entries, key=lambda e: e.accession):
            root = e.hierarchy["accession"]
            if root == e.accession and e.hierarchy["children"]:
                _write_node(e.hierarchy, fh, level=0)

    logger.info("writing protein2ipr.dat.gz")
    filepath = os.path.join(outdir, "protein2ipr.dat.gz")
    with gzip.open(filepath, "wt") as fh, Store(p_uniprot2matches) as sh:
        i = 0
        for uniprot_acc, protein_entries in sh.items():
            matches = []
            for signature_acc in sorted(protein_entries):
                try:
                    interpro_acc, name = integrated[signature_acc]
                except KeyError:
                    # Not integrated signature or InterPro entry
                    continue

                locations = protein_entries[signature_acc]

                for loc in locations:
                    matches.append((
                        uniprot_acc,
                        interpro_acc,
                        name,
                        signature_acc,
                        # We do not consider fragmented locations
                        loc["fragments"][0]["start"],
                        max(f["end"] for f in loc["fragments"])
                    ))

            for m in sorted(matches):
                fh.write('\t'.join(map(str, m)) + '\n')

            i += 1
            if not i % 10000000:
                logger.info(f"{i:>12,}")

        logger.info(f"{i:>12,}")

    logger.info("complete")


def _post_matches(matches: Sequence[dict]) -> List[Tuple[str, str, List]]:
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


def _dump_proteins(proteins_file: str, matches_file: str, signatures: dict,
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

                        for start, stop, score, aln, frags in locations:
                            lcn = doc.createElement("lcn")
                            lcn.setAttribute("start", str(start))
                            lcn.setAttribute("stop", str(stop))

                            if frags:
                                lcn.setAttribute("fragments", frags)

                            if aln:
                                lcn.setAttribute("alignment", aln)

                            lcn.setAttribute("score", str(score))
                            match.appendChild(lcn)

                        protein.appendChild(match)

                    protein.writexml(fh, addindent="  ", newl="\n")

            outqueue.put(filepath)


def export_uniparc(url: str, outdir: str, dir: Optional[str]=None,
                   processes: int=4, proteins_per_file: int=1000000):
    dt = DirectoryTree(root=dir)
    proteins_file = dt.mktemp()
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

            if not i % 1000:
                keys.append(upi)

        kvdb.sync()

    logger.info("exporting UniParc matches")
    matches_file = dt.mktemp()
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
                    logger.debug(f"{i:>15,}")

        logger.debug(f"{i:>15,}")
        size = store.merge(fn=_post_matches, processes=processes)

    logger.info("loading signatures")
    signatures = ippro.get_signatures(cur)
    cur.close()
    con.close()

    logger.info("writing XML files")
    ctx = mp.get_context(method="spawn")
    inqueue = ctx.Queue()
    outqueue = ctx.Queue()
    workers = []
    for _ in range(max(1, processes - 1)):
        p = ctx.Process(target=_dump_proteins,
                        args=(proteins_file, matches_file, signatures,
                              inqueue, outqueue))
        p.start()
        workers.append(p)

    with KVdb(proteins_file) as kvdb, Store(matches_file) as store:
        num_files = 0

        i = 0
        from_upi = None
        for upi in store:
            try:
                length, crc64 = kvdb[upi]
            except KeyError:
                continue

            i += 1
            if not i % 10000000:
                logger.debug(f"{i:>13,}")

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
        logger.debug(f"{i:>13,}")

    for _ in workers:
        inqueue.put(None)

    logger.info("creating XML archive")
    output = os.path.join(outdir, "uniparc_match.tar.gz")
    with tarfile.open(output, "w:gz") as fh:
        for _ in range(num_files):
            filepath = outqueue.get()
            fh.add(filepath, arcname=os.path.basename(filepath))
            os.remove(filepath)

    for p in workers:
        p.join()

    size += os.path.getsize(proteins_file)
    size += os.path.getsize(matches_file)
    logger.info(f"temporary files: {size/1024/1024:.0f} MB")

    dt.remove()
    logger.info("complete")


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


def export_matches(pro_url: str, p_proteins: str, p_uniprot2matches: str,
                   outdir: str):
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
