# -*- coding: utf-8 -*-

import json
from multiprocessing import Process, Queue

import MySQLdb

from i7dw import io, logger, pdbe
from i7dw.interpro import condense_locations
from i7dw.interpro import DomainArchitecture, Table
from i7dw.interpro import oracle
from .entries import get_entries, iter_sets
from .structures import iter_structures
from .taxonomy import iter_taxa
from .utils import drop_index, parse_url


def _insert(url: str, queue: Queue):
    query = """
        INSERT INTO webfront_protein (accession, identifier, organism, name,
                                      description, sequence,
                                      length, proteome, gene, go_terms,
                                      evidence_code, source_database,
                                      residues, is_fragment, structure,
                                      tax_id, extra_features, ida_id, ida,
                                      counts)
        VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s,
                %s, %s, %s, %s, %s, %s)
    """
    con = MySQLdb.connect(**parse_url(url), charset="utf8")
    with Table(con, query) as table:
        for records in iter(queue.get, None):
            for record in records:
                table.insert(record)

    con.commit()


def insert_proteins(my_url: str, ora_ippro_url: str, ora_pdbe_url: str,
                    src_comments: str, src_features: str, src_matches: str,
                    src_misc: str, src_names: str, src_proteins: str,
                    src_proteomes: str, src_residues: str, src_sequences: str,
                    chunk_size: int = 10000, processes: int=1):
    logger.info("calculating domain architectures")
    proteins = io.Store(src_proteins)
    matches = io.Store(src_matches)
    dom_arch = DomainArchitecture(get_entries(my_url))
    dom_cnts = {}

    for protein_acc, protein_info in proteins:
        dom_arch.update(matches.get(protein_acc, {}))
        dom_id = dom_arch.identifier

        try:
            dom_cnts[dom_id] += 1
        except KeyError:
            dom_cnts[dom_id] = 1

    workers = []
    queue = Queue(maxsize=processes)
    for _ in range(max(1, processes-1)):
        p = Process(target=_insert, args=(my_url, queue))
        p.start()
        workers.append(p)

    logger.info("loading data")
    # Structural features (CATH and SCOP domains)
    cath_domains = pdbe.get_cath_domains(ora_pdbe_url)
    scop_domains = pdbe.get_scop_domains(ora_pdbe_url)

    # MySQL data
    entries = get_entries(my_url)
    taxa = {}
    for taxon in iter_taxa(my_url, lineage=False):
        taxa[taxon["id"]] = json.dumps({
            "taxId": taxon["id"],
            "scientificName": taxon["scientific_name"],
            "fullName": taxon["full_name"]
        })

    structures = {}
    for s in iter_structures(my_url):
        pdbe_id = s["accession"]
        for protein_acc in s["proteins"]:
            try:
                structures[protein_acc].add(pdbe_id)
            except KeyError:
                structures[protein_acc] = {pdbe_id}

    entry_set = {}
    for s in iter_sets(my_url):
        set_acc = s["accession"]
        for entry_acc in s["members"]:
            entry_set[entry_acc] = set_acc

    con = MySQLdb.connect(**parse_url(my_url), charset="utf8")
    cur = con.cursor()
    cur.execute(
        """
        SELECT protein_acc, COUNT(*)
        FROM webfront_varsplic
        GROUP BY protein_acc
        """
    )
    isoforms = dict(cur.fetchall())

    logger.info("preparing MySQL")
    cur.execute("TRUNCATE TABLE webfront_protein")
    for idx in ("ui_webfront_protein_identifier", "i_webfront_protein_length",
                "i_webfront_protein_ida", "i_webfront_protein_fragment"):
        drop_index(cur, "webfront_protein", idx)
    cur.close()
    con.close()

    logger.info("inserting proteins")
    comments = io.Store(src_comments)
    features = io.Store(src_features)
    misc = io.Store(src_misc)
    names = io.Store(src_names)
    proteomes = io.Store(src_proteomes)
    residues = io.Store(src_residues)
    sequences = io.Store(src_sequences)
    n_proteins = 0
    chunk = []
    for protein_acc, protein_info in proteins:
        tax_id = protein_info["taxon"]
        try:
            taxon_json = taxa[tax_id]
        except KeyError:
            logger.debug(f"{protein_acc}: unknown taxon '{tax_id}'")
            continue

        try:
            sequence = sequences[protein_acc]
        except KeyError:
            logger.debug(f"{protein_acc}: no sequence")
            continue

        try:
            evidence, gene = misc[protein_acc]
        except KeyError:
            logger.debug(f"{protein_acc}: no evidence/gene")
            continue

        if not evidence:
            logger.debug(f"{protein_acc}: no evidence")
            continue

        go_terms = {}  # InterPro2GO + InterPro matches -> UniProt-GO
        protein_entries = {}
        protein_matches = matches.get(protein_acc, {})
        protein_sets = set()
        for entry_acc in protein_matches:
            e = entries[entry_acc]

            for term in e["go_terms"]:
                go_terms[term["identifier"]] = term

            try:
                set_acc = entry_set[entry_acc]
            except KeyError:
                pass
            else:
                protein_sets.add(set_acc)

            database = e["database"]
            try:
                protein_entries[database] += 1
            except KeyError:
                protein_entries[database] = 1

        protein_entries["total"] = sum(protein_entries.values())
        protein_structures = {
            "cath": {},
            "scop": {}
        }
        for pdbe_id in structures.get(protein_acc, []):
            for dom_id, dom in cath_domains.get(pdbe_id, {}).items():
                protein_structures["cath"][dom_id] = dom

            for dom_id, dom in scop_domains.get(pdbe_id, {}).items():
                protein_structures["scop"][dom_id] = dom

        upid = proteomes.get(protein_acc)

        dom_arch.update(protein_matches)
        dom_ac = dom_arch.accession
        dom_id = dom_arch.identifier
        dom_cnt = dom_cnts[dom_id]

        # Enqueue record for protein table
        chunk.append((
            protein_acc,
            protein_info["identifier"],
            taxon_json,
            names.get(protein_acc),
            json.dumps(comments.get(protein_acc, [])),
            sequence,
            protein_info["length"],
            upid,
            gene,
            json.dumps(list(go_terms.values())),
            evidence,
            "reviewed" if protein_info["is_reviewed"] else "unreviewed",
            json.dumps(residues.get(protein_acc, {})),
            1 if protein_info["is_fragment"] else 0,
            json.dumps(protein_structures),
            tax_id,
            json.dumps(features.get(protein_acc, {})),
            dom_id,
            dom_ac,
            json.dumps({
                "entries": protein_entries,
                "structures": len(structures.get(protein_acc, [])),
                "sets": len(protein_sets),
                "proteomes": 1 if upid else 0,
                "taxa": 1,
                "idas": dom_cnt,
                "isoforms": isoforms.get(protein_acc, 0)
            })
        ))

        if len(chunk) == chunk_size:
            queue.put(chunk)
            chunk = []

        n_proteins += 1
        if not n_proteins % 10000000:
            logger.info(f"{n_proteins:>12,}")

    queue.put(chunk)
    chunk = []

    logger.info(f"{n_proteins:>12,}")
    for _ in workers:
        queue.put(None)

    for p in workers:
        p.join()

    logger.info('indexing table')
    con = MySQLdb.connect(**parse_url(my_url), charset="utf8")
    cur = con.cursor()
    cur.execute(
        """
        CREATE UNIQUE INDEX ui_webfront_protein_identifier
        ON webfront_protein (identifier)
        """
    )
    cur.execute(
        """
        CREATE INDEX i_webfront_protein_database
        ON webfront_protein (source_database)
        """
    )
    cur.execute(
        """
        CREATE INDEX i_webfront_protein_ida
        ON webfront_protein (ida_id)
        """
    )
    cur.execute(
        """
        CREATE INDEX i_webfront_protein_fragment
        ON webfront_protein (is_fragment)
        """
    )
    cur.close()
    con.close()

    logger.info("complete")


def insert_isoforms(my_url: str, ora_ippro_url: str):
    entries = get_entries(my_url)
    query = """
        INSERT INTO webfront_varsplic (accession, protein_acc, length,
                                       sequence, features)
        VALUES (%s, %s, %s, %s, %s)
    """

    con = MySQLdb.connect(**parse_url(my_url), charset="utf8")
    with Table(con, query) as table:
        for isoform in oracle.get_isoforms(ora_ippro_url):
            features = {}
            to_condense = {}
            for signature_acc, locations in isoform["features"].items():
                e = entries[signature_acc]
                entry_acc = e["integrated"]
                features[signature_acc] = {
                    "accession": signature_acc,
                    "integrated": entry_acc,
                    "locations": locations,
                    "name": e["name"],
                    "type": e["type"],
                    "source_database": e["database"]
                }

                if entry_acc is None:
                    continue

                try:
                    entry = to_condense[entry_acc]
                except KeyError:
                    entry = to_condense[entry_acc] = []
                finally:
                    entry += [l["fragments"] for l in locations]

            for entry_acc in to_condense:
                locations = []
                for start, end in condense_locations(to_condense[entry_acc]):
                    locations.append({
                        "fragments": [{
                            "start": start,
                            "end": end,
                            "dc-status": "CONTINUOUS"
                        }],
                        "model_acc": None
                    })

                features[entry_acc] = {
                    "accession": entry_acc,
                    "integrated": None,
                    "locations": locations,
                    "name": entries[entry_acc]["name"],
                    "type": entries[entry_acc]["type"],
                    "source_database": entries[entry_acc]["database"]
                }

            table.insert((
                isoform["accession"],
                isoform["protein_acc"],
                isoform["length"],
                isoform["sequence"],
                json.dumps(features)
            ))

    con.commit()

    cur = con.cursor()
    cur.execute("CREATE INDEX i_webfront_varsplic "
                "ON webfront_varsplic (protein_acc)")
    cur.close()
    con.close()
