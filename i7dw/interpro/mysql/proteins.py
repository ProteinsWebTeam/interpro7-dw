# -*- coding: utf-8 -*-

import json
from typing import Optional

import MySQLdb

from i7dw import io, logger, pdbe
from i7dw.interpro import condense_locations
from i7dw.interpro import DomainArchitecture, Table
from i7dw.interpro import oracle
from .entries import get_entries, iter_sets
from .structures import iter_structures
from .taxonomy import iter_taxa
from .utils import drop_index, parse_url


def insert_proteins(my_url: str, ora_ippro_url: str, ora_pdbe_url: str,
                    src_proteins: str, src_sequences: str, src_misc: str,
                    src_names: str, src_comments: str, src_proteomes: str,
                    src_residues: str, src_features: str, src_matches: str,
                    processes: int=1, tmpdir: Optional[str]=None):
    logger.info("calculating domain architectures")
    proteins = io.Store(src_proteins)
    matches = io.Store(src_matches)
    domains = io.Store(keys=matches.keys, tmpdir=tmpdir)
    dom_arch = DomainArchitecture(get_entries(my_url))
    dom_cnts = {}
    n_proteins = 0
    for protein_acc, protein_info in proteins:
        dom_arch.update(matches.get(protein_acc, {}))
        dom_ac = dom_arch.accession
        dom_id = dom_arch.identifier
        domains[protein_acc] = (dom_ac, dom_id)

        try:
            dom_cnts[dom_id] += 1
        except KeyError:
            dom_cnts[dom_id] = 1

        n_proteins += 1
        if not n_proteins % 1000000:
            domains.sync()

    size = domains.merge(processes=processes)
    logger.info(f"  temporary files: {size/1024/1024:.0f} MB")

    logger.info("loading data")
    # Structural features (CATH and SCOP domains)
    cath_domains = pdbe.get_cath_domains(ora_pdbe_url)
    scop_domains = pdbe.get_scop_domains(ora_pdbe_url)

    # Structural predictions (ModBase and Swiss-Model models)
    predictions = oracle.get_structural_predictions(ora_ippro_url)

    # MySQL data
    entries = get_entries(my_url)
    taxa = {}
    for taxon in iter_taxa(my_url, lineage=False):
        taxa[taxon["id"]] = {
            "taxId": taxon["id"],
            "scientificName": taxon["scientific_name"],
            "fullName": taxon["full_name"]
        }

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

    logger.info("inserting proteins")
    query = """
        INSERT INTO webfront_protein (accession, identifier, organism, name,
                                      other_names, description, sequence,
                                      length, size, proteome, gene, go_terms,
                                      evidence_code, source_database,
                                      residues, is_fragment, structure,
                                      tax_id, extra_features, ida_id, ida,
                                      counts)
        VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s,
                %s, %s, %s, %s, %s, %s, %s)
    """

    with Table(con, query) as table:
        sequences = io.Store(src_sequences)
        misc = io.Store(src_misc)
        names = io.Store(src_names)
        comments = io.Store(src_comments)
        proteomes = io.Store(src_proteomes)
        residues = io.Store(src_residues)
        features = io.Store(src_features)
        n_proteins = 0

        for protein_acc, protein_info in proteins:
            tax_id = protein_info["taxon"]
            try:
                taxon = taxa[tax_id]
            except KeyError:
                continue

            try:
                sequence = sequences[protein_acc]
            except KeyError:
                continue

            try:
                evidence, gene = misc[protein_acc]
            except KeyError:
                continue

            if not evidence:
                continue

            if protein_info["length"] <= 100:
                size = "small"
            elif protein_info["length"] <= 1000:
                size = "medium"
            else:
                size = "large"

            go_terms = {}  # InterPro2GO + InterPro matches -> UniProt-GO
            protein_entries = {}
            protein_sets = set()
            for entry_acc in matches.get(protein_acc, {}):
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
                "feature": {"cath": {}, "scop": {}},
                "prediction": predictions.get(protein_acc, {})
            }
            for pdbe_id in structures.get(protein_acc, []):
                for dom_id, dom in cath_domains.get(pdbe_id, {}).items():
                    protein_structures["feature"]["cath"][dom_id] = dom

                for dom_id, dom in scop_domains.get(pdbe_id, {}).items():
                    protein_structures["feature"]["scop"][dom_id] = dom

            name, other_names = names.get(protein_acc, (None, None))
            upid = proteomes.get(protein_acc)

            dom_ac, dom_id = domains[protein_acc]
            dom_cnt = dom_cnts[dom_id]

            # Enqueue record for protein table
            table.insert((
                protein_acc,
                protein_info["identifier"],
                json.dumps(taxon),
                name,
                json.dumps(other_names),
                json.dumps(comments.get(protein_acc, [])),
                sequence,
                protein_info["length"],
                size,
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

            n_proteins += 1
            if not n_proteins % 10000000:
                logger.info('{:>12,}'.format(n_proteins))

    logger.info('{:>12,}'.format(n_proteins))
    table.close()
    con.commit()

    logger.info('indexing/analyzing table')
    cur = con.cursor()
    cur.execute(
        """
        CREATE UNIQUE INDEX ui_webfront_protein_identifier
        ON webfront_protein (identifier)
        """
    )
    cur.execute(
        """
        CREATE INDEX i_webfront_protein_length
        ON webfront_protein (length)
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
    cur.execute("ANALYZE TABLE webfront_protein")
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
