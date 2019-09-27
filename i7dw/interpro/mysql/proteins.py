# -*- coding: utf-8 -*-

import json

import MySQLdb

from i7dw import io, logger, pdbe
from i7dw.interpro import Populator, oracle, extract_loc, MIN_OVERLAP
from .entries import get_entries, iter_sets
from .structures import iter_structures
from .taxonomy import get_taxa
from .utils import parse_url


def insert_proteins(my_url: str, ora_ippro_url: str, ora_pdbe_url: str,
                    src_proteins: str, src_sequences: str, src_misc: str,
                    src_names: str, src_comments: str, src_proteomes: str,
                    src_residues: str, src_features: str, src_matches: str,
                    **kwargs):
    limit = kwargs.get("limit", 0)

    logger.info("loading data")

    # Structural features (CATH and SCOP domains)
    cath_domains = pdbe.get_cath_domains(ora_pdbe_url)
    scop_domains = pdbe.get_scop_domains(ora_pdbe_url)

    # Structural predictions (ModBase and Swiss-Model models)
    protein2predictions = oracle.get_structural_predictions(ora_ippro_url)

    # MySQL data
    taxa = get_taxa(my_url, lineage=False)
    entries = get_entries(my_url)

    protein2structures = {}
    for s in iter_structures(my_url):
        pdbe_id = s["accession"]
        for protein_acc in s["proteins"]:
            try:
                protein2structures[protein_acc].add(pdbe_id)
            except KeyError:
                protein2structures[protein_acc] = {pdbe_id}

    entry2set = {}
    for s in iter_sets(my_url):
        set_acc = s["accession"]
        for entry_acc in s["members"]:
            entry2set[entry_acc] = set_acc

    proteins = io.Store(src_proteins)
    protein2sequence = io.Store(src_sequences)
    protein2misc = io.Store(src_misc)
    protein2names = io.Store(src_names)
    protein2comments = io.Store(src_comments)
    protein2proteome = io.Store(src_proteomes)
    protein2residues = io.Store(src_residues)
    protein2features = io.Store(src_features)
    protein2matches = io.Store(src_matches)

    logger.info("truncating table")
    con = MySQLdb.connect(**parse_url(my_url), charset="utf8")
    cur = con.cursor()
    cur.execute("TRUNCATE TABLE webfront_protein")
    for index in ("ui_webfront_protein_identifier",
                  "i_webfront_protein_length",
                  "i_webfront_protein_ida"):
        try:
            cur.execute(f"DROP INDEX {index} ON webfront_protein")
        except Exception:
            pass

    logger.info("counting isoforms")
    cur.execute(
        """
        SELECT protein_acc, COUNT(*)
        FROM webfront_varsplic
        GROUP BY protein_acc
        """
    )
    isoforms = dict(cur.fetchall())
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
                %s, %s, %s, %s, %s,%s, %s)    
    """

    with Populator(con, query) as table:
        n_proteins = 0

        for protein_acc, protein in proteins:
            tax_id = protein["taxon"]
            try:
                taxon = taxa[tax_id]
            except KeyError:
                continue

            try:
                sequence = protein2sequence[protein_acc]
            except KeyError:
                continue

            try:
                evidence, gene = protein2misc[protein_acc]
            except KeyError:
                continue

            if not evidence:
                continue

            if protein["length"] <= 100:
                size = "small"
            elif protein["length"] <= 1000:
                size = "medium"
            else:
                size = "large"

            name, other_names = protein2names.get(protein_acc, (None, None))
            upid = protein2proteome.get(protein_acc)

            # InterPro2GO + InterPro matches -> UniProt-GO
            go_terms = {}
            _entries = set()
            for m in protein2matches.get(acc, []):
                method_ac = m["method_ac"]
                _entries.add(method_ac)

                if method_ac in integrated:
                    entry_ac = integrated[method_ac]

                    _entries.add(entry_ac)
                    for term in entries[entry_ac]["go_terms"]:
                        go_terms[term["identifier"]] = term

            protein2entries = {"total": len(_entries)}
            protein2sets = set()
            for entry_ac in _entries:
                db_name = entries[entry_ac]["database"]
                if db_name in protein2entries:
                    protein2entries[db_name] += 1
                else:
                    protein2entries[db_name] = 1

                if entry_ac in entry2set:
                    protein2sets.add(entry2set[entry_ac])

            cath_features = {}
            scop_features = {}
            for pdb_id in protein2pdb.get(acc, []):
                for domain_id in cath_domains.get(pdb_id, {}):
                    cath_features[domain_id] = cath_domains[pdb_id][domain_id]

                for domain_id in scop_domains.get(pdb_id, {}):
                    scop_features[domain_id] = scop_domains[pdb_id][domain_id]

            structures = {}
            if cath_features or scop_features:
                structures["feature"] = {}

                if cath_features:
                    structures["feature"]["cath"] = cath_features

                if scop_features:
                    structures["feature"]["scop"] = scop_features

            if acc in protein2predictions:
                structures["prediction"] = protein2predictions[acc]

            dom_arch = protein2ida.get(acc)
            if dom_arch:
                ida, ida_id = dom_arch
                num_ida = ida_counts[ida_id]
            else:
                ida = ida_id = None
                num_ida = 0

            # Enqueue record for protein table
            table.insert((
                acc,
                protein["identifier"],
                json.dumps(taxon),
                name,
                json.dumps(other_names),
                json.dumps(protein2comments.get(acc, [])),
                sequence,
                protein["length"],
                size,
                upid,
                gene,
                json.dumps(list(go_terms.values())),
                evidence,
                "reviewed" if protein["is_reviewed"] else "unreviewed",
                json.dumps(protein2residues.get(acc, {})),
                1 if protein["is_fragment"] else 0,
                json.dumps(structures),
                tax_id,
                json.dumps(protein2features.get(acc, {})),
                ida_id,
                ida,
                json.dumps({
                    "entries": protein2entries,
                    "structures": len(protein2pdb.get(acc, [])),
                    "sets": len(protein2sets),
                    "proteomes": 1 if upid else 0,
                    "taxa": 1,
                    "idas": num_ida,
                    "isoforms": isoforms.get(acc, 0)
                })
            ))

            n_proteins += 1
            if n_proteins == limit:
                break
            elif not n_proteins % 10000000:
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
    with Populator(con, query) as table:
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

                entry += [l["fragments"] for l in locations]

            for entry_acc in to_condense:
                start = end = None
                locations = []

                # Sort locations using their leftmost fragment
                for loc in sorted(to_condense[entry_acc], key=extract_loc):
                    # We do not consider fragmented matches
                    s = loc["fragments"][0]["start"]
                    e = loc["fragments"][0]["start"]

                    if start is None:
                        start, end = s, e
                        continue
                    elif s <= end:
                        # Locations are overlapping (at least one residue)
                        overlap = min(end, e) - max(start, s) + 1
                        shortest = min(end - start, e - s) + 1

                        if overlap >= shortest * MIN_OVERLAP:
                            # Merge
                            end = e
                            continue

                    locations.append((start, end))
                    start, end = s, e

                # Adding last location
                locations.append((start, end))

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

