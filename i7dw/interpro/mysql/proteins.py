# -*- coding: utf-8 -*-

import hashlib
import json
from typing import Dict, List, Optional, Tuple

import MySQLdb

from i7dw import io, logger, pdbe
from i7dw.interpro import Populator, oracle, extract_frag, extract_loc, MIN_OVERLAP
from .entries import get_entries, iter_sets
from .structures import iter_structures
from .taxonomy import get_taxa
from .utils import parse_url


def calculate_domain_architecture(pfam_entries: List[Tuple]) -> Tuple[str, str]:
    """
    Tuple:
      Pfam acc: str, InterPro acc: Optional[str], locations: List[Dict]

    Sort Pfam entries based on the leftmost fragment of the leftmost location
        - e[2]          -> access locations
        - e[2][0]       -> access leftmost location
        - e[2][0][0]    -> access left fragment
    """
    pfam_entries.sort(key=lambda e: extract_frag(e[2][0][0]))

    dom_arch = []
    for pfam_acc, interpro_acc, locations in pfam_entries:
        if interpro_acc:
            dom_arch.append(pfam_acc + ':' + interpro_acc)
        else:
            dom_arch.append(pfam_acc)

    dom_arch = '-'.join(dom_arch)
    dom_arch_id = hashlib.sha1(dom_arch.encode("utf-8")).hexdigest()
    return dom_arch, dom_arch_id
    

def insert_proteins(my_url: str, ora_ippro_url: str, ora_pdbe_url: str,
                    src_proteins: str, src_sequences: str, src_misc: str,
                    src_names: str, src_comments: str, src_proteomes: str,
                    src_residues: str, src_features: str, src_matches: str,
                    limit: int=0, tmpdir: Optional[str]=None):
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

            go_terms = {}  # InterPro2GO + InterPro matches -> UniProt-GO
            protein2entries = {}
            protein2sets = set()
            pfam_entries = []
            for entry_acc, locations in protein2matches.get(protein_acc, {}).items():
                e = entries[entry_acc]
                database = e["database"]
                try:
                    protein2entries[database] += 1
                except KeyError:
                    protein2entries[database] = 1

                if database == "pfam":
                    pfam_entries.append((entry_acc, e["integrated"], locations))

                for term in e["go_terms"]:
                    go_terms[term["identifier"]] = term

                try:
                    set_acc = entry2set[entry_acc]
                except KeyError:
                    pass
                else:
                    protein2sets.add(set_acc)

            protein_entries["total"] = sum(protein_entries.values())

            structures = {
                "feature": {
                    "cath": {},
                    "scop": {}
                },
                "prediction": protein2predictions.get(protein_acc, {})
            }
            cath_features = {}
            scop_features = {}
            for pdbe_id in protein2structures.get(protein_acc, []):
                for dom_id, dom in cath_domains.get(pdbe_id, {}).items():
                    structures["feature"]["cath"][dom_id] = dom

                for dom_id, dom in scop_domains.get(pdbe_id, {}).items():
                    structures["feature"]["scop"][dom_id] = dom

            dom_arch = protein2ida.get(acc)
            if dom_arch:
                ida, ida_id = dom_arch
                num_ida = ida_counts[ida_id]
            else:
                ida = ida_id = None
                num_ida = 0

            # Enqueue record for protein table
            table.insert((
                protein_acc,
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
                    e = loc["fragments"][-1]["end"]

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
