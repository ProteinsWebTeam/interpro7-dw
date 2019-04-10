import json
import time

from . import entry, structure, taxonomy
from .. import oracle
from ... import dbms, logger, pdbe
from ...io import Store


def insert(ora_ippro_uri: str, ora_pdbe_uri: str, my_uri: str,
           src_proteins: str, src_sequences: str, src_misc: str,
           src_names: str, src_comments: str, src_proteomes: str,
           src_residues: str, src_features: str, src_matches: str,
           src_idas: str, **kwargs):
    chunk_size = kwargs.get("chunk_size", 100000)
    limit = kwargs.get("limit", 0)

    logger.info("starting")

    # Structural features (CATH and SCOP domains)
    cath_domains = pdbe.get_cath_domains(ora_pdbe_uri)
    scop_domains = pdbe.get_scop_domains(ora_pdbe_uri)

    # Structural predictions (ModBase and Swiss-Model models)
    protein2predictions = oracle.get_structural_predictions(ora_ippro_uri)

    # MySQL data
    taxa = taxonomy.get_taxa(my_uri, lineage=False)
    entries = entry.get_entries(my_uri)
    integrated = {
        acc: e["integrated"]
        for acc, e in entries.items()
        if e["integrated"]
    }

    protein2pdb = {}
    for pdb_id, s in structure.get_structures(my_uri).items():
        for acc in s["proteins"]:
            if acc in protein2pdb:
                protein2pdb[acc].add(pdb_id)
            else:
                protein2pdb[acc] = {pdb_id}

    entry2set = {}
    for set_ac, s in entry.get_sets(my_uri).items():
        for acc in s["members"]:
            entry2set[acc] = set_ac

    proteins = Store(src_proteins)
    protein2sequence = Store(src_sequences)
    protein2misc = Store(src_misc)
    protein2names = Store(src_names)
    protein2comments = Store(src_comments)
    protein2proteome = Store(src_proteomes)
    protein2residues = Store(src_residues)
    protein2features = Store(src_features)
    protein2matches = Store(src_matches)
    protein2ida = Store(src_idas)

    ida_counts = {}
    for acc, dom_arch in protein2ida:
        ida, ida_id = dom_arch
        if ida_id in ida_counts:
            ida_counts[ida_id] += 1
        else:
            ida_counts[ida_id] = 1

    con, cur = dbms.connect(my_uri)

    # Truncate table *only* if no specific accessions are passed
    cur.execute("TRUNCATE TABLE webfront_protein")
    for index in ("ui_webfront_protein_identifier",
                  "i_webfront_protein_length",
                  "i_webfront_protein_ida"):
        try:
            cur.execute("DROP INDEX {} ON webfront_protein".format(index))
        except Exception:
            pass

    logger.info("inserting proteins")
    query = """
        INSERT INTO webfront_protein (
            accession, identifier, organism, name, other_names,
            description, sequence, length, size, proteome, 
            gene, go_terms, evidence_code, source_database, residues, 
            is_fragment, structure, tax_id, extra_features, ida_id, 
            ida, counts
        )
        VALUES (
            %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, 
            %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, 
            %s, %s
        )    
    """
    data = []
    n_proteins = 0
    ts = time.time()
    for acc, protein in proteins:
        tax_id = protein["taxon"]
        try:
            taxon = taxa[tax_id]
        except KeyError:
            continue

        try:
            sequence = protein2sequence[acc]
        except KeyError:
            continue

        evidence, gene = protein2misc.get(acc, (None, None))
        if not evidence:
            continue

        if protein["length"] <= 100:
            size = "small"
        elif protein["length"] <= 1000:
            size = "medium"
        else:
            size = "large"

        name, other_names = protein2names.get(acc, (None, None))
        upid = protein2proteome.get(acc)

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
        data.append((
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
            "reviewed" if protein["isReviewed"] else "unreviewed",
            json.dumps(protein2residues.get(acc, {})),
            1 if protein["isFrag"] else 0,
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
                "idas": num_ida
            })
        ))

        if len(data) == chunk_size:
            cur.executemany(query, data)
            data = []

        n_proteins += 1
        if n_proteins == limit:
            break
        elif not n_proteins % 10000000:
            logger.info('{:>12,} ({:.0f} proteins/sec)'.format(
                n_proteins, n_proteins / (time.time() - ts)
            ))

    if data:
        cur.executemany(query, data)
    con.commit()

    logger.info('{:>12,} ({:.0f} proteins/sec)'.format(
        n_proteins, n_proteins / (time.time() - ts)
    ))

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
    cur.execute("ANALYZE TABLE webfront_protein")
    cur.close()
    con.close()

    logger.info("complete")