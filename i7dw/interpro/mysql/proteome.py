import os
import json
from typing import Optional

from . import entry, structure, taxonomy, reduce
from ... import dbms, logger, uniprot
from ...io import Store


def insert_proteomes(ora_uri, my_uri, chunk_size=100000):
    proteomes = uniprot.get_proteomes(ora_uri)
    taxa = set(taxonomy.get_taxa(my_uri, lineage=False))

    data = []
    con, cur = dbms.connect(my_uri)
    for p in proteomes.values():
        if p["tax_id"] not in taxa:
            """
            If tax_id not in taxa, it's very likely that INTERPRO.ETAXI
            (source for taxonomy table) is out-of-date
            """
            logger.warning("missing taxon (ID: {})".format(p["tax_id"]))
            continue

        data.append((
            p["accession"],
            p["name"],
            1 if p["is_reference"] else 0,
            p["strain"],
            p["assembly"],
            p["tax_id"]
        ))

    for i in range(0, len(data), chunk_size):
        cur.executemany(
            """
            INSERT INTO webfront_proteome (
              accession,
              name,
              is_reference,
              strain,
              assembly,
              taxonomy_id
            ) VALUES (
              %s, %s, %s, %s, %s, %s
            )
            """,
            data[i:i+chunk_size]
        )

    cur.close()
    con.commit()
    con.close()


def get_proteomes(uri: str) -> dict:
    con, cur = dbms.connect(uri)

    cur.execute(
        """
        SELECT accession, name, is_reference, strain, assembly, taxonomy_id
        FROM webfront_proteome
        """
    )

    proteomes = {}
    for row in cur:
        proteomes[row[0]] = {
            'name': row[1],
            'is_reference': bool(row[2]),
            'strain': row[3],
            'assembly': row[4],
            'taxon': row[5]
        }

    cur.close()
    con.close()

    return proteomes


def update_counts(my_uri: str, src_proteins: str, src_proteomes:str,
                  src_matches: str, src_ida: str, processes: int=1,
                  tmpdir: Optional[str]=None, sync_frequency: int=100000):
    logger.info("starting")
    if tmpdir:
        os.makedirs(tmpdir, exist_ok=True)

    # Open existing stores containing protein-related info
    proteins = Store(src_proteins)
    protein2proteome = Store(src_proteomes)
    protein2matches = Store(src_matches)
    protein2ida = Store(src_ida)

    # Get required MySQL data
    entries = entry.get_entries(my_uri)
    protein2structures = structure.get_protein2structures(my_uri)
    entry2set = entry.get_entry2set(my_uri)
    proteomes = get_proteomes(my_uri)

    protein_counts = {}
    cnt_proteins = 0
    cnt_updates = 0
    with Store(keys=Store.chunk_keys(proteomes, 100), tmpdir=tmpdir) as xrefs:
        for protein_acc, p in proteins:
            cnt_proteins += 1
            if not cnt_proteins % 10000000:
                logger.info(f"{cnt_proteins:>12}")

            tax_id = p["taxon"]
            try:
                upid = protein2proteome[protein_acc]
            except KeyError:
                continue

            prot_entries = set()
            for match in protein2matches.get(protein_acc, []):
                method_acc = match["method_ac"]
                prot_entries.add(method_acc)

                entry_acc = entries[method_acc]["integrated"]
                if entry_acc:
                    prot_entries.add(entry_acc)

            entry_databases = {}
            entry_sets = set()
            for entry_acc in prot_entries:
                database = entries[entry_acc]["database"]
                if database in entry_databases:
                    entry_databases[database].add(entry_acc)
                else:
                    entry_databases[database] = {entry_acc}

                try:
                    set_acc = entry2set[entry_acc]
                except KeyError:
                    pass
                else:
                    entry_sets.add(set_acc)

            _xrefs = {
                "entries": entry_databases,
                "sets": entry_sets,
                "taxa": {tax_id}
            }
            try:
                ida, ida_id = protein2ida[protein_acc]
            except KeyError:
                pass
            else:
                _xrefs["domain_architectures"] = {ida}

            try:
                pdbe_ids = protein2structures[protein_acc]
            except KeyError:
                pass
            else:
                _xrefs["structures"] = pdbe_ids

            xrefs.update(upid, _xrefs)
            cnt_updates += 1
            if not cnt_updates % sync_frequency:
                xrefs.sync()

            if upid in protein_counts:
                protein_counts[upid] += 1
            else:
                protein_counts[upid] = 1

        proteins.close()
        protein2proteome.close()
        protein2matches.close()
        protein2ida.close()
        logger.info(f"{cnt_proteins:>12}")

        for upid, cnt in protein_counts.items():
            proteomes.pop(upid)
            xrefs.update(upid, {"proteins": cnt})

        # Remaining proteomes
        for upid in proteomes:
            xrefs.update(upid, {
                "proteins": [],
                "domain_architectures": [],
                "structures": [],
                "entries": {},
                "sets": [],
                "taxa": {proteomes[upid]["taxon"]},
            })

        size = xrefs.merge(processes=processes)
        logger.info("Disk usage: {:.0f}MB".format(size/1024**2))

        con, cur = dbms.connect(my_uri)
        cur.close()
        query = "UPDATE webfront_proteome SET counts = %s WHERE accession = %s"
        with dbms.Populator(con, query) as table:
            for upid, _xrefs in xrefs:
                counts = reduce(_xrefs)
                counts["entries"]["total"] = sum(counts["entries"].values())
                table.update((json.dumps(counts), upid))

        con.commit()
        con.close()

    logger.info("complete")
