import os
import time
from datetime import datetime
from multiprocessing import Process, Queue
from typing import Dict, Optional, Set

from . import mysql
from .. import dbms, logger, pdbe
from ..io import Store


def create_store(store: Store, task_queue: Queue, done_queue: Queue):
    while True:
        chunk = task_queue.get()
        if chunk is None:
            break

        for args in chunk:
            store.update_from_seq(*args)

        store.sync()

    store.dir = None            # prevent temporary directory to be removed
    done_queue.put(store.type)  # pass type to parent process


def chunk_keys(keys: list, chunk_size: int) -> list:
    return [keys[i] for i in range(0, len(keys), chunk_size)]


def export(my_uri: str, src_proteins: str, src_matches: str,
           src_proteomes: str, dst_entries: str, dst_proteomes: str,
           dst_structures: str, dst_taxa: str, flush: int=100000,
           tmpdir: str=None, limit: int=0):
    logger.info("starting")

    """
    We do not keep track of set ---> entities
    because these can be found from the union of the set's members
    """

    if tmpdir is not None:
        os.makedirs(tmpdir, exist_ok=True)

    # Get entries
    entry_database = {}
    integrated = {}
    for acc, e in mysql.entry.get_entries(my_uri).items():
        entry_database[acc] = e["database"]

        if e["integrated"]:
            integrated[acc] = e["integrated"]

    # Creating (single-threaded) stores
    entries_store = Store(dst_entries,
                          keys=chunk_keys(sorted(entry_database), 10),
                          tmpdir=tmpdir)
    proteomes_store = Store(dst_proteomes,
                            keys=chunk_keys(
                                sorted(mysql.proteome.get_proteomes(my_uri)),
                                100
                            ),
                            tmpdir=tmpdir)
    structures_store = Store(dst_structures,
                             keys=chunk_keys(sorted(
                                 mysql.structure.get_structures(my_uri)), 100
                             ),
                             tmpdir=tmpdir)
    taxa_store = Store(dst_taxa,
                       keys=chunk_keys(sorted(
                           mysql.taxonomy.get_taxa(my_uri)), 10
                       ),
                       tmpdir=tmpdir)

    done_queue = Queue()

    # Start child processes in which stores will be populated
    entries_data = []
    entries_queue = Queue(maxsize=1)
    entries_proc = Process(target=create_store,
                           args=(entries_store, entries_queue, done_queue))
    entries_proc.start()

    proteomes_data = []
    proteomes_queue = Queue(maxsize=1)
    proteomes_proc = Process(target=create_store,
                             args=(proteomes_store, proteomes_queue, done_queue))
    proteomes_proc.start()

    structures_data = []
    structures_queue = Queue(maxsize=1)
    structures_proc = Process(target=create_store,
                              args=(structures_store, structures_queue, done_queue))
    structures_proc.start()

    taxa_data = []
    taxa_queue = Queue(maxsize=1)
    taxa_proc = Process(target=create_store,
                        args=(taxa_store, taxa_queue, done_queue))
    taxa_proc.start()

    # Set members
    entry_set = {
        entry_ac: set_ac
        for set_ac, s in mysql.entry.get_sets(my_uri).items()
        for entry_ac in s["members"]
    }

    # Protein -> PDBe structure
    protein2pdb = {}
    for pdb_id, s in mysql.structure.get_structures(my_uri).items():
        for protein_ac in s["proteins"]:
            if protein_ac in protein2pdb:
                protein2pdb[protein_ac].add(pdb_id)
            else:
                protein2pdb[protein_ac] = {pdb_id}

    # Open existing stores containing protein-related info
    proteins = Store(src_proteins)
    protein2matches = Store(src_matches)
    protein2proteome = Store(src_proteomes)

    taxon2proteins = {}
    proteome2proteins = {}
    structure2proteins = {}
    entry2matches = {}

    n_proteins = 0
    ts = time.time()
    for protein_ac, protein in proteins:
        tax_id = protein["taxon"]

        protein_id = protein["identifier"]
        matches = protein2matches.get(protein_ac, [])
        upid = protein2proteome.get(protein_ac)
        protein_structures = protein2pdb.get(protein_ac, [])

        # Create domain architecture, and count protein matches
        protein_entries = {}
        dom_entries = set()
        dom_arch = []
        for m in matches:
            method_ac = m["method_ac"]
            entry_ac = integrated.get(method_ac)

            if method_ac in protein_entries:
                protein_entries[method_ac] += 1
            else:
                protein_entries[method_ac] = 1

            if entry_ac:
                if entry_ac in protein_entries:
                    protein_entries[entry_ac] += 1
                else:
                    protein_entries[entry_ac] = 1

            if entry_database[method_ac] == "pfam":
                dom_entries.add(method_ac)
                if entry_ac:
                    dom_entries.add(entry_ac)
                    dom_arch.append("{}:{}".format(method_ac, entry_ac))
                else:
                    dom_arch.append("{}".format(method_ac))

        if dom_arch:
            dom_arch = '-'.join(dom_arch)
        else:
            dom_arch = None

        # update matches and add domain architectures to entries
        protein_sets = set()
        for entry_ac, n_matches in protein_entries.items():
            if entry_ac in entry2matches:
                entry2matches[entry_ac] += n_matches
            else:
                entry2matches[entry_ac] = n_matches

            if entry_ac in dom_entries:
                # Has a domain architecture
                entries_data.append((entry_ac, "domain_architectures", dom_arch))

            if entry_ac in entry_set:
                protein_sets.add(entry_set[entry_ac])

        # Taxon ---> protein
        if tax_id in taxon2proteins:
            taxon2proteins[tax_id] += 1
        else:
            taxon2proteins[tax_id] = 1

        if upid:
            # Proteome ---> protein
            if upid in proteome2proteins:
                proteome2proteins[upid] += 1
            else:
                proteome2proteins[upid] = 1

            # Proteome <---> taxon
            proteomes_data.append((upid, "taxa", tax_id))
            taxa_data.append((tax_id, "proteomes", upid))

            # ---> Domain architecture
            if dom_arch:
                proteomes_data.append((upid, "domain_architectures", dom_arch))
                taxa_data.append((tax_id, "domain_architectures", dom_arch))

            for entry_ac in protein_entries:
                database = entry_database[entry_ac]

                # Proteome <---> entries
                proteomes_data.append((upid, "entries", database, entry_ac))
                entries_data.append((entry_ac, "proteomes", upid))

                # Entry ---> protein
                entries_data.append((entry_ac, "proteins", (protein_ac,
                                                            protein_id)))

                # Entry <---> taxon
                entries_data.append((entry_ac, "taxa", tax_id))
                taxa_data.append((tax_id, "entries", database, entry_ac))

                # Entry <---> structure
                for pdb_id in protein_structures:
                    entries_data.append((entry_ac, "structures", pdb_id))
                    structures_data.append((pdb_id, "entries", database,
                                            entry_ac))

            for set_ac in protein_sets:
                # Proteome ---> set
                proteomes_data.append((upid, "sets", set_ac))

                # Taxon --> set
                taxa_data.append((tax_id, "sets", set_ac))

            for pdb_id in protein_structures:
                # Structure ---> protein
                if pdb_id in structure2proteins:
                    structure2proteins[pdb_id] += 1
                else:
                    structure2proteins[pdb_id] = 1

                # Structure <---> taxon
                structures_data.append((pdb_id, "taxa", tax_id))
                taxa_data.append((tax_id, "structures", pdb_id))

                # Structure <---> proteome
                proteomes_data.append((upid, "structures", pdb_id))
                structures_data.append((pdb_id, "proteomes", upid))

                if dom_arch:
                    structures_data.append((pdb_id, "domain_architectures", dom_arch))

                # Structure ---> set
                for set_ac in protein_sets:
                    structures_data.append((pdb_id, "sets", set_ac))
        else:
            # ---> Domain architecture
            if dom_arch:
                taxa_data.append((tax_id, "domain_architectures", dom_arch))

            for entry_ac in protein_entries:
                database = entry_database[entry_ac]

                # Entry ---> protein
                entries_data.append((entry_ac, "proteins", (protein_ac,
                                                            protein_id)))

                # Entry <---> taxon
                entries_data.append((entry_ac, "taxa", tax_id))
                taxa_data.append((tax_id, "entries", database, entry_ac))

                # Entry <---> structure
                for pdb_id in protein_structures:
                    entries_data.append((entry_ac, "structures", pdb_id))
                    structures_data.append((pdb_id, "entries", database,
                                            entry_ac))

            for set_ac in protein_sets:
                # Taxon --> set
                taxa_data.append((tax_id, "sets", set_ac))

            for pdb_id in protein_structures:
                # Structure ---> protein
                if pdb_id in structure2proteins:
                    structure2proteins[pdb_id] += 1
                else:
                    structure2proteins[pdb_id] = 1

                # Structure <---> taxon
                structures_data.append((pdb_id, "taxa", tax_id))
                taxa_data.append((tax_id, "structures", pdb_id))

                if dom_arch:
                    structures_data.append((pdb_id, "domain_architectures", dom_arch))

                # Structure ---> set
                for set_ac in protein_sets:
                    structures_data.append((pdb_id, "sets", set_ac))

        n_proteins += 1
        if not n_proteins % flush:
            entries_queue.put(entries_data)
            entries_data = []
            proteomes_queue.put(proteomes_data)
            proteomes_data = []
            structures_queue.put(structures_data)
            structures_data = []
            taxa_queue.put(taxa_data)
            taxa_data = []

        if n_proteins == limit:
            break
        elif not n_proteins % 10000000:
            logger.info('{:>12,} ({:.0f} proteins/sec)'.format(
                n_proteins, n_proteins / (time.time() - ts)
            ))

    proteins.close()
    protein2matches.close()
    protein2proteome.close()

    # Sending remaining items to child processes
    for entry_ac, cnt in entry2matches.items():
        entries_data.append((entry_ac, "matches", cnt))
    entries_queue.put(entries_data)
    entries_data = None
    entries_queue.put(None)

    for upid, cnt in proteome2proteins.items():
        proteomes_data.append((upid, "proteins", cnt))
    proteomes_queue.put(proteomes_data)
    proteomes_data = None
    proteomes_queue.put(None)

    for pdb_id, cnt in structure2proteins.items():
        structures_data.append((pdb_id, "proteins", cnt))
    structures_queue.put(structures_data)
    structures_data = None
    structures_queue.put(None)

    for tax_id, cnt in taxon2proteins.items():
        taxa_data.append((tax_id, "proteins", cnt))
    taxa_queue.put(taxa_data)
    taxa_data = None
    taxa_queue.put(None)

    logger.info('{:>12,} ({:.0f} proteins/sec)'.format(
        n_proteins, n_proteins / (time.time() - ts)
    ))

    """
    Get type from child processes
    (order does not matter: same type for all stores)
    """
    entries_store.type = done_queue.get()
    proteomes_store.type = done_queue.get()
    structures_store.type = done_queue.get()
    taxa_store.type = done_queue.get()

    # Wait for child processes to complete
    entries_proc.join()
    proteomes_proc.join()
    structures_proc.join()
    taxa_proc.join()

    # Merge using multi-processing
    entries_store.merge(processes=4)
    logger.info("temporary files ({}): {:.0f} MB".format(
        os.path.basename(dst_entries), entries_store.size/1024/1024
    ))
    entries_store.close()

    proteomes_store.merge(processes=4)
    logger.info("temporary files ({}): {:.0f} MB".format(
        os.path.basename(dst_proteomes), proteomes_store.size/1024/1024
    ))
    proteomes_store.close()

    structures_store.merge(processes=4)
    logger.info("temporary files ({}): {:.0f} MB".format(
        os.path.basename(dst_proteomes), structures_store.size/1024/1024
    ))
    structures_store.close()

    taxa_store.merge(processes=4)
    logger.info("temporary files ({}): {:.0f} MB".format(
        os.path.basename(dst_taxa), taxa_store.size/1024/1024
    ))
    taxa_store.close()

    logger.info("complete")


class EntryMerger(object):
    def __init__(self, uri: str):
        self.uri = uri
        self.taxa = None

    def merge(self, obj: Dict) -> Dict:
        if self.taxa is None:
            self.taxa = mysql.taxonomy.get_taxa(self.uri, lineage=True)

        taxa = {}
        for tax_id in obj["taxa"]:
            for node_id in self.taxa[tax_id]["lineage"]:
                try:
                    taxa[node_id] += 1
                except KeyError:
                    taxa[node_id] = 1


def export_entries(my_uri: str, src_proteins: str, src_proteomes:str,
                   src_matches: str, src_ida: str, dst_entries: str,
                   processes: int = 1, tmpdir: Optional[str]=None):
    if tmpdir:
        os.makedirs(tmpdir, exist_ok=True)

    # Open existing stores containing protein-related info
    proteins = Store(src_proteins)
    protein2proteome = Store(src_proteomes)
    protein2matches = Store(src_matches)
    protein2ida = Store(src_ida)

    # Get required MySQL data
    entries = mysql.entry.get_entries(my_uri)
    protein2structures = get_protein2structures(my_uri)

    # Open new store
    keys = chunk_keys(keys=sorted(entries), chunk_size=100)
    entry_xrefs = Store(dst_entries, keys=keys, tmpdir=tmpdir)

    taxon_protein_counts = {}
    entry_match_counts = {}

    for protein_acc, protein in proteins:
        tax_id = protein["taxon"]
        if tax_id in taxon_protein_counts:
            taxon_protein_counts[tax_id] += 1
        else:
            taxon_protein_counts[tax_id] = 1

        matches = {}
        for match in protein2matches.get(protein_acc, []):
            method_acc = match["method_ac"]
            entry_acc = entries[method_acc]["integrated"]

            if method_acc in matches:
                matches[method_acc] += 1
            else:
                matches[method_acc] = 1

            if entry_acc:
                if entry_acc in matches:
                    matches[entry_acc] += 1
                else:
                    matches[entry_acc] = 1

        obj = {
            "proteins": {(protein_acc, protein["identifier"])},
            "taxa": {tax_id}
        }

        try:
            ida, ida_id = protein2ida[protein_acc]
        except KeyError:
            pass
        else:
            obj["domain_architectures"] = {ida}

        try:
            upid = protein2proteome[protein_acc]
        except KeyError:
            pass
        else:
            obj["proteomes"] = {upid}

        try:
            pdbe_ids = protein2structures[protein_acc]
        except KeyError:
            pass
        else:
            obj["structures"] = pdbe_ids

        for entry_acc, cnt in matches.items():
            entry_xrefs.update(entry_acc, obj)
            if entry_acc in entry_match_counts:
                entry_match_counts[entry_acc] += cnt
            else:
                entry_match_counts[entry_acc] = cnt

    for store in (proteins, protein2proteome, protein2matches, protein2ida):
        store.close()

    sets = mysql.entry.get_sets(my_uri)
    entry2set = {}
    for set_acc, s in sets.items():
        for entry_acc in s["members"]:
            entry2set[entry_acc] = set_acc

    for entry_acc, set_acc in entry2set.items():
        entry_xrefs.update(entry_acc, {"sets": {set_acc}})

    for entry_acc, cnt in entry_match_counts.items():
        entry_xrefs.update(entry_acc, {"matches": cnt})

    entry_xrefs.merge(processes=processes)
    entry_xrefs.close()


def export_taxa(my_uri: str, src_proteins: str, src_proteomes:str,
                src_matches: str, src_ida: str, dst: str,
                processes: int=1, tmpdir: Optional[str]=None,
                sync_frequency: int=100000):
    logger.info("starting")
    if tmpdir:
        os.makedirs(tmpdir, exist_ok=True)

    # Open existing stores containing protein-related info
    proteins = Store(src_proteins)
    protein2proteome = Store(src_proteomes)
    protein2matches = Store(src_matches)
    protein2ida = Store(src_ida)

    # Get required MySQL data
    entries = mysql.entry.get_entries(my_uri)
    protein2structures = get_protein2structures(my_uri)
    lineages = get_lineages(my_uri)

    # Open new store
    keys = chunk_keys(keys=sorted(lineages), chunk_size=10)
    xrefs = Store(dst, keys=keys, tmpdir=tmpdir)

    protein_counts = {}
    cnt_proteins = 0
    for protein_acc, protein in proteins:
        prot_entries = set()
        for match in protein2matches.get(protein_acc, []):
            method_acc = match["method_ac"]
            prot_entries.add(method_acc)

            entry_acc = entries[method_acc]["integrated"]
            if entry_acc:
                prot_entries.add(entry_acc)

        entry_databases = {}
        for entry_acc in prot_entries:
            database = entries[entry_acc]["database"]
            if database in entry_databases:
                entry_databases[database].add(entry_acc)
            else:
                entry_databases[database] = {entry_acc}

        obj = {"entries": entry_databases}

        try:
            ida, ida_id = protein2ida[protein_acc]
        except KeyError:
            pass
        else:
            obj["domain_architectures"] = {ida}

        try:
            upid = protein2proteome[protein_acc]
        except KeyError:
            pass
        else:
            obj["proteomes"] = {upid}

        try:
            pdbe_ids = protein2structures[protein_acc]
        except KeyError:
            pass
        else:
            obj["structures"] = pdbe_ids

        # Associate cross-references to every node in the lineage
        protein_tax_id = protein["taxon"]
        for tax_id in lineages[protein_tax_id]:
            if tax_id in protein_counts:
                protein_counts[tax_id] += 1
            else:
                protein_counts[tax_id] = 1

            xrefs.update(tax_id, obj)

        cnt_proteins += 1
        if not cnt_proteins % sync_frequency:
            xrefs.sync()

        if not cnt_proteins % 10000000:
            logger.info('{:>12,}'.format(cnt_proteins))

    logger.info('{:>12,}'.format(cnt_proteins))

    for store in (proteins, protein2proteome, protein2matches, protein2ida):
        store.close()

    for tax_id, cnt in protein_counts.items():
        xrefs.update(tax_id, {"proteins": cnt})

    xrefs.merge(processes=processes)
    xrefs.close()

    logger.info("complete")


def get_entry2set(my_uri: str) -> Dict[str, str]:
    sets = mysql.entry.get_sets(my_uri)
    entry2set = {}
    for set_acc, s in sets.items():
        for entry_acc in s["members"]:
            entry2set[entry_acc] = set_acc
    return entry2set


def export_proteomes(my_uri: str, src_proteins: str, src_proteomes:str,
                     src_matches: str, src_ida: str, dst: str,
                     processes: int=1, tmpdir: Optional[str]=None,
                     sync_frequency: int=1000000):
    logger.info("starting")
    if tmpdir:
        os.makedirs(tmpdir, exist_ok=True)

    # Open existing stores containing protein-related info
    proteins = Store(src_proteins)
    protein2proteome = Store(src_proteomes)
    protein2matches = Store(src_matches)
    protein2ida = Store(src_ida)

    # Get required MySQL data
    entries = mysql.entry.get_entries(my_uri)
    protein2structures = get_protein2structures(my_uri)
    entry2set = get_entry2set(my_uri)

    # Open new store
    xrefs = Store(filepath=dst,
                  keys=chunk_keys(
                      keys=sorted(mysql.proteome.get_proteomes(my_uri)),
                      chunk_size=100
                  ),
                  tmpdir=tmpdir)

    protein_counts = {}
    cnt_proteins = 0
    for protein_acc, protein in proteins:
        tax_id = protein["taxon"]
        try:
            upid = protein2proteome[protein_acc]
        except KeyError:
            continue

        if upid in protein_counts:
            protein_counts[upid] += 1
        else:
            protein_counts[upid] = 1

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

        obj = {
            "entries": entry_databases,
            "sets": entry_sets,
            "taxa": {tax_id}
        }

        try:
            ida, ida_id = protein2ida[protein_acc]
        except KeyError:
            pass
        else:
            obj["domain_architectures"] = {ida}

        try:
            pdbe_ids = protein2structures[protein_acc]
        except KeyError:
            pass
        else:
            obj["structures"] = pdbe_ids

        xrefs.update(upid, obj)

        cnt_proteins += 1
        if not cnt_proteins % sync_frequency:
            xrefs.sync()

        if not cnt_proteins % 10000000:
            logger.info('{:>12,}'.format(cnt_proteins))

    logger.info('{:>12,}'.format(cnt_proteins))

    for store in (proteins, protein2proteome, protein2matches, protein2ida):
        store.close()

    for upid, cnt in protein_counts.items():
        xrefs.update(upid, {"proteins": cnt})

    xrefs.merge(processes=processes)
    xrefs.close()

    logger.info("complete")


def export_structures(my_uri: str, src_proteins: str, src_proteomes:str,
                      src_matches: str, src_ida: str, dst: str,
                      processes: int=1, tmpdir: Optional[str]=None,
                      sync_frequency: int=1000000):
    logger.info("starting")
    if tmpdir:
        os.makedirs(tmpdir, exist_ok=True)

    # Open existing stores containing protein-related info
    proteins = Store(src_proteins)
    protein2proteome = Store(src_proteomes)
    protein2matches = Store(src_matches)
    protein2ida = Store(src_ida)

    # Get required MySQL data
    entries = mysql.entry.get_entries(my_uri)
    protein2structures = get_protein2structures(my_uri)
    entry2set = get_entry2set(my_uri)

    # Open new store
    xrefs = Store(filepath=dst,
                  keys=chunk_keys(
                      keys=sorted(mysql.structure.get_structures(my_uri)),
                      chunk_size=100
                  ),
                  tmpdir=tmpdir)

    protein_counts = {}
    cnt_proteins = 0
    for protein_acc, protein in proteins:
        try:
            pdbe_ids = protein2structures[protein_acc]
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

        obj = {
            "entries": entry_databases,
            "sets": entry_sets,
            "taxa": {protein["taxon"]}
        }

        try:
            upid = protein2proteome[protein_acc]
        except KeyError:
            pass
        else:
            obj["proteomes"] = {upid}

        try:
            ida, ida_id = protein2ida[protein_acc]
        except KeyError:
            pass
        else:
            obj["domain_architectures"] = {ida}

        for pdbe_id in pdbe_ids:
            xrefs.update(pdbe_id, obj)

            if pdbe_id in protein_counts:
                protein_counts[pdbe_id] += 1
            else:
                protein_counts[pdbe_id] = 1

        cnt_proteins += 1
        if not cnt_proteins % sync_frequency:
            xrefs.sync()

        if not cnt_proteins % 10000000:
            logger.info('{:>12,}'.format(cnt_proteins))

    logger.info('{:>12,}'.format(cnt_proteins))

    for store in (proteins, protein2proteome, protein2matches, protein2ida):
        store.close()

    for pdbe_id, cnt in protein_counts.items():
        xrefs.update(pdbe_id, {"proteins": cnt})

    xrefs.merge(processes=processes)
    xrefs.close()

    logger.info("complete")


def get_lineages(my_uri: str):
    lineages = {}
    for tax in mysql.taxonomy.get_taxa(my_uri, lineage=True).values():
        lineages[tax["taxId"]] = tax["lineage"]
    return lineages


def get_protein2structures(my_uri: str) -> Dict[str, Set[str]]:
    protein2pdb = {}
    for pdb_id, s in mysql.structure.get_structures(my_uri).items():
        for protein_ac in s["proteins"]:
            if protein_ac in protein2pdb:
                protein2pdb[protein_ac].add(pdb_id)
            else:
                protein2pdb[protein_ac] = {pdb_id}

    return protein2pdb


def export_goa_mappings(my_url: str, ora_url: str, outdir: str):
    databases = mysql.database.get_databases(my_url)
    interpro = databases["interpro"]
    version = interpro["version"]
    release_date = interpro["release_date"]

    con, cur = dbms.connect(ora_url)

    logger.info("exporting PDB-InterPro-GO-UniProt mapping")
    logger.debug("\tloading PDBe sequences from UniParc")
    cur.execute(
        """
        SELECT UPI, AC
        FROM UNIPARC.XREF
        WHERE DBID = 21
        AND DELETED = 'N'
        """
    )
    sequences = {}
    for upi, pdbe_acc in cur:
        if upi in sequences:
            sequences[upi]["structures"].add(pdbe_acc)
        else:
            sequences[upi] = {
                "structures": {pdbe_acc},
                "entries": set()
            }

    logger.debug("\tloading integrated signatures")
    # Only integrated signatures whose entry is checked and has GO terms
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

    logger.debug("\tloading GO terms in InterPro")
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

    logger.debug("\tloading PDBe matches")
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

    logger.debug("\tloading PDBe taxonomy")
    structures = pdbe.get_chain_taxonomy(cur)

    logger.debug("\tloading UniProt accessions")
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
    for pdbe_acc, protein_acc in cur:
        if not protein_acc:
            continue
        elif pdbe_acc in pdb2uniprot:
            pdb2uniprot[pdbe_acc].add(protein_acc)
        else:
            pdb2uniprot[pdbe_acc] = {protein_acc}

    os.makedirs(outdir, exist_ok=True)

    logger.debug("\twriting 'pdb2interpro2go.tsv'")
    dst = os.path.join(outdir, "pdb2interpro2go.tsv")
    with open(dst + ".tmp", "wt") as fh:
        fh.write("#PDBe ID\tchain\tTaxon ID\t"
                 "InterPro accession\tGO ID\tUniProt accession\n")

        for seq in sequences.values():
            for pdbe_acc in seq["structures"]:
                try:
                    s = structures[pdbe_acc]
                except KeyError:
                    # Structure does not exist in PDBe database
                    continue

                pdbe_id = s["id"]
                chain = s["chain"]
                proteins = pdb2uniprot.get(pdbe_acc, {''})

                for tax_id in s["taxa"]:
                    for entry_acc in seq["entries"]:
                        for go_id in entries[entry_acc]:
                            for protein_acc in proteins:
                                fh.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(
                                    pdbe_id, chain, tax_id, entry_acc,
                                    go_id, protein_acc
                                ))

    try:
        os.remove(dst)
    except FileNotFoundError:
        pass
    finally:
        os.rename(dst + ".tmp", dst)
        os.chmod(dst, 0o777)

    logger.info("exporting InterPro-GO-UniProt mapping")
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

    dst = os.path.join(outdir, "interpro2go2uniprot.tsv")
    with open(dst + ".tmp", "wt") as fh:
        fh.write("#InterPro accession\tGO ID\tUniProt accession\n")

        for row in cur:
            fh.write('\t'.join(row) + '\n')

    try:
        os.remove(dst)
    except FileNotFoundError:
        pass
    finally:
        os.rename(dst + ".tmp", dst)
        os.chmod(dst, 0o777)

    cur.close()
    con.close()

    dst = os.path.join(outdir, "release.txt")
    with open(dst, "wt") as fh:
        fh.write("InterPro version:        "
                 "{}\n".format(version))

        fh.write("Release date:            "
                 "{:%A, %d %B %Y}\n".format(release_date))

        fh.write("Generated on:            "
                 "{:%Y-%m-%d %H:%M}\n".format(datetime.now()))
    os.chmod(dst, 0o777)

    logger.info("complete")
