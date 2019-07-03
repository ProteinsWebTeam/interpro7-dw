import os
import time
import shutil
from datetime import datetime
from multiprocessing import Process, Queue
from typing import Dict, Optional, Set

from . import mysql
from .. import dbms, logger, pdbe
from ..io import KVdb, Store, traverse


def chunk_keys(keys: list, chunk_size: int) -> list:
    return [keys[i] for i in range(0, len(keys), chunk_size)]


def export_entries(my_uri: str, src_proteins: str, src_proteomes:str,
                   src_matches: str, src_ida: str, dst: str,
                   processes: int=1, tmpdir: Optional[str]=None,
                   sync_frequency: int=100000):
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
    xrefs = Store(filepath=dst,
                  keys=chunk_keys(
                      keys=sorted(entries),
                      chunk_size=100
                  ),
                  tmpdir=tmpdir)
    taxon_protein_counts = {}
    entry_match_counts = {}
    cnt_proteins = 0
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
            xrefs.update(entry_acc, obj)
            if entry_acc in entry_match_counts:
                entry_match_counts[entry_acc] += cnt
            else:
                entry_match_counts[entry_acc] = cnt

        cnt_proteins += 1
        if not cnt_proteins % sync_frequency:
            xrefs.sync()

    for store in (proteins, protein2proteome, protein2matches, protein2ida):
        store.close()

    sets = mysql.entry.get_sets(my_uri)
    entry2set = {}
    for set_acc, s in sets.items():
        for entry_acc in s["members"]:
            entry2set[entry_acc] = set_acc

    for entry_acc, set_acc in entry2set.items():
        xrefs.update(entry_acc, {"sets": {set_acc}})

    for entry_acc, cnt in entry_match_counts.items():
        xrefs.update(entry_acc, {"matches": cnt})

    xrefs.merge(processes=processes)
    xrefs.close()


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

    keys = chunk_keys(keys=sorted(lineages), chunk_size=10)
    with Store(keys=keys, tmpdir=tmpdir) as xrefs:
        protein_counts = {}
        cnt_proteins = 0
        for protein_acc, protein in proteins:
            protein_tax_id = protein["taxon"]
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

            xrefs.update(protein_tax_id, obj)
            for tax_id in lineages[protein_tax_id]:
                if tax_id in protein_counts:
                    protein_counts[tax_id] += 1
                else:
                    protein_counts[tax_id] = 1

            cnt_proteins += 1
            if not cnt_proteins % sync_frequency:
                xrefs.sync()
                logger.debug(f"{cnt_proteins:>12.}")

        for store in (proteins, protein2proteome, protein2matches, protein2ida):
            store.close()

        size = xrefs.merge(processes=processes)
        logger.debug(f"{cnt_proteins:>12.}")
        with KVdb(dir=tmpdir, writeback=True) as kvdb:
            for tax_id, obj in xrefs:
                # Propagate to lineage
                for node_id in lineages[tax_id]:
                    try:
                        node = kvdb[node_id]
                    except KeyError:
                        node = obj
                    else:
                        traverse(obj, node)
                    finally:
                        node["proteins"] = protein_counts[node_id]
                        kvdb[node_id] = node

                kvdb.sync()
                logger.debug(f"{tax_id:>12.}")
            size += kvdb.size
            shutil.copy(kvdb.filepath, dst)

    logger.info("Disk usage: {:.0f}MB".format(size/1024**2))


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

    for store in (proteins, protein2proteome, protein2matches, protein2ida):
        store.close()

    for upid, cnt in protein_counts.items():
        xrefs.update(upid, {"proteins": cnt})

    size = xrefs.merge(processes=processes)
    xrefs.close()

    logger.info("Disk usage: {:.0f}MB".format(size/1024**2))


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

    size = xrefs.merge(processes=processes)
    xrefs.close()

    logger.info("Disk usage: {:.0f}MB".format(size / 1024 ** 2))


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
