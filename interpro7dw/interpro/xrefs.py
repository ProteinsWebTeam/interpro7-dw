import re

import cx_Oracle

from interpro7dw import metacyc, uniprot
from interpro7dw.interpro.utils import overlaps_pdb_chain
from interpro7dw.utils import logger
from interpro7dw.utils.store import SimpleStore, SimpleStoreSorter, Store
from interpro7dw.utils.store import copy_dict, loadobj


def dump_entries(ipr_url: str, unp_url: str, proteins_file: str,
                 matches_file: str, proteomes_file: str, domorgs_file: str,
                 structures_file: str, metacyc_file: str, alphafold_file: str,
                 xrefs_file: str, **kwargs):
    """Export InterPro entries and member database signatures cross-references.
    For each entry or signature, the following information is saved:
        - proteins matched (and number of matches)
        - proteomes
        - PDBe structures
        - taxa
        - domain organisations
        - ENZYME numbers
        - MetaCyc and Reactome pathways

    :param ipr_url: InterPro Oracle connection string.
    :param unp_url: UniProt Oracle connection string.
    :param proteins_file: Store file of protein info.
    :param matches_file: Store file of protein matches.
    :param proteomes_file: Store file of protein-proteome mapping.
    :param domorgs_file: Store file of domain organisations.
    :param structures_file: File of PDBe structures.
    :param metacyc_file: MetaCyc tar archive.
    :param alphafold_file: CSV file listing AlphaFold predictions.
    :param xrefs_file: Output SimpleStore file.
    """
    buffersize = kwargs.get("buffersize", 1000000)
    tempdir = kwargs.get("tempdir")

    logger.info("loading Swiss-Prot data")
    protein2enzymes = uniprot.misc.get_swissprot2enzyme(unp_url)
    protein2reactome = uniprot.misc.get_swissprot2reactome(unp_url)

    # Create mapping protein -> structure -> chain -> locations
    logger.info("loading PDBe structures")
    protein2structures = {}
    for pdbe_id, entry in loadobj(structures_file).items():
        for protein_acc, chains in entry["proteins"].items():
            try:
                protein2structures[protein_acc][pdbe_id] = chains
            except KeyError:
                protein2structures[protein_acc] = {pdbe_id: chains}

    logger.info("iterating proteins")
    with SimpleStoreSorter(tempdir=tempdir) as stores:
        proteins = Store(proteins_file, "r")
        matches = Store(matches_file, "r")
        proteomes = Store(proteomes_file, "r")
        domorgs = Store(domorgs_file, "r")

        xrefs = {}
        num_xrefs = 0
        i = 0
        for i, (protein_acc, protein_matches) in enumerate(matches.items()):
            protein = proteins[protein_acc]
            protein_id = protein["identifier"]
            taxon_id = protein["taxid"]
            proteome_id = proteomes.get(protein_acc)
            structures = protein2structures.get(protein_acc, {})
            try:
                _, dom_id, dom_members, _ = domorgs[protein_acc]
            except KeyError:
                dom_id = None
                dom_members = []

            for entry_acc, locations in protein_matches.items():
                try:
                    entry_xrefs = xrefs[entry_acc]
                except KeyError:
                    entry_xrefs = xrefs[entry_acc] = {
                        "dom_orgs": set(),
                        "enzymes": set(),
                        "matches": 0,
                        "reactome": set(),
                        "proteins": set(),
                        "proteomes": set(),
                        "structures": set(),
                        "taxa": set(),
                    }

                entry_xrefs["matches"] += len(locations)
                entry_xrefs["proteins"].add((protein_acc, protein_id))
                entry_xrefs["taxa"].add(taxon_id)
                num_xrefs += 4

                if entry_acc in dom_members:
                    entry_xrefs["dom_orgs"].add(dom_id)
                    num_xrefs += 1

                if proteome_id:
                    entry_xrefs["proteomes"].add(proteome_id)
                    num_xrefs += 1

                for pdbe_id, chains in structures.items():
                    for chain_id, segments in chains.items():
                        if overlaps_pdb_chain(locations, segments):
                            entry_xrefs["structures"].add(pdbe_id)
                            num_xrefs += 1
                            break  # Skip other chains

                if re.fullmatch(r"IPR\d{6}", entry_acc):
                    for ecno in protein2enzymes.get(protein_acc, []):
                        entry_xrefs["enzymes"].add(ecno)
                        num_xrefs += 1

                    pathways = protein2reactome.get(protein_acc, [])
                    for pathway_id, pathway_name in pathways:
                        entry_xrefs["reactome"].add((pathway_id, pathway_name))
                        num_xrefs += 2

            if num_xrefs >= buffersize:
                stores.dump(xrefs)
                xrefs.clear()
                num_xrefs = 0

            if (i + 1) % 10000000 == 0:
                logger.info(f"{i + 1:>15,}")

        logger.info(f"{i + 1:>15,}")

        proteins.close()
        matches.close()
        proteomes.close()
        domorgs.close()

        stores.dump(xrefs)
        xrefs.clear()

        logger.info(f"temporary files: {stores.size / 1024 / 1024:.0f} MB")

        logger.info("loading MetaCyc pathways")
        ec2metacyc = metacyc.get_ec2pathways(metacyc_file)

        logger.info("loading structural models")
        con = cx_Oracle.connect(ipr_url)
        cur = con.cursor()
        cur.execute(
            """
            SELECT LOWER(ALGORITHM), METHOD_AC, COUNT(*)
            FROM INTERPRO.STRUCT_MODEL
            GROUP BY METHOD_AC, ALGORITHM
            """
        )
        struct_models = {}
        for algorithm, entry_acc, cnt in cur:
            try:
                struct_models[algorithm][entry_acc] = cnt
            except KeyError:
                struct_models[algorithm] = {entry_acc: cnt}

        cur.close()

        protein2alphafold = {}
        with open(alphafold_file, "rt") as fh:
            for line in fh:
                protein_acc, start, end, model_id = line.rstrip().split(',')
                try:
                    protein2alphafold[protein_acc] += 1
                except KeyError:
                    protein2alphafold[protein_acc] = 1

        logger.info("writing final file")
        with SimpleStore(xrefs_file) as store:
            for entry_acc, values in stores.merge():
                xrefs = {}

                for entry_xrefs in values:
                    copy_dict(entry_xrefs, xrefs, concat_or_incr=True)

                # Adds MetaCyc pathways
                pathways = set()
                for ecno in xrefs["enzymes"]:
                    for pathway_id, pathway_name in ec2metacyc.get(ecno, []):
                        pathways.add((pathway_id, pathway_name))

                xrefs["metacyc"] = pathways

                # Adds structural models
                num_struct_models = {}
                for algorithm, counts in struct_models.items():
                    num_struct_models[algorithm] = counts.get(entry_acc, 0)

                num_struct_models["alphafold"] = 0
                if re.fullmatch(r"IPR\d{6}", entry_acc):
                    for protein_acc, _ in xrefs["proteins"]:
                        try:
                            cnt = protein2alphafold[protein_acc]
                        except KeyError:
                            continue
                        else:
                            if cnt == 1:
                                num_struct_models["alphafold"] += 1

                xrefs["struct_models"] = num_struct_models
                store.add((entry_acc, xrefs))


def dump_proteomes(proteins_file: str, matches_file: str, proteomes_file: str,
                   domorgs_file: str, structures_file: str, entries_file: str,
                   xrefs_file: str, **kwargs):
    """Export proteome cross-references, that is:
        - protein matched
        - InterPro entries and member database signatures
        - clans
        - PDBe structures
        - taxa
        - domain organisations

    :param proteins_file:
    :param matches_file:
    :param proteomes_file:
    :param domorgs_file:
    :param structures_file:
    :param entries_file:
    :param xrefs_file:
    :param kwargs:
    :return:
    """
    buffersize = kwargs.get("buffersize", 1000000)
    tempdir = kwargs.get("tempdir")

    logger.info("loading entries")
    entries = loadobj(entries_file)

    # Create mapping protein -> structure -> chain -> locations
    logger.info("loading PDBe structures")
    protein2structures = {}
    for pdbe_id, entry in loadobj(structures_file).items():
        for protein_acc in entry["proteins"]:
            try:
                protein2structures[protein_acc].add(pdbe_id)
            except KeyError:
                protein2structures[protein_acc] = {pdbe_id}

    logger.info("iterating proteins")
    with SimpleStoreSorter(tempdir=tempdir) as stores:
        proteins = Store(proteins_file, "r")
        matches = Store(matches_file, "r")
        proteomes = Store(proteomes_file, "r")
        domorgs = Store(domorgs_file, "r")

        xrefs = {}
        num_xrefs = 0
        i = 0
        for i, (protein_acc, proteome_id) in enumerate(proteomes.items()):
            protein = proteins[protein_acc]
            taxon_id = protein["taxid"]

            try:
                proteome_xrefs = xrefs[proteome_id]
            except KeyError:
                proteome_xrefs = xrefs[proteome_id] = {
                    "domain_architectures": set(),
                    "entries": {},
                    "proteins": 0,
                    "sets": set(),
                    "structures": set(),
                    "taxa": set()
                }

            proteome_xrefs["proteins"] += 1

            try:
                _, dom_id, _, _ = domorgs[protein_acc]
            except KeyError:
                pass
            else:
                proteome_xrefs["domain_architectures"].add(dom_id)
                num_xrefs += 1

            for entry_acc in matches.get(protein_acc, []):
                entry = entries[entry_acc]
                try:
                    proteome_xrefs["entries"][entry.database].add(entry_acc)
                except KeyError:
                    proteome_xrefs["entries"][entry.database] = {entry_acc}

                num_xrefs += 1
                if entry.clan:
                    proteome_xrefs["sets"].add(entry.clan["accession"])
                    num_xrefs += 1

            try:
                pdbe_ids = protein2structures[protein_acc]
            except KeyError:
                pass
            else:
                proteome_xrefs["structures"] |= pdbe_ids
                num_xrefs += len(pdbe_ids)

            proteome_xrefs["taxa"].add(taxon_id)
            num_xrefs += 1

            if num_xrefs >= buffersize:
                stores.dump(xrefs)
                xrefs.clear()
                num_xrefs = 0

            if (i + 1) % 10000000 == 0:
                logger.info(f"{i + 1:>15,}")

        logger.info(f"{i + 1:>15,}")

        proteins.close()
        matches.close()
        proteomes.close()
        domorgs.close()

        stores.dump(xrefs)
        xrefs.clear()

        logger.info(f"temporary files: {stores.size / 1024 / 1024:.0f} MB")

        logger.info("writing final file")
        with SimpleStore(xrefs_file) as store:
            for proteome_id, values in stores.merge():
                xrefs = {}

                for entry_xrefs in values:
                    copy_dict(entry_xrefs, xrefs, concat_or_incr=True)

                store.add((proteome_id, xrefs))

    logger.info("done")
