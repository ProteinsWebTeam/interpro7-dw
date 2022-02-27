import os
import pickle

import cx_Oracle

from interpro7dw.utils import logger


_PDB2INTERPRO = "pdb2interpro.csv"


def export_segments(uri: str, output: str):
    logger.info("starting")
    os.makedirs(os.path.dirname(output), exist_ok=True)

    con = cx_Oracle.connect(uri)
    cur = con.cursor()
    cur.execute(
        """
        SELECT DISTINCT
          U.ACCESSION,
          E.ID,
          U.AUTH_ASYM_ID,
          U.UNP_START,
          U.UNP_END,
          U.PDB_START,
          U.PDB_END,
          U.AUTH_START,
          U.AUTH_END
          --LPAD(TRIM(TO_CHAR(S.CHECKSUM, 'XXXXXXXXXXXXXXXX')),16,'0') CRC64
        FROM PDBE.ENTRY E
        INNER JOIN SIFTS_ADMIN.SIFTS_XREF_SEGMENT U ON (
          E.ID = U.ENTRY_ID AND
          E.METHOD_CLASS IN ('nmr', 'x-ray', 'em') AND
          U.UNP_START IS NOT NULL AND
          U.UNP_END IS NOT NULL AND
          U.PDB_START IS NOT NULL AND
          U.PDB_END IS NOT NULL
        )
        INNER JOIN SIFTS_ADMIN.SPTR_DBENTRY DB
          ON U.ACCESSION = DB.ACCESSION
        --INNER JOIN SIFTS_ADMIN.SPTR_SEQUENCE S
        --  ON DB.DBENTRY_ID = S.DBENTRY_ID
        """
    )

    proteins = {}
    for rec in cur:
        protein_acc = rec[0]
        pdb_id = rec[1]
        chain_id = rec[2]
        protein_start = rec[3]
        protein_end = rec[4]
        structure_start = rec[5]
        structure_end = rec[6]
        author_structure_start = rec[7]
        author_structure_end = rec[8]

        if protein_start > protein_end:
            protein_start, protein_end = protein_end, protein_start

        try:
            structures = proteins[protein_acc]
        except KeyError:
            structures = proteins[protein_acc] = {}

        try:
            chains = structures[pdb_id]
        except KeyError:
            chains = structures[pdb_id] = {}

        segment = {
            "protein_start": protein_start,
            "protein_end": protein_end,
            "structure_start": structure_start,
            "structure_end": structure_end,
            "author_structure_start": author_structure_start,
            "author_structure_end": author_structure_end
        }

        try:
            chains[chain_id].append(segment)
        except KeyError:
            chains[chain_id] = [segment]

    cur.close()
    con.close()

    with open(output, "wb") as fh:
        pickle.dump(proteins, fh)

    logger.info("done")


def export_structures(uri: str, output: str):
    con = cx_Oracle.connect(uri)
    cur = con.cursor()

    # Retrieve citations
    logger.info("loading PDBe citations")
    cur.execute(
        """
        SELECT
          E.ID,
          C.ID,
          C.TITLE,
          C.JOURNAL_ABBREV,
          C.JOURNAL_VOLUME,
          C.PAGE_FIRST,
          C.PAGE_LAST,
          C.YEAR,
          C.DATABASE_ID_PUBMED,
          C.DATABASE_ID_DOI,
          C.CITATION_TYPE,
          A.NAME
        FROM PDBE.ENTRY E
        INNER JOIN PDBE.CITATION C
          ON E.ID = C.ENTRY_ID
        INNER JOIN PDBE.CITATION_AUTHOR A
          ON C.ENTRY_ID = A.ENTRY_ID AND C.ID = A.CITATION_ID
        ORDER BY E.ID, C.ID, A.ORDINAL
        """
    )

    entry_citations = {}
    for row in cur:
        pdb_id = row[0]

        try:
            entry = entry_citations[pdb_id]
        except KeyError:
            entry = entry_citations[pdb_id] = {}

        pub_id = row[1]
        try:
            pub = entry[pub_id]
        except KeyError:
            if row[5] and row[6]:
                pages = f"{row[5]}-{row[6]}"
            elif row[5]:
                pages = str(row[5])
            elif row[6]:
                pages = str(row[6])
            else:
                pages = None

            pub = entry[pub_id] = {
                "PMID": int(row[8]) if row[8] is not None else None,
                "volume": row[4],
                "year": row[7],
                "title": row[2],
                "raw_pages": pages,
                "ISO_journal": row[3],
                "authors": [],
                "DOI_URL": row[9],
                "type": row[10]
            }

        pub["authors"].append(row[11])

    # Retrieve secondary structures
    logger.info("loading secondary structures")
    cur.execute(
        """
        SELECT SS.ENTRY_ID, SS.STRUCT_ASYM_ID, SS.ELEMENT_TYPE,
          R1.UNP_SEQ_ID AS POS_FROM, R1.CHEM_COMP_ID AS RES_FROM,
          R2.UNP_SEQ_ID AS POS_TO, R2.CHEM_COMP_ID AS RES_TO
        FROM (
          SELECT ENTRY_ID, STRUCT_ASYM_ID, ELEMENT_TYPE,
            RESIDUE_BEG_ID, RESIDUE_END_ID
          FROM PDBE.SS_HELIX
          UNION ALL
          SELECT ENTRY_ID, STRUCT_ASYM_ID, ELEMENT_TYPE,
            RESIDUE_BEG_ID, RESIDUE_END_ID
          FROM PDBE.SS_STRAND
        ) SS
        INNER JOIN SIFTS_ADMIN.SIFTS_XREF_RESIDUE R1
          ON (SS.ENTRY_ID=R1.ENTRY_ID
            AND SS.STRUCT_ASYM_ID=R1.STRUCT_ASYM_ID
            AND SS.RESIDUE_BEG_ID=R1.ID
            AND R1.CANONICAL_ACC=1
            AND R1.OBSERVED='Y'
            AND R1.UNP_SEQ_ID IS NOT NULL)
        INNER JOIN SIFTS_ADMIN.SIFTS_XREF_RESIDUE R2
          ON (SS.ENTRY_ID=R2.ENTRY_ID
            AND SS.STRUCT_ASYM_ID=R2.STRUCT_ASYM_ID
            AND SS.RESIDUE_END_ID=R2.ID
            AND R2.CANONICAL_ACC=1
            AND R2.OBSERVED='Y'
            AND R2.UNP_SEQ_ID IS NOT NULL)
        """
    )

    entry_sec_structures = {}
    for row in cur:
        pdb_id = row[0]
        try:
            chains = entry_sec_structures[pdb_id]
        except KeyError:
            chains = entry_sec_structures[pdb_id] = {}

        chain_id = row[1]
        try:
            chain = chains[chain_id]
        except KeyError:
            chain = chains[chain_id] = {}

        elem_type = row[2]
        try:
            fragments = chain[elem_type]
        except KeyError:
            fragments = chain[elem_type] = []

        fragments.append({
            # add the type of secondary structure to the fragment
            "shape": elem_type,
            "start": row[3],
            "end": row[5],
            # "res_start": row[4],
            # "res_end": row[6],
        })

    # Sort chains by fragment
    for pdb_id, chains in entry_sec_structures.items():
        sorted_chains = []

        for chain_id in sorted(chains):
            locations = []

            for elem_type, fragments in chains[chain_id].items():
                fragments.sort(key=lambda f: (f["start"], f["end"]))
                locations.append({
                    "fragments": fragments
                })

            sorted_chains.append({
                "accession": chain_id,
                "locations": locations
            })

        entry_sec_structures[pdb_id] = sorted_chains

    logger.info("loading PDBe entries")
    cur.execute(
        """
        SELECT ID, TITLE, METHOD_CLASS, RESOLUTION, FIRST_REV_DATE
        FROM PDBE.ENTRY
        """
    )

    entries = {}
    for row in cur:
        pdb_id = row[0]
        entries[pdb_id] = {
            "id": pdb_id,
            "date": row[4],
            "name": row[1],
            "resolution": row[3],
            "evidence": row[2],
            "citations": entry_citations.get(pdb_id),
            "secondary_structures": entry_sec_structures.get(pdb_id)
        }

    cur.close()
    con.close()

    with open(output, "wb") as fh:
        pickle.dump({
            "entries": entries,
            "cath": get_cath_domains(uri),
            "scop": get_cath_domains(uri),
            "taxonomy": get_chain_taxonomy(uri)
        }, fh)

    logger.info("done")


def get_cath_domains(uri: str) -> dict:
    con = cx_Oracle.connect(uri)
    cur = con.cursor()
    cur.execute(
        """
        SELECT
          XS.ACCESSION, XS.ENTRY_ID, XS.AUTH_ASYM_ID, EC.DOMAIN,
          -- Superfamily info
          EC.ACCESSION, CD.HOMOL,
          -- UniProt range
          XS.UNP_START, XS.UNP_END
          ---- PDB range
          -- CS.BEG_SEQ, CS.END_SEQ,
          ---- PDBe range
          -- EC."START", EC.END
        FROM SIFTS_ADMIN.SIFTS_XREF_SEGMENT XS
          INNER JOIN SIFTS_ADMIN.ENTITY_CATH EC
            ON (XS.ENTRY_ID = EC.ENTRY_ID
                AND XS.ENTITY_ID = EC.ENTITY_ID
                AND XS.AUTH_ASYM_ID = EC.AUTH_ASYM_ID)
          INNER JOIN SIFTS_ADMIN.CATH_DOMAIN CD
            ON (EC.ENTRY_ID = CD.ENTRY
                AND EC.DOMAIN = CD.DOMAIN
                AND EC.AUTH_ASYM_ID = CD.AUTH_ASYM_ID)
          -- INNER JOIN SIFTS_ADMIN.CATH_SEGMENT CS
          --   ON (EC.ENTRY_ID = CS.ENTRY
          --       AND EC.DOMAIN = CS.DOMAIN
          --       AND EC.AUTH_ASYM_ID = CS.AUTH_ASYM_ID)
        WHERE XS.UNP_START IS NOT NULL
              AND XS.UNP_END IS NOT NULL
        """
    )

    domains = {}
    for row in cur:
        uniprot_acc = row[0]
        pdb_id = row[1]
        chain = row[2]
        domain_id = row[3]
        superfamily_id = row[4]
        superfamily_description = row[5]
        unp_start = int(row[6])
        unp_end = int(row[7])

        if unp_start > unp_end:
            unp_start, unp_end = unp_end, unp_start

        try:
            protein = domains[uniprot_acc]
        except KeyError:
            protein = domains[uniprot_acc] = {}

        try:
            dom = protein[domain_id]
        except KeyError:
            dom = protein[domain_id] = {
                "id": domain_id,
                "pdb_id": pdb_id,
                "chain": chain,
                "superfamily": {
                    "id": superfamily_id,
                    "description": superfamily_description
                },
                "locations": []
            }
        finally:
            dom["locations"].append({
                "start": unp_start,
                "end": unp_end
            })

    cur.close()
    con.close()

    for protein_domains in domains.values():
        for domain in protein_domains.values():
            domain["locations"].sort(key=lambda x: (x["start"], x["end"]))

    return domains


def get_scop_domains(uri: str) -> dict:
    con = cx_Oracle.connect(uri)
    cur = con.cursor()
    cur.execute(
        """
        SELECT
          XS.ACCESSION, XS.ENTRY_ID, XS.AUTH_ASYM_ID, ES.SCOP_ID,
          SD.SCCS, ES.SCOP_SUPERFAMILY,
          SD.BEG_SEQ, SD.END_SEQ
        FROM SIFTS_ADMIN.SIFTS_XREF_SEGMENT XS
          INNER JOIN SIFTS_ADMIN.ENTITY_SCOP ES
            ON (XS.ENTRY_ID = ES.ENTRY_ID
                AND XS.ENTITY_ID = ES.ENTITY_ID)
          INNER JOIN SIFTS_ADMIN.SCOP_DOMAIN SD
            ON (ES.ENTRY_ID = SD.ENTRY
                AND ES.AUTH_ASYM_ID = SD.AUTH_ASYM_ID
                AND ES.SCOP_ID = SD.SCOP_ID)
        WHERE SD.BEG_SEQ IS NOT NULL
              AND SD.END_SEQ IS NOT NULL
        """
    )

    domains = {}
    for row in cur:
        uniprot_acc = row[0]
        pdb_id = row[1]
        chain = row[2]
        domain_id = row[3]
        sccs = row[4]
        superfamily_description = row[5]
        unp_start = int(row[6])
        unp_end = int(row[7])

        if unp_start > unp_end:
            unp_start, unp_end = unp_end, unp_start

        try:
            protein = domains[uniprot_acc]
        except KeyError:
            protein = domains[uniprot_acc] = {}

        try:
            dom = protein[domain_id]
        except KeyError:
            dom = protein[domain_id] = {
                "id": domain_id,
                "pdb_id": pdb_id,
                "chain": chain,
                "superfamily": {
                    "id": sccs,
                    "description": superfamily_description
                },
                "locations": []
            }
        finally:
            dom["locations"].append({
                "start": unp_start,
                "end": unp_end
            })

    cur.close()
    con.close()

    for protein_domains in domains.values():
        for domain in protein_domains.values():
            domain["locations"].sort(key=lambda x: (x["start"], x["end"]))

    return domains


def get_chain_taxonomy(uri: str) -> dict:
    con = cx_Oracle.connect(uri)
    cur = con.cursor()
    cur.execute(
        """
        SELECT DISTINCT 
          ASYM.ENTRY_ID, 
          ASYM.AUTH_ASYM_ID, 
          SRC.TAX_ID
        FROM PDBE.STRUCT_ASYM ASYM
        INNER JOIN PDBE.ENTITY_SRC SRC
          ON ASYM.ENTRY_ID = SRC.ENTRY_ID AND ASYM.ENTITY_ID = SRC.ENTITY_ID
        """
    )

    structures = {}
    for pdb_id, chain, tax_id in cur:
        pdb_acc = pdb_id + '_' + chain

        if pdb_acc in structures:
            s = structures[pdb_acc]
        else:
            s = structures[pdb_acc] = {
                "id": pdb_id,
                "chain": chain,
                "taxa": set()
            }
        s["taxa"].add(tax_id)

    cur.close()
    con.close()

    return structures
