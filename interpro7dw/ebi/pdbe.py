# -*- coding: utf-8 -*-

from typing import Tuple

import cx_Oracle

from interpro7dw.utils import datadump
from interpro7dw.ebi.interpro.utils import repr_fragment


def export_structures(url: str, output: str):
    con = cx_Oracle.connect(url)
    cur = con.cursor()

    # Retrieve citations
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
        FROM ENTRY@PDBE_LIVE E
        INNER JOIN CITATION@PDBE_LIVE C
          ON E.ID = C.ENTRY_ID
        INNER JOIN CITATION_AUTHOR@PDBE_LIVE A
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
    cur.execute(
        """
        SELECT SS.ENTRY_ID, SS.STRUCT_ASYM_ID, SS.ELEMENT_TYPE,
          R1.UNP_SEQ_ID AS POS_FROM, R1.CHEM_COMP_ID AS RES_FROM,
          R2.UNP_SEQ_ID AS POS_TO, R2.CHEM_COMP_ID AS RES_TO
        FROM (
          SELECT ENTRY_ID, STRUCT_ASYM_ID, ELEMENT_TYPE,
            RESIDUE_BEG_ID, RESIDUE_END_ID
          FROM PDBE.SS_HELIX@PDBE_LIVE
          UNION ALL
          SELECT ENTRY_ID, STRUCT_ASYM_ID, ELEMENT_TYPE,
            RESIDUE_BEG_ID, RESIDUE_END_ID
          FROM PDBE.SS_STRAND@PDBE_LIVE
        ) SS
        INNER JOIN SIFTS_ADMIN.SIFTS_XREF_RESIDUE@PDBE_LIVE R1
          ON (SS.ENTRY_ID=R1.ENTRY_ID
            AND SS.STRUCT_ASYM_ID=R1.STRUCT_ASYM_ID
            AND SS.RESIDUE_BEG_ID=R1.ID
            AND R1.CANONICAL_ACC=1
            AND R1.OBSERVED='Y'
            AND R1.UNP_SEQ_ID IS NOT NULL)
        INNER JOIN SIFTS_ADMIN.SIFTS_XREF_RESIDUE@PDBE_LIVE R2
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
    for pdb_id, dict_chains in entry_sec_structures.items():
        list_chains = []

        for chain_id in sorted(dict_chains):
            locations = []

            for elem_type, fragments in dict_chains[chain_id].items():
                fragments.sort(key=repr_fragment)
                locations.append({
                    "fragments": fragments
                })

            list_chains.append({
                "accession": chain_id,
                "locations": locations
            })

        entry_sec_structures[pdb_id] = list_chains

    """
    Retrieve PDBe entries with the proteins they are associated to
    (CRC64 not stored in hexadecimal, so need to convert)
    """
    cur.execute(
        """
        SELECT DISTINCT
          E.ID,
          E.TITLE,
          E.METHOD_CLASS,
          E.RESOLUTION,
          E.FIRST_REV_DATE,
          U.ACCESSION,
          U.AUTH_ASYM_ID,
          U.UNP_START,
          U.UNP_END,
          U.PDB_START,
          U.PDB_END,
          U.AUTH_START,
          U.AUTH_END
        FROM PDBE.ENTRY@PDBE_LIVE E
        INNER JOIN SIFTS_ADMIN.SIFTS_XREF_SEGMENT@PDBE_LIVE U ON (
          E.ID = U.ENTRY_ID AND
          E.METHOD_CLASS IN ('nmr', 'x-ray', 'em') AND
          U.UNP_START IS NOT NULL AND
          U.UNP_END IS NOT NULL AND
          U.PDB_START IS NOT NULL AND
          U.PDB_END IS NOT NULL
        )
        INNER JOIN SIFTS_ADMIN.SPTR_DBENTRY@PDBE_LIVE DB
          ON U.ACCESSION = DB.ACCESSION
        INNER JOIN SIFTS_ADMIN.SPTR_SEQUENCE@PDBE_LIVE S
          ON DB.DBENTRY_ID = S.DBENTRY_ID
        INNER JOIN INTERPRO.PROTEIN P ON (
          U.ACCESSION = P.PROTEIN_AC AND
          P.CRC64 = LPAD(TRIM(TO_CHAR(S.CHECKSUM, 'XXXXXXXXXXXXXXXX')),16,'0')
        )
        """
    )

    entries = {}
    for row in cur:
        pdb_id = row[0]
        try:
            entry = entries[pdb_id]
        except KeyError:
            entry = entries[pdb_id] = {
                "id": pdb_id,
                "date": row[4],
                "name": row[1],
                "resolution": row[3],
                "evidence": row[2],
                "proteins": {},
                "citations": entry_citations.get(pdb_id),
                "secondary_structures": entry_sec_structures.get(pdb_id)
            }

        protein_ac = row[5]
        try:
            chains = entry["proteins"][protein_ac]
        except KeyError:
            chains = entry["proteins"][protein_ac] = {}

        chain_id = row[6]
        try:
            chain = chains[chain_id]
        except KeyError:
            chain = chains[chain_id] = []

        unp_start = row[7]
        unp_end = row[8]
        if unp_start > unp_end:
            unp_start, unp_end = unp_end, unp_start

        chain.append({
            "protein_start": unp_start,
            "protein_end": unp_end,
            "structure_start": row[9],
            "structure_end": row[10],
            "author_structure_start": row[11],
            "author_structure_end": row[12]
        })

    cur.close()
    con.close()

    # Sort chains by fragment
    for entry in entries.values():
        for chains in entry["proteins"].values():
            for fragments in chains.values():
                fragments.sort(key=_repr_protein)

    datadump(output, entries)


def _repr_protein(fragment: dict) -> Tuple[int, int]:
    return fragment["protein_start"], fragment["protein_end"]


def get_chain_taxonomy(cur: cx_Oracle.Cursor) -> dict:
    cur.execute(
        """
        SELECT DISTINCT 
          ASYM.ENTRY_ID, 
          ASYM.AUTH_ASYM_ID, 
          SRC.TAX_ID
        FROM PDBE.STRUCT_ASYM@PDBE_LIVE ASYM
        INNER JOIN PDBE.ENTITY_SRC@PDBE_LIVE SRC
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

    return structures


def get_scop_domains(url: str) -> dict:
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    cur.execute(
        """
        SELECT
          ES.ENTRY_ID, SC.FAMILY_ID, ES.SCOP_SUPERFAMILY, ES.SCOP_FOLD,
          ES.SCOP_FAMILY, ES.SCOP_CLASS, SC.SCCS, ES.AUTH_ASYM_ID,
          ES.SCOP_ID, ES."START", ES.END, SC.BEG_SEQ, SC.END_SEQ
        FROM SIFTS_ADMIN.ENTITY_SCOP@PDBE_LIVE ES
        INNER JOIN SIFTS_ADMIN.SCOP_CLASS@PDBE_LIVE SC
          ON ES.SUNID = SC.SUNID
        WHERE ES."START" IS NOT NULL
        AND ES.END IS NOT NULL
        AND SC.BEG_SEQ IS NOT NULL
        AND SC.END_SEQ IS NOT NULL
        """
    )

    domains = {}
    for row in cur:
        pdb_id = row[0]
        family_id = row[1]
        superfamily_desc = row[2]
        fold_desc = row[3]
        family_desc = row[4]
        class_desc = row[5]
        sccs = row[6]  # SCOP Concise Classification String
        chain_id = row[7]
        scop_id = row[8]

        # PDB locations
        start = int(row[9])
        end = int(row[10])

        # UniProt locations
        seq_start = int(row[11])
        seq_end = int(row[12])

        if pdb_id in domains:
            s = domains[pdb_id]
        else:
            s = domains[pdb_id] = {}

        if family_id in s:
            fam = s[family_id]
        else:
            fam = s[family_id] = {
                "id": family_id,
                "family": family_desc,
                "superfamily": superfamily_desc,
                "fold": fold_desc,
                "class": class_desc,
                "sccs": sccs,
                "mappings": {}
            }

        fam["mappings"][scop_id] = {
            "chain": chain_id,
            "start": start,
            "end": end,
            "seq_start": seq_start,
            "seq_end": seq_end
        }

    cur.close()
    con.close()

    # Transform dictionary to fit format expected by InterPro7 API
    structures = {}
    for pdb_id, families in domains.items():
        structures[pdb_id] = {}
        for fam in families.values():
            for scop_id, mapping in fam["mappings"].items():
                structures[pdb_id][scop_id] = {
                    "class_id": scop_id,
                    "domain_id": fam["sccs"],
                    "coordinates": [{
                        "start": mapping["seq_start"],
                        "end": mapping["seq_end"]
                    }],

                    # Bonus (not used by API/client as of Nov 2018)
                    "description": fam["family"],
                    "chain": mapping["chain"]
                }

    return structures


def get_cath_domains(url: str) -> dict:
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    cur.execute(
        """
        SELECT
          EC.ENTRY_ID, EC.ACCESSION, CD.HOMOL, CD.TOPOL, CD.ARCH, CD.CLASS,
          -- CD.NAME,
          EC.DOMAIN, EC.AUTH_ASYM_ID, EC."START", EC.END,
          CS.BEG_SEQ, CS.END_SEQ
        FROM SIFTS_ADMIN.ENTITY_CATH@PDBE_LIVE EC
          INNER JOIN SIFTS_ADMIN.CATH_DOMAIN@PDBE_LIVE CD ON (
            EC.ENTRY_ID = CD.ENTRY
            AND EC.ACCESSION = CD.CATHCODE
            AND EC.DOMAIN = CD.DOMAIN
            AND EC.AUTH_ASYM_ID = CD.AUTH_ASYM_ID
          )
          INNER JOIN SIFTS_ADMIN.CATH_SEGMENT@PDBE_LIVE CS ON (
            EC.ENTRY_ID = CS.ENTRY
            AND EC.DOMAIN = CS.DOMAIN
            AND EC.AUTH_ASYM_ID = CS.AUTH_ASYM_ID
          )
          WHERE EC."START" IS NOT NULL
          AND EC.END IS NOT NULL
          AND CS.BEG_SEQ IS NOT NULL
          AND CS.END_SEQ IS NOT NULL
        """
    )

    domains = {}
    for row in cur:
        pdb_id = row[0]
        cath_id = row[1]
        homology = row[2]
        topology = row[3]
        architecture = row[4]
        _class = row[5]
        # name = row[6].read()  # CLOB
        domain_id = row[6]
        chain_id = row[7]

        # PDB locations
        start = int(row[8])
        end = int(row[9])

        # UniProt locations
        seq_start = int(row[10])
        seq_end = int(row[11])

        if pdb_id in domains:
            s = domains[pdb_id]
        else:
            s = domains[pdb_id] = {}

        if cath_id in s:
            fam = s[cath_id]
        else:
            fam = s[cath_id] = {
                "id": cath_id,
                "homology": homology,
                "topology": topology,
                "architecture": architecture,
                "class": _class,
                # "name": name,
                "mappings": {}
            }

        fam["mappings"][domain_id] = {
            "chain": chain_id,
            "start": start,
            "end": end,
            "seq_start": seq_start,
            "seq_end": seq_end
        }

    cur.close()
    con.close()

    # Transform dictionary to fit format expected by InterPro7 API
    structures = {}
    for pdb_id, families in domains.items():
        structures[pdb_id] = {}
        for fam in families.values():
            for domain_id, mapping in fam["mappings"].items():
                structures[pdb_id][domain_id] = {
                    "class_id": domain_id,
                    "domain_id": fam["id"],
                    "coordinates": [{
                        "start": mapping["seq_start"],
                        "end": mapping["seq_end"]
                    }],

                    # Bonus (not used by API/client as of Nov 2018)
                    "description": fam["topology"],
                    "chain": mapping["chain"]
                }

    return structures
