# -*- coding: utf-8 -*-

from typing import Tuple

import cx_Oracle

from interpro7dw.ebi.interpro.utils import repr_fragment


def get_structures(url: str):
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
                "citations": entry_citations.get(pdb_id, {}),
                "secondary_structures": entry_sec_structures.get(pdb_id, [])
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
                fragments.sort(key=repr_protein)


def repr_protein(fragment: dict) -> Tuple[int, int]:
    return fragment["protein_start"], fragment["protein_end"]
