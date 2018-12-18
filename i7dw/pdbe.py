#!/usr/bin/env python
# -*- coding: utf-8 -*-

from . import dbms


def get_structures(uri: str) -> dict:
    con, cur = dbms.connect(uri)

    """
    Filters:
        - only nrm/x-ray
        - check for CRC64 mismatches (not stored in hexa so need to convert)
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
          E.METHOD_CLASS IN ('nmr', 'x-ray') AND
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

    structures = {}
    for row in cur:
        pdbe_id = row[0]
        if pdbe_id in structures:
            s = structures[pdbe_id]
        else:
            s = structures[pdbe_id] = {
                "id": pdbe_id,
                "date": row[4],
                "name": row[1],
                "resolution": row[3],
                "evidence": row[2],
                "proteins": {},
                "citations": {}
            }

        protein_ac = row[5]
        if protein_ac in s["proteins"]:
            chains = s["proteins"][protein_ac]
        else:
            chains = s["proteins"][protein_ac] = {}

        chain_id = row[6]
        if chain_id in chains:
            chain = chains[chain_id]
        else:
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
        WHERE E.METHOD_CLASS IN ('nmr', 'x-ray')
        ORDER BY E.ID, C.ID, A.ORDINAL
        """
    )

    for row in cur:
        pdbe_id = row[0]

        if pdbe_id not in structures:
            continue
        else:
            citations = structures[pdbe_id]["citations"]

        pub_id = row[1]
        if pub_id not in citations:
            if row[5] is None:
                pages = None
            elif row[6] is None:
                pages = str(row[5])
            else:
                pages = "{}-{}".format(row[5], row[6])

            citations[pub_id] = {
                "authors": [],
                "DOI_URL": row[9],
                "ISO_journal": row[3],
                "raw_pages": pages,
                "PMID": int(row[8]) if row[8] is not None else None,
                "title": row[2],
                "type": row[10],
                "volume": row[4],
                "year": row[7]
            }

        citations[pub_id]["authors"].append(row[11])

    cur.close()
    con.close()

    for s in structures.values():
        for chains in s["proteins"].values():
            for chain_id in chains:
                chains[chain_id].sort(key=lambda x: (
                    x["protein_start"],
                    x["protein_end"]
                ))

    return structures


def get_scop_domains(uri: str) -> dict:
    con, cur = dbms.connect(uri)
    cur.execute(
        """
        SELECT
          ES.ENTRY_ID, SC.FAMILY_ID, ES.SCOP_SUPERFAMILY, ES.SCOP_FOLD,
          ES.SCOP_FAMILY, ES.SCOP_CLASS, SC.SCCS, ES.AUTH_ASYM_ID,
          ES.SCOP_ID, ES."START", ES.END, SC.BEG_SEQ, SC.END_SEQ
        FROM SIFTS_ADMIN.ENTITY_SCOP ES
        INNER JOIN SIFTS_ADMIN.SCOP_CLASS SC
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


def get_cath_domains(uri: str) -> dict:
    con, cur = dbms.connect(uri)
    cur.execute(
        """
        SELECT
          EC.ENTRY_ID, EC.ACCESSION, CD.HOMOL, CD.TOPOL, CD.ARCH, CD.CLASS,
          CD.NAME, EC.DOMAIN, EC.AUTH_ASYM_ID, EC."START", EC.END,
          CS.BEG_SEQ, CS.END_SEQ
        FROM SIFTS_ADMIN.ENTITY_CATH EC
          INNER JOIN SIFTS_ADMIN.CATH_DOMAIN CD ON (
            EC.ENTRY_ID = CD.ENTRY
            AND EC.ACCESSION = CD.CATHCODE
            AND EC.DOMAIN = CD.DOMAIN
            AND EC.AUTH_ASYM_ID = CD.AUTH_ASYM_ID
          )
          INNER JOIN SIFTS_ADMIN.CATH_SEGMENT CS ON (
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
        name = row[6].read()  # CLOB
        domain_id = row[7]
        chain_id = row[8]

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

        if cath_id in s:
            fam = s[cath_id]
        else:
            fam = s[cath_id] = {
                "id": cath_id,
                "homology": homology,
                "topology": topology,
                "architecture": architecture,
                "class": _class,
                "name": name,
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
