#!/usr/bin/env python
# -*- coding: utf-8 -*-

from i7dw import dbms


def get_structures(uri, citations=True, fragments=False, by_protein=False):
    con, cur = dbms.connect(uri)

    # PDBe does not store CRC64 in hexa, hence we have to convert it for the join
    # cur.execute(
    #     """
    #     SELECT DISTINCT
    #       E.ID,
    #       E.TITLE,
    #       E.METHOD_CLASS,
    #       E.RESOLUTION,
    #       E.FIRST_REV_DATE,
    #       U.ACCESSION,
    #       U.AUTH_ASYM_ID,
    #       U.UNP_START,
    #       U.UNP_END
    #     FROM PDBE.ENTRY@PDBE_LIVE E
    #     INNER JOIN SIFTS_ADMIN.SIFTS_XREF_SEGMENT@PDBE_LIVE U ON (
    #       E.ID = U.ENTRY_ID AND
    #       E.METHOD_CLASS IN ('nmr', 'x-ray') AND
    #       U.UNP_START IS NOT NULL AND
    #       U.UNP_END IS NOT NULL AND
    #       U.PDB_START IS NOT NULL AND
    #       U.PDB_END IS NOT NULL AND
    #       U.UNP_END - U.UNP_START > 10
    #     )
    #     INNER JOIN SIFTS_ADMIN.SPTR_DBENTRY@PDBE_LIVE DB ON U.ACCESSION = DB.ACCESSION
    #     INNER JOIN SIFTS_ADMIN.SPTR_SEQUENCE@PDBE_LIVE S ON DB.DBENTRY_ID = S.DBENTRY_ID
    #     INNER JOIN INTERPRO.PROTEIN P ON (
    #       U.ACCESSION = P.PROTEIN_AC AND
    #       P.CRC64 = LPAD(TRIM(TO_CHAR(S.CHECKSUM, 'XXXXXXXXXXXXXXXX')),16,'0')
    #     )
    #     LEFT OUTER JOIN INTERPRO.PROTEIN_ACCPAIR PS ON U.ACCESSION = PS.SECONDARY_AC
    #     ORDER BY E.ID, U.AUTH_ASYM_ID, U.UNP_START, U.UNP_END
    #     """
    # )

    cur.execute(
        """
        SELECT E.ENTRY_ID, E.TITLE, E.METHOD, E.RESOLUTION, X.FIRST_REV_DATE, E.SPTR_AC, E.CHAIN, E.BEG_SEQ, E.END_SEQ
        FROM INTERPRO.UNIPROT_PDBE E
        INNER JOIN PDBE.ENTRY@PDBE_LIVE X ON E.ENTRY_ID = X.ID
        ORDER BY E.ENTRY_ID, E.CHAIN, E.BEG_SEQ, E.END_SEQ
        """
    )

    structures = {}
    for row in cur:
        pdb_ac = row[0]
        if pdb_ac in structures:
            s2 = structures[pdb_ac]
        else:
            s2 = structures[pdb_ac] = {
                'accession': pdb_ac,
                'date': row[4],
                'name': row[1],
                'resolution': row[3],
                'evidence': row[2],
                'proteins': {},
                'citations': {}
            }

        protein_ac = row[5]
        if protein_ac in s2['proteins']:
            p = s2['proteins'][protein_ac]
        else:
            p = s2['proteins'][protein_ac] = {} if fragments else set()

        chain = row[6]
        if fragments:
            if chain not in p:
                p[chain] = []
            p[chain].append({'start': row[7], 'end': row[8]})
        else:
            p.add(chain)

    if citations:
        # Get citations for PDBe structures
        cur.execute(
            """
            SELECT
              LOWER(E.ID),
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
            INNER JOIN CITATION@PDBE_LIVE C ON E.ID = C.ENTRY_ID
            INNER JOIN CITATION_AUTHOR@PDBE_LIVE A ON C.ENTRY_ID = A.ENTRY_ID AND C.ID = A.CITATION_ID
            WHERE E.METHOD_CLASS IN ('nmr', 'x-ray')
            ORDER BY E.ID, C.ID, A.ORDINAL
            """
        )

        for row in cur:
            pdb_ac = row[0]

            if pdb_ac not in structures:
                continue
            else:
                citations = structures[pdb_ac]['citations']

            pub_id = row[1]

            if pub_id not in citations:
                if row[5] is None:
                    pages = None
                elif row[6] is None:
                    pages = str(row[5])
                else:
                    pages = '{}-{}'.format(row[5], row[6])

                citations[pub_id] = {
                    'authors': [],
                    'DOI_URL': row[9],
                    'ISO_journal': row[3],
                    'raw_pages': pages,
                    'PMID': int(row[8]) if row[8] is not None else None,
                    'title': row[2],
                    'type': row[10],
                    'volume': row[4],
                    'year': row[7]
                }

            citations[pub_id]['authors'].append(row[11])

    if by_protein:
        proteins = {}

        for pdb_ac in structures:
            s = structures[pdb_ac]

            for protein_ac in s['proteins']:
                if protein_ac in proteins:
                    p = proteins[protein_ac]
                else:
                    p = proteins[protein_ac] = {}

                p[pdb_ac] = s

        return proteins
    else:
        return list(structures.values())
