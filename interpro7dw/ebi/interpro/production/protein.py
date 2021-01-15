# -*- coding: utf-8 -*-

from typing import List, Optional, Sequence

import cx_Oracle

from interpro7dw import logger
from interpro7dw.utils import DirectoryTree, DumpFile, Store
from interpro7dw.ebi.interpro.utils import DC_STATUSES
from interpro7dw.ebi.interpro.utils import condense_locations
from interpro7dw.ebi.interpro.utils import repr_fragment


def chunk_proteins(url: str, keyfile: str, chunk_size: int = 50000):
    logger.info("loading")
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    cur.execute(
        """
        SELECT PROTEIN_AC
        FROM INTERPRO.PROTEIN
        """
    )

    accessions = [acc for acc, in cur]
    cur.close()
    con.close()

    logger.info("splitting into chunks")
    Store.dump_keys(Store.chunk(accessions, chunk_size), keyfile)
    logger.info("complete")


def export_features(url: str, keyfile: str, output: str,
                    processes: int = 1, tmpdir: Optional[str] = None):
    logger.info("starting")
    with Store(output, Store.load_keys(keyfile), tmpdir) as store:
        con = cx_Oracle.connect(url)
        cur = con.cursor()
        cur.execute(
            """
            SELECT FM.PROTEIN_AC, FM.METHOD_AC, LOWER(DB.DBSHORT),
                   FM.POS_FROM, FM.POS_TO, FM.SEQ_FEATURE
            FROM INTERPRO.FEATURE_MATCH FM
            INNER JOIN INTERPRO.CV_DATABASE DB ON FM.DBCODE = DB.DBCODE
            """
        )

        i = 0
        for row in cur:
            protein_acc = row[0]
            signature_acc = row[1]
            database = row[2]
            pos_start = row[3]
            pos_end = row[4]
            seq_feature = row[5]

            if database == "mobidblt" and seq_feature is None:
                seq_feature = "Consensus Disorder Prediction"

            store.update(protein_acc, {
                signature_acc: {
                    "database": database,
                    "locations": [(pos_start, pos_end, seq_feature)]
                }
            }, replace=True)

            i += 1
            if not i % 1000000:
                store.sync()

                if not i % 100000000:
                    logger.info(f"{i:>13,}")

        cur.close()
        con.close()

        logger.info(f"{i:>13,}")
        size = store.merge(processes=processes)
        logger.info(f"temporary files: {size/1024/1024:.0f} MB")


def export_matches(url: str, keyfile: str, output: str,
                   processes: int = 1, tmpdir: Optional[str] = None):
    logger.info("starting")
    with Store(output, Store.load_keys(keyfile), tmpdir) as store:
        con = cx_Oracle.connect(url)
        cur = con.cursor()
        cur.execute(
            """
            SELECT M.PROTEIN_AC, M.METHOD_AC, M.MODEL_AC, M.POS_FROM, 
                   M.POS_TO, M.FRAGMENTS, M.SCORE, E.ENTRY_AC
            FROM INTERPRO.MATCH M
            LEFT OUTER JOIN (
              SELECT E.ENTRY_AC, EM.METHOD_AC
              FROM INTERPRO.ENTRY E
              INNER JOIN INTERPRO.ENTRY2METHOD EM
                ON E.ENTRY_AC = EM.ENTRY_AC
              WHERE E.CHECKED = 'Y'
            ) E ON M.METHOD_AC = E.METHOD_AC
            """
        )

        i = 0
        for row in cur:
            if row[5]:
                fragments = []
                for frag in row[5].split(','):
                    # Format: START-END-STATUS
                    s, e, t = frag.split('-')
                    fragments.append({
                        "start": int(s),
                        "end": int(e),
                        "dc-status": DC_STATUSES[t]
                    })
            else:
                fragments = [{
                    "start": row[3],
                    "end": row[4],
                    "dc-status": DC_STATUSES['S']  # Continuous
                }]

            store.append(row[0], (
                row[1],     # signature
                row[2],     # model
                row[6],     # score
                fragments,
                row[7]      # InterPro entry
            ))

            i += 1
            if not i % 1000000:
                store.sync()

                if not i % 100000000:
                    logger.info(f"{i:>13,}")

        cur.close()
        con.close()

        logger.info(f"{i:>13,}")
        size = store.merge(fn=_post_matches, processes=processes)
        logger.info(f"temporary files: {size/1024/1024:.0f} MB")


def _post_matches(matches: Sequence[tuple]) -> dict:
    entries = {}
    signatures = {}
    for signature_acc, model, score, fragments, entry_acc in matches:
        fragments.sort(key=repr_fragment)

        try:
            s = signatures[signature_acc]
        except KeyError:
            s = signatures[signature_acc] = []

        s.append({
            "fragments": fragments,
            "model": model or signature_acc,
            "score": score
        })

        if entry_acc:
            try:
                entries[entry_acc].append(fragments)
            except KeyError:
                entries[entry_acc] = [fragments]

    for entry_acc, locations in entries.items():
        condensed = []
        for start, end in condense_locations(locations):
            condensed.append({
                "fragments": [{
                    "start": start,
                    "end": end,
                    "dc-status": DC_STATUSES['S']
                }],
                "model": None,
                "score": None
            })

        entries[entry_acc] = condensed

    for signature_acc, locations in signatures.items():
        # Sort locations using their leftmost fragment (fragments are sorted)
        locations.sort(key=lambda l: repr_fragment(l["fragments"][0]))
        entries[signature_acc] = locations

    return entries


def export_proteins(url: str, keyfile: str, output: str,
                    processes: int = 1, tmpdir: Optional[str] = None):
    logger.info("starting")
    with Store(output, Store.load_keys(keyfile), tmpdir) as store:
        con = cx_Oracle.connect(url)
        cur = con.cursor()
        cur.execute(
            """
            SELECT 
              PROTEIN_AC, NAME, DBCODE, LEN, FRAGMENT, 
              TO_CHAR(TAX_ID), CRC64
            FROM INTERPRO.PROTEIN
            """
        )

        i = 0
        for row in cur:
            store[row[0]] = {
                "identifier": row[1],
                "reviewed": row[2] == 'S',
                "length": row[3],
                "fragment": row[4] == 'Y',
                "taxid": row[5],
                "crc64": row[6]
            }

            i += 1
            if not i % 1000000:
                store.sync()

                if not i % 10000000:
                    logger.info(f"{i:>12,}")

        cur.close()
        con.close()

        logger.info(f"{i:>12,}")
        size = store.merge(processes=processes)
        logger.info(f"temporary files: {size/1024/1024:.0f} MB")


def export_residues(url: str, dt: DirectoryTree) -> List[str]:
    files = []

    con = cx_Oracle.connect(url)
    cur = con.cursor()
    cur.execute(
        """
        SELECT S.PROTEIN_AC, S.METHOD_AC, M.NAME, LOWER(D.DBSHORT),
               S.DESCRIPTION, S.RESIDUE, S.RESIDUE_START, S.RESIDUE_END
        FROM INTERPRO.SITE_MATCH S
        INNER JOIN INTERPRO.METHOD M ON S.METHOD_AC = M.METHOD_AC
        INNER JOIN INTERPRO.CV_DATABASE D ON M.DBCODE = D.DBCODE
        """
    )

    i = 0
    proteins = {}
    for row in cur:
        protein_acc = row[0]
        signature_acc = row[1]
        signature_name = row[2]
        database = row[3]
        description = row[4]
        residue = row[5]
        pos_start = row[6]
        pos_end = row[7]

        try:
            entries = proteins[protein_acc]
        except KeyError:
            entries = proteins[protein_acc] = {}

        try:
            entry = entries[signature_acc]
        except KeyError:
            entry = entries[signature_acc] = {
                "name": signature_name,
                "database": database,
                "descriptions": {}
            }

        try:
            fragments = entry["descriptions"][description]
        except KeyError:
            fragments = entry["descriptions"][description] = []

        fragments.append((residue, pos_start, pos_end))
        i += 1
        if not i % 1000000:
            files.append(dt.mktemp())
            with DumpFile(files[-1], compress=True) as df:
                for protein_acc in sorted(proteins):
                    df.dump((protein_acc, proteins[protein_acc]))

            proteins = {}

            if not i % 100000000:
                logger.info(f"{i:>15,}")

    logger.info(f"{i:>15,}")
    cur.close()
    con.close()

    files.append(dt.mktemp())
    with DumpFile(files[-1], compress=True) as df:
        for protein_acc in sorted(proteins):
            df.dump((protein_acc, proteins[protein_acc]))

    return files


def export_sequences(url: str, keyfile: str, output: str,
                     processes: int = 1, tmpdir: Optional[str] = None):
    logger.info("starting")
    with Store(output, Store.load_keys(keyfile), tmpdir) as store:
        con = cx_Oracle.connect(url)
        cur = con.cursor()
        cur.execute(
            """
            SELECT /*+ PARALLEL */ UX.AC, UP.SEQ_SHORT, UP.SEQ_LONG
            FROM UNIPARC.XREF UX
            INNER JOIN UNIPARC.PROTEIN UP ON UX.UPI = UP.UPI
            WHERE UX.DBID IN (2, 3)
            AND UX.DELETED = 'N'
            """
        )

        i = 0
        for row in cur:
            store[row[0]] = row[2].read() if row[2] is not None else row[1]

            i += 1
            if not i % 1000000:
                store.sync()

                if not i % 10000000:
                    logger.info(f"{i:>12,}")

        cur.close()
        con.close()

        logger.info(f"{i:>12,}")
        size = store.merge(processes=processes)
        logger.info(f"temporary files: {size/1024/1024:.0f} MB")


def get_isoforms(url: str) -> dict:
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    cur.execute(
        """
        SELECT EM.METHOD_AC, EM.ENTRY_AC 
        FROM INTERPRO.ENTRY2METHOD EM
        INNER JOIN INTERPRO.ENTRY E ON EM.ENTRY_AC = E.ENTRY_AC
        WHERE E.CHECKED = 'Y'
        """
    )
    integrated = dict(cur.fetchall())

    cur.execute(
        """
        SELECT V.PROTEIN_AC, V.VARIANT, V.LENGTH, V.CRC64, 
               P.SEQ_SHORT, P.SEQ_LONG
        FROM INTERPRO.VARSPLIC_MASTER V
        INNER JOIN UNIPARC.PROTEIN P 
          ON V.CRC64 = P.CRC64
        """
    )

    isoforms = {}
    for row in cur:
        variant_acc = row[0] + '-' + str(row[1])
        isoforms[variant_acc] = {
            "protein_acc": row[0],
            "length": row[2],
            "crc64": row[3],
            "sequence": row[5].read() if row[5] is not None else row[4],
            "matches": []
        }

    # PROTEIN_AC is actually PROTEIN-VARIANT (e.g. Q13733-1)
    cur.execute(
        """
        SELECT PROTEIN_AC, METHOD_AC, MODEL_AC, SCORE, POS_FROM, POS_TO, 
               FRAGMENTS
        FROM INTERPRO.VARSPLIC_MATCH M
        """
    )

    for row in cur:
        try:
            isoform = isoforms[row[0]]
        except KeyError:
            continue

        if row[6]:
            fragments = []
            for frag in row[6].split(','):
                # Format: START-END-STATUS
                s, e, t = frag.split('-')
                fragments.append({
                    "start": int(s),
                    "end": int(e),
                    "dc-status": DC_STATUSES[t]
                })
        else:
            fragments = [{
                "start": row[4],
                "end": row[5],
                "dc-status": DC_STATUSES['S']  # Continuous
            }]

        signature_acc = row[1]
        isoform["matches"].append((
            signature_acc,                  # signature
            row[2],                         # model
            row[3],                         # score
            fragments,
            integrated.get(signature_acc)   # InterPro entry
        ))

    cur.close()
    con.close()

    for variant in isoforms.values():
        variant["matches"] = _post_matches(variant["matches"])

    return isoforms
