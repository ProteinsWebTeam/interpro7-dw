# -*- coding: utf-8 -*-

from typing import Optional, Sequence

import cx_Oracle

from interpro7dw import logger
from interpro7dw.utils import Store
from interpro7dw.ebi.interpro.utils import DC_STATUSES
from interpro7dw.ebi.interpro.utils import condense_locations
from interpro7dw.ebi.interpro.utils import repr_fragment


def chunk_proteins(url: str, keyfile: str, chunk_size: int=50000):
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
                    dir: Optional[str]=None, processes: int=1):
    logger.info("starting")
    with Store(output, Store.load_keys(keyfile), dir) as store:
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

            store.append(protein_acc, (signature_acc, database, pos_start,
                                       pos_end, seq_feature))

            i += 1
            if not i % 1000000:
                store.sync()

                if not i % 100000000:
                    logger.info(f"{i:>13,}")

        cur.close()
        con.close()

        logger.info(f"{i:>13,}")
        size = store.merge(fn=_post_features, processes=processes)
        logger.info(f"temporary files: {size/1024/1024:.0f} MB")


def _post_features(matches: Sequence[dict]) -> dict:
    entries = {}
    for acc, database, pos_start, pos_end, seq_feature in matches:
        try:
            obj = entries[acc]
        except KeyError:
            obj = entries[acc] = {
                "accession": acc,
                "source_database": database,
                "locations": []
            }
        finally:
            obj["locations"].append({
                "fragments": [{
                    "start": pos_start,
                    "end": pos_end,
                    "seq_feature": seq_feature
                }]
            })

    for obj in entries.values():
        # Only one fragment per location
        obj["locations"].sort(key=lambda l: repr_fragment(l["fragments"][0]))

    return entries


def export_matches(url: str, keyfile: str, output: str,
                   dir: Optional[str]=None, processes: int=1):
    logger.info("starting")
    with Store(output, Store.load_keys(keyfile), dir) as store:
        con = cx_Oracle.connect(url)
        cur = con.cursor()
        cur.execute(
            """
            SELECT M.PROTEIN_AC, M.METHOD_AC, M.MODEL_AC, M.POS_FROM, 
                   M.POS_TO, FRAGMENTS, E.ENTRY_AC
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
                row[1],
                row[2],
                row[6],
                fragments
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
    for signature_acc, model, entry_acc, fragments in matches:
        fragments.sort(key=repr_fragment)

        try:
            s = signatures[signature_acc]
        except KeyError:
            s = signatures[signature_acc] = []

        s.append({
            "fragments": fragments,
            "model_acc": model or signature_acc
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
                "model_acc": None
            })

        entries[entry_acc] = condensed

    for signature_acc, locations in signatures.items():
        # Sort locations using their leftmost fragment (fragments are sorted)
        locations.sort(key=lambda l: repr_fragment(l["fragments"][0]))
        entries[signature_acc] = locations

    return entries


def export_proteins(url: str, keyfile: str, output: str,
                    dir: Optional[str]=None, processes: int=1):
    logger.info("starting")
    with Store(output, Store.load_keys(keyfile), dir) as store:
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


def export_residues(url: str, keyfile: str, output: str,
                    dir: Optional[str]=None, processes: int=1):
    logger.info("starting")
    with Store(output, Store.load_keys(keyfile), dir) as store:
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
        for row in cur:
            protein_acc = row[0]
            signature_acc = row[1]
            signature_name = row[2]
            database = row[3]
            description = row[4]
            residue = row[5]
            pos_start = row[6]
            pos_end = row[7]

            store.append(protein_acc, (signature_acc, signature_name,
                                       database, description, residue,
                                       pos_start, pos_end))

            i += 1
            if not i % 1000000:
                store.sync()

                if not i % 100000000:
                    logger.info(f"{i:>13,}")

        cur.close()
        con.close()

        logger.info(f"{i:>13,}")
        size = store.merge(fn=_post_residues, processes=processes)
        logger.info(f"temporary files: {size/1024/1024:.0f} MB")


def _post_residues(matches: Sequence[dict]) -> dict:
    entries = {}
    for acc, name, database, descr, residue, pos_start, pos_end in matches:
        try:
            entry = entries[acc]
        except KeyError:
            entry = entries[acc] = {
                "accession": acc,
                "name": name,
                "source_database": database,
                "locations": {}
            }

        try:
            fragments = entry["locations"][descr]
        except KeyError:
            fragments = entry["locations"][descr] = []
        finally:
            fragments.append({
                "residues": residue,
                "start": pos_start,
                "end": pos_end
            })

    for entry in entries.values():
        locations = []
        for descr, fragments in entry["locations"].items():
            locations.append({
                "description": descr,
                "fragments": sorted(fragments, key=repr_fragment)
            })

        locations.sort(key=lambda l: repr_fragment(l["fragments"][0]))
        entry["locations"] = locations

    return entries


def export_sequences(url: str, keyfile: str, output: str,
                     dir: Optional[str]=None, processes: int=1):
    logger.info("starting")
    with Store(output, Store.load_keys(keyfile), dir) as store:
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


def get_isoforms(url: str):
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
        SELECT V.PROTEIN_AC, V.VARIANT, V.LENGTH, P.SEQ_SHORT, P.SEQ_LONG
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
            "sequence": row[4].read() if row[4] is not None else row[3],
            "matches": []
        }

    # PROTEIN_AC is actually PROTEIN-VARIANT (e.g. Q13733-1)
    cur.execute(
        """
        SELECT PROTEIN_AC, METHOD_AC, MODEL_AC, POS_FROM, POS_TO, FRAGMENTS
        FROM INTERPRO.VARSPLIC_MATCH M
        """
    )

    for row in cur:
        try:
            isoform = isoforms[row[0]]
        except KeyError:
            continue

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

        signature_acc = row[1]
        isoform["matches"].append((
            signature_acc,
            row[2] or signature_acc,
            integrated.get(signature_acc),
            fragments
        ))

    cur.close()
    con.close()

    for accession, variant in isoforms.items():
        yield (
            accession,
            variant["protein_acc"],
            variant["length"],
            variant["sequence"],
            _post_matches(variant["matches"])
        )
