# -*- coding: utf-8 -*-

import json
from typing import Optional

import cx_Oracle

from i7dw import io, logger
from i7dw.interpro import MIN_OVERLAP, extract_frag
from .utils import DC_STATUSES


def chunk_proteins(url: str, dst: str, order_by: bool=True,
                   chunk_size: int=100000):
    chunks = []

    con = cx_Oracle.connect(url)
    cur = con.cursor()
    if order_by:
        cur.execute(
            """
            SELECT PROTEIN_AC
            FROM INTERPRO.PROTEIN
            ORDER BY PROTEIN_AC
            """
        )

        cnt = 0
        for row in cur:
            cnt += 1
            if cnt % chunk_size == 1:
                chunks.append(row[0])
    else:
        cur.execute(
            """
            SELECT PROTEIN_AC
            FROM INTERPRO.PROTEIN
            """
        )

        proteins = [row[0] for row in cur]
        proteins.sort()

        for i in range(0, len(proteins), chunk_size):
            chunks.append(proteins[i])

    cur.close()
    con.close()

    with open(dst, "wt") as fh:
        json.dump(chunks, fh)


def export_matches(url: str, src: str, dst: str,
                   tmpdir: Optional[str]=None, processes: int=1,
                   sync_frequency: int=1000000):
    logger.info("starting")

    with open(src, "rt") as fh:
        keys = json.load(fh)

    with io.Store(dst, keys, tmpdir) as store:
        con = cx_Oracle.connect(url)
        cur = con.cursor()

        cur.execute("SELECT METHOD_AC, ENTRY_AC FROM INTERPRO.ENTRY2METHOD")
        entries = dict(cur.fetchall())

        cur.execute(
            """
            SELECT
              PROTEIN_AC, METHOD_AC, MODEL_AC, POS_FROM, POS_TO, FRAGMENTS
            FROM INTERPRO.MATCH
            """
        )

        i = 0
        for row in cur:
            protein_acc = row[0]
            method_acc = row[1]
            model_acc = row[2]
            pos_start = row[3]
            pos_end = row[4]
            fragments_str = row[5]

            if fragments_str is None:
                fragments = [{
                    "start": pos_start,
                    "end": pos_end,
                    "dc-status": "CONTINUOUS"
                }]
            else:
                fragments = []
                for frag in fragments_str.split(','):
                    # Format: START-END-STATUS
                    s, e, t = frag.split('-')
                    fragments.append({
                        "start": int(s),
                        "end": int(e),
                        "dc-status": DC_STATUSES[t]
                    })

            store.update(protein_acc, {
                method_acc: {
                    "condense": False,
                    "locations": [{
                        "model": model_acc if model_acc != method_acc else None,
                        "fragments": fragments
                    }]
                }
            })

            try:
                entry_acc = entries[method_acc]
            except KeyError:
                pass
            else:
                store.update(protein_acc, {
                    entry_acc: {
                        "condense": True,
                        "locations": [{
                            "model": None,
                            "fragments": fragments
                        }]
                    }
                })

            i += 1
            if sync_frequency and not i % sync_frequency:
                store.sync()

            if not i % 100000000:
                logger.info(f"{i:>15,}")

        cur.close()
        con.close()

        logger.info(f"{i:>15,}")
        store.merge(func=sort_matches, processes=processes)
        logger.info(f"temporary files: {store.size/1024/1024:.0f} MB")


def sort_matches(item: dict) -> dict:
    for entry in item.values():
        for loc in entry["locations"]:
            loc["fragments"].sort(key=extract_frag)

        # Fragments are sorted, so use the leftmost fragment to sort locations
        entry["locations"].sort(key=lambda l: extract_frag(l["fragments"][0]))

        if entry["condense"]:
            start = end = None
            locations = []

            for loc in entry["locations"]:
                # We do not consider fragmented matches
                s = loc["fragments"][0]["start"]
                e = loc["fragments"][0]["start"]

                if start is None:
                    start, end = s, e
                    continue
                elif s <= end:
                    # Locations are overlapping (at least one residue)
                    overlap = min(end, e) - max(start, s) + 1
                    shortest = min(end - start, e - s) + 1

                    if overlap >= shortest * MIN_OVERLAP:
                        # Merge
                        end = e
                        continue

                locations.append((start, end))
                start, end = s, e

            # Adding last location
            locations.append((start, end))

            # Replacing current locations with condensed
            entry["locations"] = [
                {
                    "model": None,
                    "fragments": [{
                        "start": start,
                        "end": end,
                        "dc-status": "CONTINUOUS"
                    }]
                }
                for start, end in locations
            ]

    return item


def export_features(url: str, src: str, dst: str,
                    tmpdir: Optional[str]=None, processes: int=1,
                    sync_frequency: int=1000000):
    logger.info("starting")

    with open(src, "rt") as fh:
        keys = json.load(fh)

    with io.Store(dst, keys, tmpdir) as store:
        con = cx_Oracle.connect(url)
        cur = con.cursor()
        cur.execute(
            """
            SELECT
              FM.PROTEIN_AC, FM.METHOD_AC, LOWER(DB.DBSHORT),
              FM.POS_FROM, FM.POS_TO, FM.SEQ_FEATURE
            FROM INTERPRO.FEATURE_MATCH FM
            INNER JOIN INTERPRO.CV_DATABASE DB ON FM.DBCODE = DB.DBCODE
            """
        )

        i = 0
        for protein_acc, method_acc, database, start, end, seq_feature in cur:
            store.update(
                protein_acc,
                {
                    method_acc: {
                        "accession": method_acc,
                        "source_database": database,
                        "locations": [{
                            "fragments": [{
                                "start": start,
                                "end": end,
                                "seq_feature": seq_feature
                            }]
                        }]
                    }
                }
            )

            i += 1
            if sync_frequency and not i % sync_frequency:
                store.sync()

            if not i % 100000000:
                logger.info(f"{i:>15,}")

        cur.close()
        con.close()

        logger.info(f"{i:>15,}")
        store.merge(func=sort_features, processes=processes)
        logger.info(f"temporary files: {store.size/1024/1024:.0f} MB")


def sort_features(item: dict) -> dict:
    for entry in item.values():
        # Only one fragment per location
        entry["locations"].sort(key=lambda l: extract_frag(l["fragments"][0]))

    return item


def export_residues(url: str, src: str, dst: str,
                    tmpdir: Optional[str]=None, processes: int=1,
                    sync_frequency: int=1000000):
    logger.info("starting")

    with open(src, "rt") as fh:
        keys = json.load(fh)

    with io.Store(dst, keys, tmpdir) as store:
        con = cx_Oracle.connect(url)
        cur = con.cursor()
        cur.execute(
            """
            SELECT
              S.PROTEIN_AC, S.METHOD_AC, M.NAME, LOWER(D.DBSHORT),
              S.DESCRIPTION, S.RESIDUE, S.RESIDUE_START, S.RESIDUE_END
            FROM INTERPRO.SITE_MATCH S
            INNER JOIN INTERPRO.METHOD M ON S.METHOD_AC = M.METHOD_AC
            INNER JOIN INTERPRO.CV_DATABASE D ON M.DBCODE = D.DBCODE
            """
        )

        i = 0
        for row in cur:
            protein_acc = row[0]
            method_acc = row[1]
            method_name = row[2]
            database = row[3]
            description = row[4]
            residue = row[5]
            start = row[6]
            end = row[7]

            store.update(
                protein_acc,
                {
                    method_acc: {
                        "accession": method_acc,
                        "name": method_name,
                        "source_database": database,
                        "locations": {
                            description: [{
                                "residues": residue,
                                "start": start,
                                "end": end
                            }]
                        }
                    }
                }
            )

            i += 1
            if sync_frequency and not i % sync_frequency:
                store.sync()

            if not i % 100000000:
                logger.info(f"{i:>15,}")

        cur.close()
        con.close()

        logger.info(f"{i:>15,}")
        store.merge(func=sort_residues, processes=processes)
        logger.info(f"temporary files: {store.size/1024/1024:.0f} MB")


def sort_residues(item: dict) -> dict:
    for entry in item.values():
        locations = []

        for description, fragments in entry["locations"].items():
            locations.append({
                "description": description,
                "fragments": sorted(fragments, key=extract_frag)
            })

        entry["locations"] = sorted(locations, key=lambda l: extract_frag(l["fragments"][0]))

    return item


def export_proteins(url: str, src: str, dst: str,
                    tmpdir: Optional[str]=None, processes: int=1,
                    sync_frequency: int=1000000):
    logger.info("starting")

    with open(src, "rt") as fh:
        keys = json.load(fh)

    with io.Store(dst, keys, tmpdir) as store:
        con = cx_Oracle.connect(url)
        cur = con.cursor()
        # TODO: JOIN with TAXONOMY.V_PUBLIC_NODE@SWPREAD instead of ETAXI
        cur.execute(
            """
            SELECT
              PROTEIN_AC, TO_CHAR(TAX_ID), NAME, DBCODE, FRAGMENT, LEN
            FROM INTERPRO.PROTEIN P
            WHERE TAX_ID IN (SELECT TAX_ID FROM INTERPRO.ETAXI)
            """
        )

        i = 0
        for acc, tax_id, name, dbcode, frag, length in cur:
            store[acc] = {
                "identifier": name,
                "is_reviewed": dbcode == 'S',
                "is_fragment": frag == 'Y',
                "length": length,
                "taxon": tax_id
            }

            i += 1
            if sync_frequency and not i % sync_frequency:
                store.sync()

            if not i % 10000000:
                logger.info(f"{i:>12,}")

        cur.close()
        con.close()

        logger.info(f"{i:>12,}")
        store.merge(processes=processes)
        logger.info(f"temporary files: {store.size/1024/1024:.0f} MB")


def export_sequences(url: str, src: str, dst: str,
                     tmpdir: Optional[str]=None, processes: int=1,
                     sync_frequency: int=1000000):
    logger.info("starting")

    with open(src, "rt") as fh:
        keys = json.load(fh)

    with io.Store(dst, keys, tmpdir) as store:
        con = cx_Oracle.connect(url)
        cur = con.cursor()
        cur.execute(
            """
            SELECT
              UX.AC, UP.SEQ_SHORT, UP.SEQ_LONG
            FROM UNIPARC.XREF UX
            INNER JOIN UNIPARC.PROTEIN UP ON UX.UPI = UP.UPI
            WHERE UX.DBID IN (2, 3)
            AND UX.DELETED = 'N'
            """
        )

        i = 0
        for acc, seq_short, seq_long in cur:
            if seq_long is not None:
                store[acc] = seq_long.read()
            else:
                store[acc] = seq_short

            i += 1
            if sync_frequency and not i % sync_frequency:
                store.sync()

            if not i % 10000000:
                logger.info(f"{i:>12,}")

        cur.close()
        con.close()

        logger.info(f"{i:>12,}")
        store.merge(processes=processes)
        logger.info(f"temporary files: {store.size/1024/1024:.0f} MB")
