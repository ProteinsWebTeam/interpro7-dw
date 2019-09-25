# -*- coding: utf-8 -*-

import json

import cx_Oracle

from i7dw import io, logger
from i7dw.interpro import repr_frag
from . import DC_STATUSES


def chunk_proteins(url: str, dst: str, order_by: bool = True,
                   chunk_size: int = 100000):
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


def export_protein2matches(url, src, dst, tmpdir=None, processes=1,
                           sync_frequency=1000000):
    logger.info("starting")

    with open(src, "rt") as fh:
        keys = json.load(fh)

    with io.Store(dst, keys, tmpdir) as store:
        con = cx_Oracle.connect(url)
        cur = con.cursor()
        cur.execute(
            """
            SELECT
              M.PROTEIN_AC, M.METHOD_AC, M.MODEL_AC,
              M.POS_FROM, M.POS_TO, M.FRAGMENTS
            FROM INTERPRO.MATCH M
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
                fragments.sort(key=repr_frag)

            store.append(protein_acc, {
                "method_ac": method_acc,
                "model_ac": model_acc if model_acc != method_acc else None,
                "fragments": fragments
            })

            i += 1
            if sync_frequency and not i % sync_frequency:
                store.sync()

            if not i % 100000000:
                logger.info("{:>15,}".format(i))

        cur.close()
        con.close()
        logger.info("{:>15,}".format(i))
        store.merge(func=sort_matches, processes=processes)
        logger.info(
            "temporary files: {:.0f} MB".format(store.size / 1024 / 1024))


def sort_matches(matches: list) -> list:
    return sorted(matches, key=lambda m: repr_frag(m["fragments"][0]))


def export_protein2features(url, src, dst, tmpdir=None, processes=1,
                            sync_frequency=1000000):
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
                            "start": start,
                            "end": end,
                            "seq_feature": seq_feature
                        }]
                    }
                }
            )

            i += 1
            if sync_frequency and not i % sync_frequency:
                store.sync()

            if not i % 100000000:
                logger.info("{:>15,}".format(i))

        cur.close()
        con.close()
        logger.info("{:>15,}".format(i))
        store.merge(func=sort_feature_locations, processes=processes)
        logger.info(
            "temporary files: {:.0f} MB".format(store.size / 1024 / 1024))


def sort_feature_locations(item: dict) -> dict:
    for method in item.values():
        locations = []
        for loc in sorted(method["locations"], key=repr_frag):
            locations.append({"fragments": [loc]})
        method["locations"] = locations

    return item


def export_protein2residues(url, src, dst, tmpdir=None, processes=1,
                            sync_frequency=1000000):
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
                            description: {
                                "description": description,
                                "fragments": [{
                                    "residues": residue,
                                    "start": start,
                                    "end": end
                                }]
                            }
                        }
                    }
                }
            )

            i += 1
            if sync_frequency and not i % sync_frequency:
                store.sync()

            if not i % 100000000:
                logger.info("{:>15,}".format(i))

        cur.close()
        con.close()
        logger.info("{:>15,}".format(i))
        store.merge(func=sort_residues, processes=processes)
        logger.info(
            "temporary files: {:.0f} MB".format(store.size / 1024 / 1024))


def sort_residues(item: dict) -> dict:
    for method in item.values():
        locations = []
        for loc in method["locations"].values():
            loc["fragments"].sort(key=repr_frag)
            locations.append(loc)

        locations.sort(key=lambda m: repr_frag(m["fragments"][0]))
        method["locations"] = locations

    return item


def export_proteins(url, src, dst, tmpdir=None, processes=1,
                    sync_frequency=1000000):
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
              PROTEIN_AC,
              TO_CHAR(TAX_ID),
              NAME,
              DBCODE,
              FRAGMENT,
              LEN
            FROM INTERPRO.PROTEIN P
            WHERE TAX_ID IN (
                SELECT TAX_ID
                FROM INTERPRO.ETAXI
            )
            """
        )

        i = 0
        for acc, tax_id, name, dbcode, frag, length in cur:
            store[acc] = {
                "taxon": tax_id,
                "identifier": name,
                "is_reviewed": dbcode == 'S',
                "is_fragment": frag == 'Y',
                "length": length
            }

            i += 1
            if sync_frequency and not i % sync_frequency:
                store.sync()

            if not i % 10000000:
                logger.info("{:>12,}".format(i))

        cur.close()
        con.close()
        logger.info("{:>12,}".format(i))
        store.merge(processes=processes)
        logger.info(
            "temporary files: {:.0f} MB".format(store.size / 1024 / 1024))


def export_sequences(url, src, dst, tmpdir=None, processes=1,
                     sync_frequency=1000000):
    logger.info("starting")

    with open(src, "rt") as fh:
        keys = json.load(fh)

    with io.Store(dst, keys, tmpdir) as store:
        con = cx_Oracle.connect(url)
        cur = con.cursor()
        cur.execute(
            """
            SELECT
              UX.AC,
              UP.SEQ_SHORT,
              UP.SEQ_LONG
            FROM UNIPARC.XREF UX
            INNER JOIN UNIPARC.PROTEIN UP
              ON UX.UPI = UP.UPI
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
                logger.info("{:>12,}".format(i))

        cur.close()
        con.close()
        logger.info("{:>12,}".format(i))
        store.merge(processes=processes)
        logger.info("temporary files: {:.0f} MB".format(store.size / 1024 / 1024))