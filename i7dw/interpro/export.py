#!/usr/bin/env python
# -*- coding: utf-8 -*-

import json
import logging

from .. import dbms, io

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s: %(levelname)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)


def chunk_proteins(uri: str, dst: str, order_by: bool=True,
                   chunk_size: int=100000):
    chunks = []
    con, cur = dbms.connect(uri)

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


def export_protein2matches(uri, src, dst, tmpdir=None, processes=1,
                           sync_frequency=1000000):
    logging.info("starting")

    with open(src, "rt") as fh:
        keys = json.load(fh)

    with io.Store(dst, keys, tmpdir) as store:
        con, cur = dbms.connect(uri)
        cur.execute(
            """
            SELECT
              M.PROTEIN_AC, M.METHOD_AC, M.MODEL_AC, NULL,
              M.POS_FROM, M.POS_TO, M.FRAGMENTS
            FROM INTERPRO.MATCH M
            UNION ALL
            SELECT
              FM.PROTEIN_AC, FM.METHOD_AC, NULL, FM.SEQ_FEATURE,
              FM.POS_FROM, FM.POS_TO, NULL
            FROM INTERPRO.FEATURE_MATCH FM
            WHERE FM.DBCODE = 'g'
            """
        )


        i = 0
        for row in cur:
            protein_acc = row[0]
            method_acc = row[1]
            model_acc = row[2]
            seq_feature = row[3]
            pos_start = row[4]
            pos_end = row[5]
            fragments_str = row[6]

            if fragments_str is None:
                fragments = [{"start": pos_start, "end": pos_end}]
            else:
                # Discontinuous domains
                fragments = []
                for frag in fragments_str.split(','):
                    """
                    Format: START-END-TYPE
                    Types:
                        * S: Continuous single chain domain
                        * N: N-terminal discontinuous
                        * C: C-terminal discontinuous
                        * NC: N and C -terminal discontinuous
                    """
                    s, e, t = frag.split('-')
                    s = int(s)
                    e = int(e)
                    if s < e:
                        fragments.append({
                            "start": s,
                            "end": e
                        })

                if not fragments:
                    # Fallback to match start/end positions
                    fragments.append({"start": pos_start, "end": pos_end})

            if model_acc == method_acc:
                model_acc = None

            store.append(protein_acc, {
                "method_ac": method_acc,
                "model_ac": model_acc,
                "seq_feature": seq_feature,
                "fragments": fragments
            })

            i += 1
            if sync_frequency and not i % sync_frequency:
                store.sync()

            if not i % 100000000:
                logging.info("{:>15,}".format(i))

        cur.close()
        con.close()
        logging.info("{:>15,}".format(i))
        store.merge(func=sort_matches, processes=processes)
        logging.info("temporary files: {:,} bytes".format(store.size))


def sort_fragments(fragments: list) -> tuple:
    start = end = None
    for f in fragments:
        if start is None or f["start"] < start:
            start = f["start"]

        if end is None or f["end"] < end:
            end = f["end"]

    return start, end


def sort_matches(matches: list) -> list:
    return sorted(matches, key=lambda m: sort_fragments(m["fragments"]))


def export_protein2features(uri, src, dst, tmpdir=None, processes=1,
                            sync_frequency=1000000):
    logging.info("starting")

    with open(src, "rt") as fh:
        keys = json.load(fh)

    with io.Store(dst, keys, tmpdir) as store:
        con, cur = dbms.connect(uri)
        cur.execute(
            """
            SELECT
              FM.PROTEIN_AC, FM.METHOD_AC, LOWER(DB.DBSHORT),
              FM.POS_FROM, FM.POS_TO
            FROM INTERPRO.FEATURE_MATCH FM
            INNER JOIN INTERPRO.CV_DATABASE DB ON FM.DBCODE = DB.DBCODE
            WHERE FM.DBCODE != 'g'
            """
        )

        i = 0
        for protein_acc, method_acc, database, start, end in cur:
            store.update(
                protein_acc,
                {
                    method_acc: {
                        "accession": method_acc,
                        "source_database": database,
                        "locations": [{"start": start, "end": end}]
                    }
                }
            )

            i += 1
            if sync_frequency and not i % sync_frequency:
                store.sync()

            if not i % 100000000:
                logging.info("{:>15,}".format(i))

        cur.close()
        con.close()
        logging.info("{:>15,}".format(i))
        store.merge(func=sort_feature_locations, processes=processes)
        logging.info("temporary files: {:,} bytes".format(store.size))


def sort_feature_locations(item: dict) -> dict:
    for method in item.values():
        method["locations"] = [{
            "fragments": sorted(method["locations"],
                                key=lambda x: (x["start"], x["end"]))
        }]
    return item


def export_protein2residues(uri, src, dst, tmpdir=None, processes=1,
                            sync_frequency=1000000):
    logging.info("starting")

    with open(src, "rt") as fh:
        keys = json.load(fh)

    with io.Store(dst, keys, tmpdir) as store:
        con, cur = dbms.connect(uri)
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
                logging.info("{:>15,}".format(i))

        cur.close()
        con.close()
        logging.info("{:>15,}".format(i))
        store.merge(func=sort_residues, processes=processes)
        logging.info("temporary files: {:,} bytes".format(store.size))


def sort_residues(item: dict) -> dict:
    for method in item.values():
        locations = []
        for loc in method["locations"].values():
            loc["fragments"].sort(key=lambda x: (x["start"], x["end"]))
            locations.append(loc)

        method["locations"] = locations

    return item


def export_proteins(uri, src, dst, tmpdir=None, processes=1,
                    sync_frequency=1000000):
    logging.info("starting")

    with open(src, "rt") as fh:
        keys = json.load(fh)

    with io.Store(dst, keys, tmpdir) as store:
        con, cur = dbms.connect(uri)
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
                "isReviewed": dbcode == 'S',
                "isFrag": frag == 'Y',
                "length": length
            }

            i += 1
            if sync_frequency and not i % sync_frequency:
                store.sync()

            if not i % 10000000:
                logging.info("{:>12,}".format(i))

        cur.close()
        con.close()
        logging.info("{:>12,}".format(i))
        store.merge(processes=processes)
        logging.info("temporary files: {:,} bytes".format(store.size))


def export_sequences(uri, src, dst, tmpdir=None, processes=1,
                     sync_frequency=1000000):
    logging.info("starting")

    with open(src, "rt") as fh:
        keys = json.load(fh)

    with io.Store(dst, keys, tmpdir) as store:
        con, cur = dbms.connect(uri)
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
                logging.info("{:>12,}".format(i))

        cur.close()
        con.close()
        logging.info("{:>12,}".format(i))
        store.merge(processes=processes)
        logging.info("temporary files: {:,} bytes".format(store.size))
