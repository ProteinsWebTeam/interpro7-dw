import hashlib
import json
from typing import Optional, Tuple

from . import mysql
from .. import dbms, io, logger


DC_STATUSES = {
    # Continuous single chain domain
    "S": "CONTINUOUS",
    # N terminus discontinuous
    "N": "N_TERMINAL_DISC",
    # C terminus discontinuous
    "C": "C_TERMINAL_DISC",
    # N and C terminus discontinuous
    "NC": "NC_TERMINAL_DISC"
}


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
    logger.info("starting")

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
                "seq_feature": seq_feature,
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
        logger.info("temporary files: {:.0f} MB".format(store.size/1024/1024))


def repr_frag(f: dict) -> Tuple[int, int]:
    return f["start"], f["end"]


def sort_matches(matches: list) -> list:
    return sorted(matches, key=lambda m: repr_frag(m["fragments"][0]))


def export_protein2features(uri, src, dst, tmpdir=None, processes=1,
                            sync_frequency=1000000):
    logger.info("starting")

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
                logger.info("{:>15,}".format(i))

        cur.close()
        con.close()
        logger.info("{:>15,}".format(i))
        store.merge(func=sort_feature_locations, processes=processes)
        logger.info("temporary files: {:.0f} MB".format(store.size/1024/1024))


def sort_feature_locations(item: dict) -> dict:
    for method in item.values():
        method["locations"].sort(key=repr_frag)

    return item


def export_protein2residues(uri, src, dst, tmpdir=None, processes=1,
                            sync_frequency=1000000):
    logger.info("starting")

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
                logger.info("{:>15,}".format(i))

        cur.close()
        con.close()
        logger.info("{:>15,}".format(i))
        store.merge(func=sort_residues, processes=processes)
        logger.info("temporary files: {:.0f} MB".format(store.size/1024/1024))


def sort_residues(item: dict) -> dict:
    for method in item.values():
        locations = []
        for loc in method["locations"].values():
            loc["fragments"].sort(key=repr_frag)
            locations.append(loc)

        locations.sort(key=lambda m: repr_frag(m["fragments"][0]))
        method["locations"] = locations

    return item


def export_proteins(uri, src, dst, tmpdir=None, processes=1,
                    sync_frequency=1000000):
    logger.info("starting")

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
                logger.info("{:>12,}".format(i))

        cur.close()
        con.close()
        logger.info("{:>12,}".format(i))
        store.merge(processes=processes)
        logger.info("temporary files: {:.0f} MB".format(store.size/1024/1024))


def export_sequences(uri, src, dst, tmpdir=None, processes=1,
                     sync_frequency=1000000):
    logger.info("starting")

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
                logger.info("{:>12,}".format(i))

        cur.close()
        con.close()
        logger.info("{:>12,}".format(i))
        store.merge(processes=processes)
        logger.info("temporary files: {:.0f} MB".format(store.size/1024/1024))


def export_ida(my_uri: str, src_matches: str, dst_ida: str,
               tmpdir: Optional[str]=None, processes: int=1,
               sync_frequency: int=1000000):

    logger.info("starting")
    pfam_entries = {}
    for e in mysql.entry.get_entries(my_uri).values():
        if e["database"] == "pfam":
            pfam_ac = e["accession"]
            interpro_ac = e["integrated"]
            pfam_entries[pfam_ac] = interpro_ac

    with io.Store(src_matches) as src, io.Store(dst_ida, src.keys, tmpdir) as dst:
        i = 0

        for acc, matches in src:
            dom_arch = []
            for m in matches:
                method_ac = m["method_ac"]

                if method_ac in pfam_entries:

                    interpro_ac = pfam_entries[method_ac]
                    if interpro_ac:
                        dom_arch.append("{}:{}".format(method_ac, interpro_ac))
                    else:
                        dom_arch.append("{}".format(method_ac))

            if dom_arch:
                ida = '-'.join(dom_arch)
                ida_id = hashlib.sha1(ida.encode("utf-8")).hexdigest()
                dst[acc] = (ida, ida_id)

            i += 1
            if sync_frequency and not i % sync_frequency:
                dst.sync()

            if not i % 10000000:
                logger.info("{:>12,}".format(i))

        logger.info("{:>12,}".format(i))
        dst.merge(processes=processes)
        logger.info("temporary files: {:.0f} MB".format(dst.size/1024/1024))


def export_isoforms(uri, src, dst, tmpdir=None, processes=1,
                    sync_frequency=1000000):
    logger.info("starting")

    with open(src, "rt") as fh:
        keys = json.load(fh)

    with io.Store(dst, keys, tmpdir) as store:
        con, cur = dbms.connect(uri)
        cur.execute(
            """
            SELECT
              XV.PROTEIN_AC,
              XV.VARIANT,
              P.LEN,
              P.SEQ_SHORT,
              P.SEQ_LONG,
              MA.METHOD_AC,
              MA.SEQ_START,
              MA.SEQ_END,
              MA.FRAGMENTS,
              MA.MODEL_AC
            FROM (
              SELECT
                SUBSTR(AC, 1, INSTR(AC, '-') - 1) AS PROTEIN_AC,
                SUBSTR(AC, INSTR(AC, '-') + 1) AS VARIANT,
                UPI
              FROM UNIPARC.XREF
              WHERE DBID IN (24, 25) 
                AND DELETED = 'N'
            ) XV
            INNER JOIN UNIPARC.PROTEIN P
              ON XV.UPI = P.UPI
            INNER JOIN UNIPARC.XREF XP
              ON XV.PROTEIN_AC = XP.AC
                AND XP.DBID IN (2, 3)
                AND XP.DELETED = 'N'
            INNER JOIN IPRSCAN.MV_IPRSCAN MA
              ON XV.UPI = MA.UPI
            INNER JOIN INTERPRO.IPRSCAN2DBCODE I2D
              ON MA.ANALYSIS_ID = I2D.IPRSCAN_SIG_LIB_REL_ID
            WHERE XV.UPI != XP.UPI
            """
            # WHERE I2D.DBCODE NOT IN ('g', 'j', 'n', 'q', 's', 'v', 'x')
        )

        i = 0
        for row in cur:
            if row[8] is None:
                fragments = [{
                    "start": row[6],
                    "end": row[7],
                    "dc-status": "CONTINUOUS"
                }]
            else:
                fragments = []
                for frag in row[8].split(','):
                    # Format: START-END-STATUS
                    s, e, t = frag.split('-')
                    fragments.append({
                        "start": int(s),
                        "end": int(e),
                        "dc-status": DC_STATUSES[t]
                    })
                fragments.sort(key=repr_frag)

            isoform = int(row[1])
            store.update(
                row[0],
                {
                    isoform: {
                        "isoform": isoform,
                        "length": row[2],
                        "sequence": row[3] if row[3] else row[4].read(),
                        "locations": [{
                            "method_acc": row[5],
                            "fragments": fragments,
                            "model_ac": row[9]
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
            store.merge(func=sort_isoforms, processes=processes)
            logger.info("temporary files: {:.0f} MB".format(store.size/1024/1024))


def sort_isoforms(item: dict) -> list:
    for isoform in item.values():
        isoform["locations"].sort(key=lambda l: repr_frag(l["fragments"][0]))

    return sorted(item.values(), key=lambda x: x["isoform"])
