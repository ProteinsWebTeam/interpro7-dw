import json
import pickle

import oracledb

from interpro7dw import pfam
from interpro7dw.utils import logger
from interpro7dw.utils.oracle import lob_as_str
from interpro7dw.utils.store import BasicStore


def get_clans(cur: oracledb.Cursor) -> dict[str, dict]:
    cur.execute(
        """
        SELECT
          C.CLAN_AC, C.NAME, C.DESCRIPTION, D.DBSHORT, CM.MEMBER_AC,
          M.NAME, M.DESCRIPTION, CM.LEN, CM.SCORE
        FROM INTERPRO.CLAN C
        INNER JOIN INTERPRO.CV_DATABASE D
          ON C.DBCODE = D.DBCODE
        INNER JOIN INTERPRO.CLAN_MEMBER CM
          ON C.CLAN_AC = CM.CLAN_AC
        INNER JOIN INTERPRO.METHOD M 
          ON CM.MEMBER_AC = M.METHOD_AC 
        """
    )

    clans = {}
    for row in cur:
        clan_acc = row[0]
        clan_name = row[1]
        clan_desc = row[2]
        database = row[3]
        member_acc = row[4]
        member_name = row[5]
        member_desc = row[6]
        seq_length = row[7]
        score = row[8]

        try:
            c = clans[clan_acc]
        except KeyError:
            c = clans[clan_acc] = {
                "accession": clan_acc,
                "name": clan_name,
                "description": clan_desc,
                "database": database,
                "members": []
            }
        finally:
            c["members"].append((member_acc, member_name, member_desc, score,
                                 seq_length))

    return clans


def iter_alignments(cur: oracledb.Cursor):
    # Fetching DOMAINS (LOB object) as a string
    cur.outputtypehandler = lob_as_str
    cur.execute(
        """
        SELECT QUERY_AC, TARGET_AC, EVALUE, DOMAINS
        FROM INTERPRO.CLAN_MATCH
        """
    )

    for query, target, evalue, json_domains in cur:
        domains = []
        for start, end in json.loads(json_domains):
            domains.append({
                "start": start,
                "end": end
            })

        yield query, target, evalue, domains


def export_clans(ipr_uri: str, pfam_uri: str, clans_file: str,
                 alignments_file: str, **kwargs):
    threshold = kwargs.get("threshold", 1e-2)

    logger.info("loading clans")
    con = oracledb.connect(ipr_uri)
    cur = con.cursor()
    clans = get_clans(cur)

    clan_links = {}
    entry2clan = {}
    for accession, clan in clans.items():
        clan_links[accession] = {}
        for member_acc, _, _, _, seq_length in clan["members"]:
            entry2clan[member_acc] = (accession, seq_length)

    logger.info("exporting alignments")
    with BasicStore(alignments_file, "w") as store:
        alignments = iter_alignments(cur)

        for i, (query, target, evalue, domains) in enumerate(alignments):
            if evalue > threshold:
                continue

            try:
                query_clan_acc, seq_length = entry2clan[query]
            except KeyError:
                continue

            try:
                target_clan_acc, _ = entry2clan[target]
            except KeyError:
                target_clan_acc = None

            store.write((query_clan_acc, query, target, target_clan_acc,
                         evalue, seq_length, json.dumps(domains)))

            if query_clan_acc == target_clan_acc:
                # Query and target from the same clan: update clan's links
                links = clan_links[query_clan_acc]

                if query > target:
                    query, target = target, query

                try:
                    targets = links[query]
                except KeyError:
                    links[query] = {target: evalue}
                else:
                    if target not in targets or evalue < targets[target]:
                        targets[target] = evalue

            if (i + 1) % 1e7 == 0:
                logger.info(f"{i + 1:>15,}")

        logger.info(f"{i + 1:>15,}")

    cur.close()
    con.close()

    logger.info("loading additional details for Pfam clans")
    pfam_clans = pfam.get_clans(pfam_uri)

    logger.info("finalizing")
    for clan_acc, clan in clans.items():
        nodes = []
        for member_acc, member_name, member_desc, score, _ in clan["members"]:
            nodes.append({
                "accession": member_acc,
                "short_name": member_name,
                "name": member_desc,
                "type": "entry",
                "score": score
            })

        links = []
        for query_acc, targets in clan_links[clan_acc].items():
            for target_acc, score in targets.items():
                links.append({
                    "source": query_acc,
                    "target": target_acc,
                    "score": score
                })

        clan["relationships"] = {
            "nodes": nodes,
            "links": links
        }

        if clan_acc in pfam_clans:
            # Replace `description`, add `authors` and `literature`
            clan.update(pfam_clans[clan_acc])

    with open(clans_file, "wb") as fh:
        pickle.dump(clans, fh)

    logger.info("complete")
