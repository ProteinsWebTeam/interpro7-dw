import hashlib
from datetime import datetime

from interpro7dw.utils import logger
from interpro7dw.utils.store import SimpleStore, Store, dumpobj

import cx_Oracle


ENTRY_DATABASES = [
    'B',    # SFLD
    'F',    # PRINTS
    'H',    # Pfam
    'I',    # InterPro
    'J',    # CDD
    'M',    # PROSITE profiles
    'N',    # TIGRFAMs
    'P',    # PROSITE patterns
    'Q',    # HAMAP
    'R',    # SMART
    'U',    # PIRSF
    'V',    # PANTHER
    'X',    # CATH-Gene3D
    'Y',    # SUPERFAMILY
]
FEATURE_DATABASES = [
    'g',    # MobiDB Lite
    'j',    # Phobius
    'n',    # Signal Euk
    'q',    # TMHMM
    's',    # SignalP Gram positive
    'v',    # SignalP Gram negative
    'x',    # COILS
]
SEQUENCE_DATABASES = [
    'S',    # Swiss-Prot
    'T',    # TrEMBL
    'u',    # UniProtKB
]


def dump_databases(url: str, version: str, date: str, file: str,
                   update: bool = False):
    """Exports information on databases/data sources used in InterPro.

    :param url: The Oracle connection string.
    :param version: The version of the upcoming InterPro release.
    :param date: The date of the upcoming InterPro release (YYYY-MM-DD).
    :param file: The output file.
    :param update: If True, update the production table.
    """
    con = cx_Oracle.connect(url)
    cur = con.cursor()

    cur.execute("SELECT COUNT(*) FROM INTERPRO.ENTRY WHERE CHECKED = 'Y'")
    num_interpro_entries, = cur.fetchone()

    cur.execute("SELECT VERSION FROM INTERPRO.DB_VERSION WHERE DBCODE = 'I'")
    prod_version, = cur.fetchone()

    if prod_version == version:
        # DB_VERSION is already up-to-date
        use_db_version = True
    elif update:
        # DB_VERSION is outdated, but will be up-to-date
        use_db_version = True
        cur.execute(
            """
            UPDATE INTERPRO.DB_VERSION
            SET VERSION = :1,
                FILE_DATE = :2,
                ENTRY_COUNT = :3
            WHERE DBCODE = 'I'
            """, (version, datetime.strptime(date, "%Y-%m-%d"),
                  num_interpro_entries)
        )
        con.commit()
    else:
        # DB_VERSION is outdated and will stay outdated
        # This run is a test done on the production database (SRSLY?!!111)
        use_db_version = False

    """
    Using RN=2 to join with the second most recent action in DB_VERSION_AUDIT
    (the most recent is the same record as in DB_VERSION)
    """
    cur.execute(
        """
        SELECT
          DB.DBCODE, LOWER(DB.DBSHORT), DB.DBSHORT, DB.DBNAME, 
          DB.DESCRIPTION, V.VERSION, V.FILE_DATE, V.ENTRY_COUNT, VA.VERSION, 
          VA.FILE_DATE
        FROM INTERPRO.CV_DATABASE DB
        LEFT OUTER JOIN INTERPRO.DB_VERSION V ON DB.DBCODE = V.DBCODE
        LEFT OUTER JOIN (
          SELECT
            DBCODE, VERSION, FILE_DATE,
            ROW_NUMBER() OVER (
              PARTITION BY DBCODE ORDER BY TIMESTAMP DESC
            ) RN
          FROM INTERPRO.DB_VERSION_AUDIT
          WHERE ACTION = 'U'
        ) VA ON DB.DBCODE = VA.DBCODE AND VA.RN = 2
        """
    )

    databases = []
    for rec in cur:
        code = rec[0]
        identifier = rec[1]
        short_name = rec[2]
        name = rec[3]
        description = rec[4]
        release_version = rec[5]
        release_date = rec[6]
        num_entries = rec[7]
        prev_release_version = rec[8]
        prev_release_date = rec[9]

        if code in ENTRY_DATABASES:
            db_type = "entry"
        elif code in FEATURE_DATABASES:
            db_type = "feature"
        elif code in SEQUENCE_DATABASES:
            if code == 'S':
                identifier = "reviewed"
            elif code == 'T':
                identifier = "unreviewed"

            db_type = "protein"
        else:
            db_type = "other"

        if code == 'I' and not use_db_version:
            # DB_VERSION is outdated:
            # it contains info for the live release (soon to be 'previous')
            num_entries = num_interpro_entries
            prev_release_version = release_version
            prev_release_date = release_date
            release_version = version
            release_date = datetime.strptime(date, "%Y-%m-%d")

        databases.append((
            identifier,
            name,
            short_name,
            description,
            db_type,
            num_entries,
            release_version,
            release_date,
            prev_release_version,
            prev_release_date
        ))

    cur.close()
    con.close()

    dumpobj(databases, file)


def dump_domain_organisation(url: str, proteins_src: str, matches_str,
                             domorg_dst: str, **kwargs):
    """Calculates and exports the domain architectures/organisations of
    UniProt entries based on the Pfam matches.

    :param url: The Oracle connection string.
    :param proteins_src: The file containing protein information.
    :param matches_str: The file containing protein matches.
    :param domorg_dst: The output domain organisation file.
    """
    tempdir = kwargs.get("tempdir")
    workers = kwargs.get("workers", 1)

    # Loads Pfam signatures, and the InterPro entries they are integrated in
    logger.info("loading Pfam signatures")
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    cur = cur.execute(
        """
        SELECT M.METHOD_AC, EM.ENTRY_AC
        FROM INTERPRO.METHOD M
        LEFT OUTER JOIN (
            SELECT EM.METHOD_AC, E.ENTRY_AC
            FROM INTERPRO.ENTRY2METHOD EM
            INNER JOIN INTERPRO.ENTRY E ON EM.ENTRY_AC = E.ENTRY_AC
            WHERE E.CHECKED = 'Y'
        ) EM ON M.METHOD_AC = EM.METHOD_AC
        WHERE M.DBCODE = 'H'  -- Pfam
        """
    )

    pfam_signatures = {}
    for pfam_acc, interpro_acc in cur:
        pfam_signatures[pfam_acc] = interpro_acc

    cur.close()
    con.close()

    logger.info("iterating protein matches")
    all_domains = {}
    with SimpleStore(domorg_dst) as tmp:
        with Store(proteins_src, "r") as st1, Store(matches_str, "r") as st2:
            keys = st1.file_keys

            for i, (protein_acc, matches) in enumerate(st2.items()):
                locations = []
                for entry_acc in matches:
                    try:
                        interpro_acc = pfam_signatures[entry_acc]
                    except KeyError:
                        # Not a Pfam match
                        continue

                    for loc in matches[entry_acc]:
                        locations.append({
                            "pfam": entry_acc,
                            "interpro": interpro_acc,
                            # We do not consider fragmented locations
                            "start": loc["fragments"][0]["start"],
                            "end": max(f["end"] for f in loc["fragments"])
                        })

                if (i + 1) % 10000000 == 0:
                    logger.info(f"{i + 1:>15,}")

                if not locations:
                    continue  # No Pfam matches: no domain organisation

                domains = []
                members = set()
                for loc in sorted(locations, key=lambda l: (l["start"],
                                                            l["end"])):
                    pfam_acc = loc["pfam"]
                    interpro_acc = loc["interpro"]

                    if interpro_acc:
                        domains.append(f"{pfam_acc}:{interpro_acc}")
                        members.add(interpro_acc)
                    else:
                        domains.append(pfam_acc)

                    members.add(pfam_acc)

                dom_str = '-'.join(domains)
                dom_id = hashlib.sha1(dom_str.encode("utf-8")).hexdigest()
                tmp.add((protein_acc, dom_str, dom_id, members))

                # string (YYYY-MM-DD) which is enough to compare dates
                date = st1[protein_acc]["date"]

                # Selects the oldest protein to represent
                # this domain organisation.
                try:
                    other_date, _ = all_domains[dom_id]
                except KeyError:
                    all_domains[dom_id] = (date, protein_acc)
                else:
                    if date < other_date:
                        all_domains[dom_id] = (date, protein_acc)

            logger.info(f"{i + 1:>15,}")

        logger.info("exporting domain organisations")
        with Store(domorg_dst, mode="w", keys=keys, tempdir=tempdir) as st:
            for i, (protein_acc, dom_str, dom_id, members) in enumerate(tmp):
                _, repr_acc = all_domains[dom_id]

                st.add(protein_acc, (dom_str, dom_id, members, repr_acc))

                if (i + 1) % 10000000 == 0:
                    logger.info(f"{i + 1:>15,}")

            logger.info(f"{i + 1:>15,}")

            st.merge(workers, apply=st.get_first)
            logger.info(f"temporary files: {st.size / 1024 / 1024:.0f} MB")

    logger.info("done")


def dump_similar_entries(url: str, matches_src: str, relationships_dst: str,
                         min_similarity: float = 0.75):
    con = cx_Oracle.connect(url)
    cur = con.cursor()
    cur.execute(
        """
        SELECT E.ENTRY_AC, LOWER(ET.ABBREV)
        FROM INTERPRO.ENTRY E
        INNER JOIN INTERPRO.CV_ENTRY_TYPE ET ON E.ENTRY_TYPE = ET.CODE
        WHERE E.CHECKED = 'Y'
        """
    )
    entries = dict(cur.fetchall())
    cur.close()
    con.close()

    num_proteins = {}  # number of proteins matched by entry
    num_overlaps = {}  # number of proteins where two entries overlap >= 50%

    with Store(matches_src, "r") as store:
        for i, (protein_acc, matches) in store.items():
            entries = {}
            for entry_acc, locations in matches.items():
                if entry_acc in entries:
                    entries[entry_acc] = []

                    for loc in locations:
                        # InterPro locations have one fragment only
                        entries[entry_acc].append((
                            loc["fragments"][0]["start"],
                            loc["fragments"][0]["end"],
                        ))

            # Evaluate how entries overlap
            for entry_acc, locations in entries.items():
                try:
                    num_proteins[entry_acc] += 1
                except KeyError:
                    num_proteins[entry_acc] = 1

                for other_acc, other_locations in entries.items():
                    if other_acc >= entry_acc:
                        continue

                    try:
                        entry_overlaps = num_overlaps[entry_acc]
                    except KeyError:
                        entry_overlaps = num_overlaps[entry_acc] = {}

                    try:
                        overlaps = entry_overlaps[other_acc]
                    except KeyError:
                        overlaps = entry_overlaps[other_acc] = [0, 0]

                    flag = 0
                    for start1, end1 in locations:
                        length1 = end1 - start1 + 1

                        for start2, end2 in other_locations:
                            length2 = end2 - start2 + 1
                            overlap = min(end1, end2) - max(start1, start2) + 1

                            if not flag & 1 and overlap >= length1 * 0.5:
                                # 1st time fragments overlap >= 50% of entry 1
                                flag |= 1
                                overlaps[0] += 1

                            if not flag & 2 and overlap >= length2 * 0.5:
                                # 1st time fragments overlap >= 50% of entry 2
                                flag |= 2
                                overlaps[1] += 1

                        if flag == 3:
                            # Both cases already happened
                            break

            if (i + 1) % 10000000 == 0:
                logger.info(f"{i + 1:>15,}")

    logger.info(f"{i + 1:>15,}")

    supfam = "homologous_superfamily"
    types = (supfam, "domain", "family", "repeat")
    overlapping_entries = []
    for entry_acc, overlaps in num_overlaps.items():
        entry_cnt = num_proteins[entry_acc]

        for other_acc, (cnt1, cnt2) in overlaps.items():
            other_cnt = num_proteins[other_acc]

            # Independent coefficients
            coef1 = cnt1 / (entry_cnt + other_cnt - cnt1)
            coef2 = cnt2 / (entry_cnt + other_cnt - cnt2)

            # Final coefficient (average of independent coefficients)
            coef = (coef1 + coef2) * 0.5

            # Containment indices
            cont1 = cnt1 / entry_cnt
            cont2 = cnt2 / other_cnt

            if all(e < min_similarity for e in (coef, cont1, cont2)):
                continue

            # Entries are deemed similar
            type1 = entries[entry_acc]
            type2 = entries[other_acc]
            if ((type1 == supfam and type2 in types)
                    or (type2 == supfam and type1 in types)):
                overlapping_entries.append((entry_acc, other_acc))

    dumpobj(overlapping_entries, relationships_dst)
