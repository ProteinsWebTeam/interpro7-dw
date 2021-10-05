import hashlib
from datetime import datetime

from interpro7dw import uniprot
from interpro7dw.interpro.utils import overlaps_pdb_chain
from interpro7dw.utils import logger
from interpro7dw.utils.store import SimpleStore, SimpleStoreSorter, Store
from interpro7dw.utils.store import copy_dict, dumpobj, loadobj

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


def dump_domain_organisation(url: str, proteins_src: str, matches_src: str,
                             domorgs_dst: str, **kwargs):
    """Calculates and exports the domain architectures/organisations of
    UniProt entries based on the Pfam matches.

    :param url: The Oracle connection string.
    :param proteins_src: The file containing protein information.
    :param matches_src: The file containing protein matches.
    :param domorgs_dst: The output domain organisation file.
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
    with SimpleStore(tempdir=tempdir) as tmp:
        with Store(proteins_src, "r") as st1, Store(matches_src, "r") as st2:
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

        size = tmp.size

        logger.info("exporting domain organisations")
        with Store(domorgs_dst, mode="w", keys=keys, tempdir=tempdir) as st:
            for i, (protein_acc, dom_str, dom_id, members) in enumerate(tmp):
                _, repr_acc = all_domains[dom_id]

                st.add(protein_acc, (dom_str, dom_id, members, repr_acc))

                if (i + 1) % 10000000 == 0:
                    logger.info(f"{i + 1:>15,}")

            logger.info(f"{i + 1:>15,}")

            size += st.size
            st.merge(workers, apply=st.get_first)

        logger.info(f"temporary files: {size / 1024 / 1024:.0f} MB")

    logger.info("done")


def dump_xrefs(url: str, proteins_src: str, matches_src: str,
               proteomes_src: str, domorgs_src: str, structures_src: str,
               xrefs_dst: str, **kwargs):
    """Export InterPro entries and member database signatures with proteins
    they match, and from this, assign proteomes, structures, and taxa to them.

    :param url: UniProt Oracle connection string.
    :param proteins_src: Store file of protein info.
    :param matches_src: Store file of protein matches.
    :param proteomes_src: Store file of protein-proteome mapping.
    :param domorgs_src: Store file of domain organisations.
    :param structures_src: File of PDBe structures.
    :param xrefs_dst: Output SimpleStore file.
    """
    buffersize = kwargs.get("buffersize", 1000000)
    tempdir = kwargs.get("tempdir")

    logger.info("loading data from UniProt database")
    protein2enzymes = uniprot.misc.get_swissprot2enzyme(url)
    protein2reactome = uniprot.misc.get_swissprot2reactome(url)

    # Create mapping protein -> structure -> chain -> locations
    logger.info("loading PDBe structures")
    protein2structures = {}
    for pdbe_id, entry in loadobj(structures_src).items():
        for protein_acc, chains in entry["proteins"].items():
            try:
                protein2structures[protein_acc][pdbe_id] = chains
            except KeyError:
                protein2structures[protein_acc] = {pdbe_id: chains}

    logger.info("iterating proteins")
    with SimpleStoreSorter(tempdir=tempdir) as stores:
        proteins = Store(proteins_src, "r")
        matches = Store(matches_src, "r")
        proteomes = Store(proteomes_src, "r")
        domorgs = Store(domorgs_src, "r")

        xrefs = {}
        num_xrefs = 0
        i = 0
        for i, (protein_acc, protein_matches) in enumerate(matches.items()):
            protein = proteins[protein_acc]
            protein_id = protein["identifier"]
            taxon_id = protein["taxid"]
            proteome_id = proteomes.get(protein_acc)
            structures = protein2structures.get(protein_acc, {})
            try:
                _, dom_id, dom_members, _ = domorgs[protein_acc]
            except KeyError:
                dom_id = None
                dom_members = []

            for entry_acc, locations in protein_matches.items():
                try:
                    entry_xrefs = xrefs[entry_acc]
                except KeyError:
                    entry_xrefs = xrefs[entry_acc] = {
                        "dom_orgs": set(),
                        "enzymes": set(),
                        "matches": 0,
                        "reactome": set(),
                        "proteins": set(),
                        "proteomes": set(),
                        "structures": set(),
                        "taxa": set(),
                    }

                entry_xrefs["matches"] += len(locations)
                entry_xrefs["proteins"].add((protein_acc, protein_id))
                entry_xrefs["taxa"].add(taxon_id)
                num_xrefs += 4

                if entry_acc in dom_members:
                    entry_xrefs["dom_orgs"].add(dom_id)
                    num_xrefs += 1

                if proteome_id:
                    entry_xrefs["proteomes"].add(proteome_id)
                    num_xrefs += 1

                for pdbe_id, chains in structures.items():
                    for chain_id, segments in chains.items():
                        if overlaps_pdb_chain(locations, segments):
                            entry_xrefs["structures"].add(pdbe_id)
                            num_xrefs += 1
                            break  # Skip other chains

                for ecno in protein2enzymes.get(protein_acc, []):
                    entry_xrefs["enzymes"].add(ecno)
                    num_xrefs += 1

                pathways = protein2reactome.get(protein_acc)
                for pathway_id, pathway_name in pathways:
                    entry_xrefs["reactome"].add((pathway_id, pathway_name))
                    num_xrefs += 2

            if num_xrefs >= buffersize:
                stores.dump(xrefs)
                xrefs.clear()
                num_xrefs = 0

            if (i + 1) % 10000000 == 0:
                logger.info(f"{i + 1:>15,}")

        logger.info(f"{i + 1:>15,}")

        proteins.close()
        matches.close()
        proteomes.close()

        stores.dump(xrefs)
        xrefs.clear()

        logger.info(f"temporary files: {stores.size / 1024 / 1024:.0f} MB")

        with SimpleStore(xrefs_dst) as store:
            for entry_acc, values in stores.merge():
                logger.info(entry_acc)
                xrefs = {}

                for entry_xrefs in values:
                    copy_dict(entry_xrefs, xrefs, concat_or_incr=True)

                store.add((entry_acc, xrefs))


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


def get_signatures(cur: cx_Oracle.Cursor) -> dict:
    cur.execute(
        """
        SELECT M.METHOD_AC, M.NAME, DB.DBSHORT, EVI.ABBREV,
               E2M.ENTRY_AC, E2M.NAME, E2M.ABBREV, E2M.PARENT_AC
        FROM INTERPRO.METHOD M
        INNER JOIN  INTERPRO.CV_DATABASE DB
          ON M.DBCODE = DB.DBCODE
        INNER JOIN  INTERPRO.IPRSCAN2DBCODE I2D
          ON M.DBCODE = I2D.DBCODE
        INNER JOIN INTERPRO.CV_EVIDENCE EVI
          ON I2D.EVIDENCE = EVI.CODE
        LEFT OUTER JOIN (
          SELECT E2M.METHOD_AC, E.ENTRY_AC, E.NAME, ET.ABBREV, E2E.PARENT_AC
          FROM INTERPRO.ENTRY E
          INNER JOIN INTERPRO.ENTRY2METHOD E2M
            ON E.ENTRY_AC = E2M.ENTRY_AC
          INNER JOIN INTERPRO.CV_ENTRY_TYPE ET
            ON E.ENTRY_TYPE = ET.CODE
          LEFT OUTER JOIN INTERPRO.ENTRY2ENTRY E2E
            ON E.ENTRY_AC = E2E.ENTRY_AC
          WHERE E.CHECKED = 'Y'
        ) E2M
          ON M.METHOD_AC = E2M.METHOD_AC
        """
    )
    signatures = {}
    for row in cur:
        signatures[row[0]] = {
            "accession": row[0],
            "name": row[1] or row[0],
            "database": row[2],
            "evidence": row[3],
            "interpro": {
                "id": row[4],
                "name": row[5],
                "type": row[6],
                "parent_id": row[7],
            } if row[4] else None
        }
    return signatures
