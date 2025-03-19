import json
import tarfile
from collections.abc import Callable
from pathlib import Path

import oracledb

from interpro7dw import uniprot
from ..oracle.entries import REPR_DOM_TYPES, REPR_DOM_DATABASES
from ..oracle.entries import REPR_FAM_TYPES, REPR_FAM_DATABASES
from interpro7dw.utils import logger


def package_data(ipr_uri: str, goa_uri: str, data_dir: str, ipr_version: str,
                 outdir: str):
    logger.info("Exporting JSON files")
    outdir = Path(outdir)

    pathways_file = outdir / "pathways.json"
    entry2pathways_file = outdir / "pathways.ipr.json"
    terms_file = outdir / "goterms.json"
    entry2terms_file = outdir / "goterms.ipr.json"
    entries_file = outdir / "entries.json"
    database_file = outdir / "database.json"

    con = oracledb.connect(ipr_uri)
    cur = con.cursor()

    _export_pathways(cur, pathways_file, entry2pathways_file)
    _export_go_terms(cur, goa_uri, terms_file, entry2terms_file)
    _export_entries(cur, entries_file, database_file, ipr_version)

    cur.execute(
        """
        SELECT LOWER(D.DBSHORT), V.VERSION
        FROM INTERPRO.IPRSCAN2DBCODE I2D
        INNER JOIN INTERPRO.CV_DATABASE D on I2D.DBCODE = D.DBCODE
        INNER JOIN INTERPRO.DB_VERSION V ON D.DBCODE = V.DBCODE
        """
    )
    versions = dict(cur.fetchall())
    cur.close()
    con.close()

    data_dir = Path(data_dir)
    output = Path(outdir) / "interproscan-data.tar.gz"
    prefix = f"interproscan-data-{ipr_version}/"
    logger.info("Building the InterPro archive")
    with tarfile.open(str(output), "w:gz") as tar:
        for file in [pathways_file, entry2pathways_file, terms_file,
                     entry2terms_file, entries_file]:
            tar.add(file, arcname=f"{prefix}xrefs/{file.name}")

        logger.info("Archiving AntiFam")
        pkg_antifam(data_dir, versions["antifam"], tar, prefix=prefix)

        logger.info("Archiving CATH (Gene3D + FunFam)")
        pkg_cath(data_dir, versions["cathgene3d"], tar, prefix=prefix)

        logger.info("Archiving CDD")
        pkg_cdd(data_dir, versions["cdd"], tar, prefix=prefix)

        logger.info("Archiving HAMAP")
        pkg_hamap(data_dir, versions["hamap"], tar, prefix=prefix)

        logger.info("Archiving NCBIfam")
        pkg_ncbifam(data_dir, versions["ncbifam"], tar, prefix=prefix)

        logger.info("Archiving PANTHER")
        pkg_panther(data_dir, versions["panther"], tar, prefix=prefix)

        logger.info("Archiving Pfam")
        pkg_pfam(data_dir, versions["pfam"], tar, prefix=prefix)

        logger.info("Archiving PIRSF")
        pkg_pirsf(data_dir, versions["pirsf"], tar, prefix=prefix)

        logger.info("Archiving PIRSR")
        pkg_pirsr(data_dir, versions["pirsr"], tar, prefix=prefix)

        logger.info("Archiving PRINTS")
        pkg_prints(data_dir, versions["prints"], tar, prefix=prefix)

        logger.info("Archiving PROSITE (Patterns + Profiles)")
        pkg_prosite(data_dir, versions["prosite"], tar, prefix=prefix)

        logger.info("Archiving SFLD")
        pkg_sfld(data_dir, versions["sfld"], tar, prefix=prefix)

        logger.info("Archiving SMART")
        pkg_smart(data_dir, versions["smart"], tar, prefix=prefix)

        logger.info("Archiving SUPERFAMILY")
        pkg_superfamily(data_dir, versions["ssf"], tar, prefix=prefix)

    logger.info("Building InterPro Member archives")
    logger.info("Archiving AntiFam")
    build_member_archive("AntiFam", versions['antifam'], outdir, data_dir, pkg_antifam)

    logger.info("Archiving CATH (Gene3D + FunFam)")
    build_member_archive("Cath", versions['cathgene3d'], outdir, data_dir, pkg_cath)

    logger.info("Archiving CDD")
    build_member_archive("CDD", versions['cdd'], outdir, data_dir, pkg_cdd)

    logger.info("Archiving HAMAP")
    build_member_archive("HAMAP", versions['hamap'], outdir, data_dir, pkg_hamap)

    logger.info("Archiving NCBIfam")
    build_member_archive("NCBIfam", versions['ncbifam'], outdir, data_dir, pkg_ncbifam)

    logger.info("Archiving PANTHER")
    build_member_archive("PANTHER", versions['panther'], outdir, data_dir, pkg_panther)

    logger.info("Archiving Pfam")
    build_member_archive("Pfam", versions['pfam'], outdir, data_dir, pkg_pfam)

    logger.info("Archiving PIRSF")
    build_member_archive("PIRSF", versions['pirsf'], outdir, data_dir, pkg_pirsf)

    logger.info("Archiving PIRSR")
    build_member_archive("PIRSR", versions['pirsr'], outdir, data_dir, pkg_pirsr)

    logger.info("Archiving PRINTS")
    build_member_archive("PRINTS", versions['prints'], outdir, data_dir, pkg_prints)

    logger.info("Archiving PROSITE (Patterns + Profiles)")
    build_member_archive("PROSITE", versions['prosite'], outdir, data_dir, pkg_prosite)

    logger.info("Archiving SFLD")
    build_member_archive("SFLD", versions['sfld'], outdir, data_dir, pkg_sfld)

    logger.info("Archiving SMART")
    build_member_archive("SMART", versions['smart'], outdir, data_dir, pkg_smart)

    logger.info("Archiving SUPERFAMILY")
    build_member_archive("SUPERFAMILY", versions['ssf'], outdir, data_dir, pkg_superfamily)

    logger.info("Done")


def _export_pathways(
        cur: oracledb.Cursor,
        pathways_file: Path,
        entry2pathways_file: Path
):
    cur.execute(
        """
        SELECT ENTRY_AC, DBCODE, AC, NAME
        FROM INTERPRO.ENTRY2PATHWAY
        """
    )

    pathways = {}
    interpro2pathways = {}
    for entry_acc, dbcode, pathway_id, pathway_name in cur:

        try:
            interpro2pathways[entry_acc].append(pathway_id)
        except KeyError:
            interpro2pathways[entry_acc] = [pathway_id]

        pathways[pathway_id] = [dbcode, pathway_name]

    with pathways_file.open("wt") as fh:
        json.dump(pathways, fh)

    with entry2pathways_file.open("wt") as fh:
        json.dump(interpro2pathways, fh)


def _export_go_terms(
        cur: oracledb.Cursor,
        goa_uri: str,
        terms_file: Path,
        entry2terms_file: Path
):
    version = uniprot.goa.get_timestamp(goa_uri).strftime("%Y-%m-%d")
    terms = {}
    for go_id, (name, aspect, _, _) in uniprot.goa.get_terms(goa_uri).items():
        terms[go_id] = [name, aspect]

    cur.execute(
        """
        SELECT ENTRY_AC, GO_ID
        FROM INTERPRO.INTERPRO2GO
        WHERE ENTRY_AC IN (
            SELECT ENTRY_AC
            FROM INTERPRO.ENTRY
            WHERE CHECKED = 'Y'
        )
        """
    )

    interpro2go = {}
    for entry_acc, go_id in cur:
        if go_id not in terms:
            logger.error(f"{entry_acc}: term {go_id} not found")
            continue
        elif entry_acc in interpro2go:
            interpro2go[entry_acc].append(go_id)
        else:
            interpro2go[entry_acc] = [go_id]

    with terms_file.open("wt") as fh:
        json.dump({
            "version": version,
            "terms": terms
        }, fh)

    with entry2terms_file.open("wt") as fh:
        json.dump(interpro2go, fh)


def _export_entries(cur: oracledb.Cursor, entries_file: Path, database_file: Path, ipr_version: str):
    cur.execute("SELECT CODE, ABBREV FROM INTERPRO.CV_ENTRY_TYPE")
    types = dict(cur.fetchall())

    cur.execute(
        """
        SELECT D.DBCODE, LOWER(D.DBSHORT), D.DBNAME, V.VERSION
        FROM INTERPRO.CV_DATABASE D
        INNER JOIN INTERPRO.DB_VERSION V ON D.DBCODE = V.DBCODE
        WHERE D.DBCODE IN (
            SELECT DBCODE FROM INTERPRO.IPRSCAN2DBCODE
        ) 
           OR D.DBCODE = 'I'
        """
    )
    databases = {row[0]: row[1:] for row in cur.fetchall()}

    cur.execute(
        r"""
        SELECT E.ENTRY_AC, E.SHORT_NAME, E.NAME, E.ENTRY_TYPE, 'I', NULL
        FROM INTERPRO.ENTRY E
        WHERE E.CHECKED = 'Y'
        UNION ALL
        SELECT M.METHOD_AC, M.NAME, M.DESCRIPTION, M.SIG_TYPE, M.DBCODE, 
               EM.ENTRY_AC
        FROM INTERPRO.METHOD M
        LEFT OUTER JOIN (
            SELECT E.ENTRY_AC, EM.METHOD_AC
            FROM INTERPRO.ENTRY E
            INNER JOIN INTERPRO.ENTRY2METHOD EM
                ON E.ENTRY_AC = EM.ENTRY_AC
            WHERE E.CHECKED = 'Y'
        ) EM ON M.METHOD_AC = EM.METHOD_AC
        UNION ALL
        SELECT FM.METHOD_AC, FM.NAME, FM.DESCRIPTION, 'G', FM.DBCODE, NULL
        FROM INTERPRO.FEATURE_METHOD FM
        WHERE FM.DBCODE IN ('a', 'f', 'g', 'j', 'n', 'q', 's', 'v', 'x')
        """
    )

    entries = {}
    for row in cur.fetchall():
        dbshort, dbname, dbversion = databases[row[4]]

        if (types[row[3]].lower() in REPR_DOM_TYPES and
                dbshort in REPR_DOM_DATABASES):
            repr_type = "domain"
            repr_index = REPR_DOM_DATABASES.index(dbshort)
        elif (types[row[3]].lower() in REPR_FAM_TYPES and
                dbshort in REPR_FAM_DATABASES):
            repr_type = "family"
            repr_index = REPR_FAM_DATABASES.index(dbshort)
        else:
            repr_type = None
            repr_index = 0

        entries[row[0]] = {
            "name": row[1],
            "description": row[2],
            "type": types[row[3]],
            "integrated": row[5],
            "representative": {
                "type": repr_type,
                "index": repr_index
            },
            "database": dbname
        }

    with entries_file.open("wt") as fh:
        json.dump(entries, fh)

    databases = {n: v for _, n, v in databases.values()}
    with database_file.open("wt") as fh:
        json.dump(databases, fh)


def build_member_archive(member: str, version: str, outdir: Path, data_dir: Path,
                        pkg_func: Callable[[Path, str, tarfile.TarFile], None]):
    member_output = outdir / f"{member}-{version}.tar.gz"
    with tarfile.open(str(member_output), "w:gz") as member_tar:
        pkg_func(data_dir, version, member_tar)


def pkg_antifam(root: Path, version: str, tar: tarfile.TarFile,
                prefix: str = ""):
    path = root / "antifam" / version / "AntiFam.hmm"
    tar.add(path, arcname=f"{prefix}antifam/AntiFam.hmm")


def pkg_cath(root: Path, version: str, tar: tarfile.TarFile, prefix: str = ""):
    base = root / "cath-gene3d" / version

    path = base / "gene3d_main.hmm"
    tar.add(path, arcname=f"{prefix}cath/gene3d/{path.name}")

    path = base / "discontinuous" / "discontinuous_regs.pkl.py3"
    tar.add(path, arcname=f"{prefix}cath/gene3d/{path.name}")

    path = base / "model_to_family_map.tsv"
    tar.add(path, arcname=f"{prefix}cath/gene3d/{path.name}")

    path = base / "funfam" / "models"
    for child in path.iterdir():
        tar.add(child, arcname=f"{prefix}cath/funfam/{child.name}")


def pkg_cdd(root: Path, version: str, tar: tarfile.TarFile, prefix: str = ""):
    path = root / "cdd" / version / "data"
    tar.add(path, arcname=f"{prefix}cdd/data")

    path = root / "cdd" / version / "db"
    tar.add(path, arcname=f"{prefix}cdd/db")


def pkg_hamap(root: Path, version: str, tar: tarfile.TarFile, prefix: str = ""):
    path = root / "hamap" / version / "hamap.prf"
    tar.add(path, arcname=f"{prefix}hamap/{path.name}")

    path = root / "hamap" / version / "hamap.hmm.lib"
    tar.add(path, arcname=f"{prefix}hamap/{path.name}")

    path = root / "hamap" / version / "profiles"
    tar.add(path, arcname=f"{prefix}hamap/{path.name}")


def pkg_ncbifam(root: Path, version: str, tar: tarfile.TarFile,
                prefix: str = ""):
    path = root / "ncbifam" / version / "ncbifam.hmm"
    tar.add(path, arcname=f"{prefix}ncbifam/{path.name}")


def pkg_panther(root: Path, version: str, tar: tarfile.TarFile,
                prefix: str = ""):
    path = root / "panther" / version / "famhmm"
    tar.add(path, arcname=f"{prefix}panther/{path.name}")

    path = root / "panther" / version / "PAINT_Annotations"
    for child in path.iterdir():
        if child.suffix == ".json":
            name = f"{prefix}panther/{child.parent.name}/{child.name}"
            tar.add(child, arcname=name)

    path = root / "panther" / version / "Tree_MSF"
    tar.add(path, arcname=f"{prefix}panther/{path.name}")


def pkg_pfam(root: Path, version: str, tar: tarfile.TarFile, prefix: str = ""):
    path = root / "pfam" / version / "pfam_a.hmm"
    tar.add(path, arcname=f"{prefix}pfam/{path.name}")

    path = root / "pfam" / version / "pfam_clans"
    tar.add(path, arcname=f"{prefix}pfam/{path.name}")

    path = root / "pfam" / version / "pfam_a.dat"
    tar.add(path, arcname=f"{prefix}pfam/{path.name}")

    path = root / "pfam" / version / "pfam_a.seed"
    tar.add(path, arcname=f"{prefix}pfam/{path.name}")


def pkg_pirsf(root: Path, version: str, tar: tarfile.TarFile, prefix: str = ""):
    path = root / "pirsf" / version / "pirsf.dat"
    tar.add(path, arcname=f"{prefix}pirsf/{path.name}")

    path = root / "pirsf" / version / "sf_hmm_all"
    tar.add(path, arcname=f"{prefix}pirsf/{path.name}")


def pkg_pirsr(root: Path, version: str, tar: tarfile.TarFile, prefix: str = ""):
    path = root / "pirsr" / version / "sr_hmm_all"
    tar.add(path, arcname=f"{prefix}pirsr/{path.name}")

    path = root / "pirsr" / version / "sr_uru.json"
    tar.add(path, arcname=f"{prefix}pirsr/{path.name}")


def pkg_prints(root: Path, version: str, tar: tarfile.TarFile,
               prefix: str = ""):
    members = [
        ("FingerPRINTShierarchy21Feb2012", "FingerPRINTShierarchy.db"),
        ("prints42_0.pval_blos62", "prints.pval"),
    ]

    for src, dst in members:
        path = root / "prints" / version / src
        tar.add(path, arcname=f"{prefix}prints/{dst}")


def pkg_prosite(root: Path, version: str, tar: tarfile.TarFile,
                prefix: str = ""):
    path = root / "prosite" / version / "evaluator.dat"
    tar.add(path, arcname=f"{prefix}prosite/{path.name}")

    path = root / "prosite" / version / "prosite_patterns.dat"
    tar.add(path, arcname=f"{prefix}prosite/{path.name}")

    path = root / "prosite" / version / "prosite_profiles"
    tar.add(path, arcname=f"{prefix}prosite/{path.name}")

    path = root / "prosite" / version / "skip_flagged_profiles.txt"
    tar.add(path, arcname=f"{prefix}prosite/{path.name}")


def pkg_sfld(root: Path, version: str, tar: tarfile.TarFile, prefix: str = ""):
    members = [
        "sfld.hmm",
        "sfld_sites.annot",
        "sfld_hierarchy.txt",
    ]

    for member in members:
        path = root / "sfld" / version / member
        tar.add(path, arcname=f"{prefix}sfld/{path.name}")


def pkg_smart(root: Path, version: str, tar: tarfile.TarFile, prefix: str = ""):
    members = [
        "smart.HMMs",
        "smart.HMMs.bin",
    ]

    for member in members:
        path = root / "smart" / version / member
        tar.add(path, arcname=f"{prefix}smart/{path.name}")


def pkg_superfamily(root: Path, version: str, tar: tarfile.TarFile,
                    prefix: str = ""):
    members = [
        "hmmlib_1.75",
        "hmmlib_1.75.h3f",
        "hmmlib_1.75.h3i",
        "hmmlib_1.75.h3m",
        "hmmlib_1.75.h3p",
        "self_hits.tab",
        "dir.cla.scop.txt_1.75",
        "dir.des.scop.txt_1.75",
        "model.tab",
        "pdbj95d",
        "LICENSE",
    ]

    for member in members:
        path = root / "superfamily" / version / member
        tar.add(path, arcname=f"{prefix}superfamily/{path.name}")
