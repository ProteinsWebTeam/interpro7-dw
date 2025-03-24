import json
import tarfile
from collections.abc import Callable
from pathlib import Path

import oracledb

from interpro7dw import uniprot
from ..oracle.entries import REPR_DOM_TYPES, REPR_DOM_DATABASES
from ..oracle.entries import REPR_FAM_TYPES, REPR_FAM_DATABASES
from interpro7dw.utils import logger


def package_data(
    ipr_uri: str, goa_uri: str, data_dir: str, ipr_version: str, outdir: str
):
    """
    Package InterProScan data files into gzip-compressed tar archives
    :param ipr_uri: InterPro Oracle connection string
    :param goa_uri: GOA Oracle connection string
    :param data_dir: Path to root directory where member database data files are stored
    :param ipr_version: InterPro release version
    :param outdir: Output directory for archives
    """
    logger.info("Exporting JSON files")
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    pathways_file = outdir / "pathways.json"
    entry2pathways_file = outdir / "pathways.ipr.json"
    go_terms_file = outdir / "goterms.json"
    entry2go_terms_file = outdir / "goterms.ipr.json"
    entries_file = outdir / "entries.json"
    databases_file = outdir / "database.json"

    con = oracledb.connect(ipr_uri)
    cur = con.cursor()

    _export_pathways(cur, pathways_file, entry2pathways_file)
    _export_go_terms(cur, goa_uri, go_terms_file, entry2go_terms_file)
    _export_entries(cur, entries_file, databases_file)

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

    logger.info("Creating InterPro archive")
    output = Path(outdir) / "interpro" / f"{ipr_version}.tar.gz"
    with tarfile.open(str(output), "w:gz") as tar:
        for file in [
            pathways_file,
            entry2pathways_file,
            go_terms_file,
            entry2go_terms_file,
            entries_file,
            databases_file,
        ]:
            tar.add(file, arcname=f"interpro/{ipr_version}/{file.name}")

    data_dir = Path(data_dir)

    logger.info("Creating AntiFam archive")
    build_member_archive("antifam", versions["antifam"], data_dir, outdir, pkg_antifam)

    logger.info("Creating CATH (Gene3D + FunFam) archive")
    build_member_archive("cath", versions["cathgene3d"], data_dir, outdir, pkg_cath)

    logger.info("Creating CDD archive")
    build_member_archive("cdd", versions["cdd"], data_dir, outdir, pkg_cdd)

    logger.info("Creating HAMAP archive")
    build_member_archive("hamap", versions["hamap"], data_dir, outdir, pkg_hamap)

    logger.info("Creating NCBIfam archive")
    build_member_archive("ncbifam", versions["ncbifam"], data_dir, outdir, pkg_ncbifam)

    logger.info("Creating PANTHER archive")
    build_member_archive("panther", versions["panther"], data_dir, outdir, pkg_panther)

    logger.info("Creating Pfam archive")
    build_member_archive("pfam", versions["pfam"], data_dir, outdir, pkg_pfam)

    logger.info("Creating PIRSF archive")
    build_member_archive("pirsf", versions["pirsf"], data_dir, outdir, pkg_pirsf)

    logger.info("Creating PIRSR archive")
    build_member_archive("pirsr", versions["pirsr"], data_dir, outdir, pkg_pirsr)

    logger.info("Creating PRINTS archive")
    build_member_archive("prints", versions["prints"], data_dir, outdir, pkg_prints)

    logger.info("Creating PROSITE (Patterns + Profiles) archive")
    build_member_archive("prosite", versions["prosite"], data_dir, outdir, pkg_prosite)

    logger.info("Creating SFLD archive")
    build_member_archive("sfld", versions["sfld"], data_dir, outdir, pkg_sfld)

    logger.info("Creating SMART archive")
    build_member_archive("smart", versions["smart"], data_dir, outdir, pkg_smart)

    logger.info("Creating SUPERFAMILY archive")
    build_member_archive(
        "superfamily", versions["ssf"], data_dir, outdir, pkg_superfamily
    )

    logger.info("Done")


def _export_pathways(
    cur: oracledb.Cursor, pathways_file: Path, entry2pathways_file: Path
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
    cur: oracledb.Cursor, goa_uri: str, terms_file: Path, entry2terms_file: Path
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
        json.dump({"version": version, "terms": terms}, fh)

    with entry2terms_file.open("wt") as fh:
        json.dump(interpro2go, fh)


def _export_entries(cur: oracledb.Cursor, entries_file: Path, databases_file: Path):
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

        if types[row[3]].lower() in REPR_DOM_TYPES and dbshort in REPR_DOM_DATABASES:
            repr_type = "domain"
            repr_index = REPR_DOM_DATABASES.index(dbshort)
        elif types[row[3]].lower() in REPR_FAM_TYPES and dbshort in REPR_FAM_DATABASES:
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
            "representative": {"type": repr_type, "index": repr_index},
            "database": dbname,
        }

    with entries_file.open("wt") as fh:
        json.dump(entries, fh)

    databases = {n: v for _, n, v in databases.values()}
    with databases_file.open("wt") as fh:
        json.dump(databases, fh)


def build_member_archive(
    member: str,
    version: str,
    indir: Path,
    outdir: Path,
    fn: Callable[[Path, str, tarfile.TarFile], None],
):
    output = outdir / member / f"{member}-{version}.tar.gz"
    with tarfile.open(str(output), "w:gz") as tar:
        fn(indir, version, tar)


def pkg_antifam(root: Path, version: str, tar: tarfile.TarFile):
    path = root / "antifam" / version / "AntiFam.hmm"
    tar.add(path, arcname=f"antifam/{version}/AntiFam.hmm")


def pkg_cath(root: Path, version: str, tar: tarfile.TarFile):
    base = root / "cath-gene3d" / version

    path = base / "gene3d_main.hmm"
    tar.add(path, arcname=f"cath/{version}/gene3d/{path.name}")

    path = base / "discontinuous" / "discontinuous_regs.pkl.py3"
    tar.add(path, arcname=f"cath/{version}/gene3d/{path.name}")

    path = base / "model_to_family_map.tsv"
    tar.add(path, arcname=f"cath/{version}/gene3d/{path.name}")

    path = base / "funfam" / "models"
    for child in path.iterdir():
        tar.add(child, arcname=f"cath/{version}/funfam/{child.name}")


def pkg_cdd(root: Path, version: str, tar: tarfile.TarFile):
    path = root / "cdd" / version / "data"
    tar.add(path, arcname=f"cdd/{version}/data")

    path = root / "cdd" / version / "db"
    tar.add(path, arcname=f"cdd/{version}/db")


def pkg_hamap(root: Path, version: str, tar: tarfile.TarFile):
    path = root / "hamap" / version / "hamap.prf"
    tar.add(path, arcname=f"hamap/{version}/{path.name}")

    path = root / "hamap" / version / "hamap.hmm.lib"
    tar.add(path, arcname=f"hamap/{version}/{path.name}")

    path = root / "hamap" / version / "profiles"
    tar.add(path, arcname=f"hamap/{version}/{path.name}")


def pkg_ncbifam(root: Path, version: str, tar: tarfile.TarFile):
    path = root / "ncbifam" / version / "ncbifam.hmm"
    tar.add(path, arcname=f"ncbifam/{version}/{path.name}")


def pkg_panther(root: Path, version: str, tar: tarfile.TarFile):
    path = root / "panther" / version / "famhmm"
    tar.add(path, arcname=f"panther/{version}/{path.name}")

    path = root / "panther" / version / "PAINT_Annotations"
    for child in path.iterdir():
        if child.suffix == ".json":
            name = f"panther/{version}/{child.parent.name}/{child.name}"
            tar.add(child, arcname=name)

    path = root / "panther" / version / "Tree_MSF"
    tar.add(path, arcname=f"panther/{version}/{path.name}")


def pkg_pfam(root: Path, version: str, tar: tarfile.TarFile):
    path = root / "pfam" / version / "pfam_a.hmm"
    tar.add(path, arcname=f"pfam/{version}/{path.name}")

    path = root / "pfam" / version / "pfam_clans"
    tar.add(path, arcname=f"pfam/{version}/{path.name}")

    path = root / "pfam" / version / "pfam_a.dat"
    tar.add(path, arcname=f"pfam/{version}/{path.name}")

    path = root / "pfam" / version / "pfam_a.seed"
    tar.add(path, arcname=f"pfam/{version}/{path.name}")


def pkg_pirsf(root: Path, version: str, tar: tarfile.TarFile):
    path = root / "pirsf" / version / "pirsf.dat"
    tar.add(path, arcname=f"pirsf/{version}/{path.name}")

    path = root / "pirsf" / version / "sf_hmm_all"
    tar.add(path, arcname=f"pirsf/{version}/{path.name}")


def pkg_pirsr(root: Path, version: str, tar: tarfile.TarFile):
    path = root / "pirsr" / version / "sr_hmm_all"
    tar.add(path, arcname=f"pirsr/{version}/{path.name}")

    path = root / "pirsr" / version / "sr_uru.json"
    tar.add(path, arcname=f"pirsr/{version}/{path.name}")


def pkg_prints(root: Path, version: str, tar: tarfile.TarFile):
    members = [
        ("FingerPRINTShierarchy21Feb2012", "FingerPRINTShierarchy.db"),
        ("prints42_0.pval_blos62", "prints.pval"),
    ]

    for src, dst in members:
        path = root / "prints" / version / src
        tar.add(path, arcname=f"prints/{version}/{dst}")


def pkg_prosite(root: Path, version: str, tar: tarfile.TarFile):
    path = root / "prosite" / version / "evaluator.dat"
    tar.add(path, arcname=f"prosite/{version}/{path.name}")

    path = root / "prosite" / version / "prosite_patterns.dat"
    tar.add(path, arcname=f"prosite/{version}/{path.name}")

    path = root / "prosite" / version / "prosite_profiles"
    tar.add(path, arcname=f"prosite/{version}/{path.name}")

    path = root / "prosite" / version / "skip_flagged_profiles.txt"
    tar.add(path, arcname=f"prosite/{version}/{path.name}")


def pkg_sfld(root: Path, version: str, tar: tarfile.TarFile):
    members = [
        "sfld.hmm",
        "sfld_sites.annot",
        "sfld_hierarchy.txt",
    ]

    for member in members:
        path = root / "sfld" / version / member
        tar.add(path, arcname=f"sfld/{version}/{path.name}")


def pkg_smart(root: Path, version: str, tar: tarfile.TarFile):
    members = [
        "smart.HMMs",
        "smart.HMMs.bin",
    ]

    for member in members:
        path = root / "smart" / version / member
        tar.add(path, arcname=f"smart/{version}/{path.name}")


def pkg_superfamily(root: Path, version: str, tar: tarfile.TarFile):
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
        tar.add(path, arcname=f"superfamily/{version}/{path.name}")
