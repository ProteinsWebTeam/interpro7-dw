import json
import shutil
import tarfile
from pathlib import Path
from tempfile import mkstemp

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
    logger.info("Starting")
    data_dir = Path(data_dir)
    outdir = Path(outdir)

    if data_dir.absolute() == outdir.absolute():
        # Ensure we don't delete the source data directory
        raise ValueError(f"`outdir` cannot be the same as `data_dir`")

    try:
        shutil.rmtree(outdir)
    except FileNotFoundError:
        pass

    outdir.mkdir(parents=True)

    # Get member database versions
    con = oracledb.connect(ipr_uri)
    cur = con.cursor()
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

    full_archive = outdir / f"interproscan-data-{ipr_version}.tar.gz"
    with tarfile.open(str(full_archive), "w:gz") as tar1:
        tasks = [
            ("AntiFam", versions["antifam"], "antifam", iter_antifam),
            ("CATH", versions["cathgene3d"], "cath", iter_cath),
            ("CDD", versions["cdd"], "cdd", iter_cdd),
            ("HAMAP", versions["hamap"], "hamap", iter_hamap),
            ("NCBIFAM", versions["ncbifam"], "ncbifam", iter_ncbifam),
            ("PANTHER", versions["panther"], "panther", iter_panther),
            ("Pfam", versions["pfam"], "pfam", iter_pfam),
            ("PIRSF", versions["pirsf"], "pirsf", iter_pirsf),
            ("PIRSR", versions["pirsr"], "pirsr", iter_pirsr),
            ("PRINTS", versions["prints"], "prints", iter_prints),
            ("PROSITE", versions["prosite"], "prosite", iter_prosite),
            ("SFLD", versions["sfld"], "sfld", iter_sfld),
            ("SMART", versions["smart"], "smart", iter_smart),
            ("SUPERFAMILY", versions["ssf"], "superfamily", iter_superfamily),
        ]

        for db_name, version, subdir, fn in tasks:
            logger.info(f"Archiving {db_name} {version}")

            archive = outdir / subdir / f"{subdir}-{version}.tar.gz"
            archive.parent.mkdir(parents=True, exist_ok=True)
            with tarfile.open(str(archive), "w:gz") as tar2:
                for path, name in fn(data_dir, version):
                    tar1.add(path, arcname=name)
                    tar2.add(path, arcname=name)

        logger.info(f"Archiving InterPro {ipr_version}")
        archive = outdir / "interpro" / f"interpro-{ipr_version}.tar.gz"
        archive.parent.mkdir(parents=True, exist_ok=True)
        with tarfile.open(str(archive), "w:gz") as tar2:
            for path, name in iter_interpro(ipr_uri, goa_uri, ipr_version):
                tar1.add(path, arcname=name)
                tar2.add(path, arcname=name)
                path.unlink()

    logger.info("Done")


def _export_pathways(cur: oracledb.Cursor) -> tuple[Path, Path]:
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

    fd, pathways_file = mkstemp()
    with open(fd, "wt") as fh:
        json.dump(pathways, fh)

    fd, entry2pathways_file = mkstemp()
    with open(fd, "wt") as fh:
        json.dump(interpro2pathways, fh)

    return Path(pathways_file), Path(entry2pathways_file)


def _export_go_terms(cur: oracledb.Cursor, goa_uri: str) -> tuple[Path, Path]:
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

    fd, terms_file = mkstemp()
    with open(fd, "wt") as fh:
        json.dump({"version": version, "terms": terms}, fh)

    fd, entry2terms_file = mkstemp()
    with open(fd, "wt") as fh:
        json.dump(interpro2go, fh)

    return Path(terms_file), Path(entry2terms_file)


def _export_entries(cur: oracledb.Cursor) -> tuple[Path, Path]:
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

    fd, entries_file = mkstemp()
    with open(fd, "wt") as fh:
        json.dump(entries, fh)

    fd, databases_file = mkstemp()
    with open(fd, "wt") as fh:
        # Full databasename (e.g. InterPro, Pfam, CATH-Gene3D) -> version
        json.dump({n: v for _, n, v in databases.values()}, fh)

    return Path(entries_file), Path(databases_file)


def iter_antifam(root: Path, version: str) -> list[tuple[Path, str]]:
    return [
        (root / "antifam" / version / "AntiFam.hmm", f"antifam/{version}/AntiFam.hmm")
    ]


def iter_cath(root: Path, version: str) -> list[tuple[Path, str]]:
    base = root / "cath-gene3d" / version

    path = base / "gene3d_main.hmm"
    members = [(path, f"cath/{version}/gene3d/{path.name}")]

    path = base / "discontinuous" / "discontinuous_regs.pkl.py3"
    members.append((path, f"cath/{version}/gene3d/{path.name}"))

    path = base / "model_to_family_map.tsv"
    members.append((path, f"cath/{version}/gene3d/{path.name}"))

    path = base / "funfam" / "models"
    for child in path.iterdir():
        members.append((child, f"cath/{version}/funfam/{child.name}"))

    return members


def iter_cdd(root: Path, version: str) -> list[tuple[Path, str]]:
    return [
        (root / "cdd" / version / "data", f"cdd/{version}/data"),
        (root / "cdd" / version / "db", f"cdd/{version}/db"),
    ]


def iter_hamap(root: Path, version: str) -> list[tuple[Path, str]]:
    members = []
    for member in ["hamap.prf", "hamap.hmm.lib", "profiles"]:
        path = root / "hamap" / version / member
        members.append((path, f"hamap/{version}/{path.name}"))

    return members


def iter_interpro(ipr_uri: str, goa_uri: str, version: str) -> list[tuple[Path, str]]:
    con = oracledb.connect(ipr_uri)
    cur = con.cursor()
    pathways_file, entry2pathways_file = _export_pathways(cur)
    terms_file, entry2terms_file = _export_go_terms(cur, goa_uri)
    entries_file, databases_file = _export_entries(cur)

    return [
        (pathways_file, f"interpro/{version}/pathways.json"),
        (entry2pathways_file, f"interpro/{version}/pathways.ipr.json"),
        (terms_file, f"interpro/{version}/goterms.json"),
        (entry2terms_file, f"interpro/{version}/goterms.ipr.json"),
        (entries_file, f"interpro/{version}/entries.json"),
        (databases_file, f"interpro/{version}/database.json"),
    ]


def iter_ncbifam(root: Path, version: str) -> list[tuple[Path, str]]:
    return [
        (root / "ncbifam" / version / "ncbifam.hmm", f"ncbifam/{version}/ncbifam.hmm")
    ]


def iter_panther(root: Path, version: str) -> list[tuple[Path, str]]:
    members = [(root / "panther" / version / "famhmm", f"panther/{version}/famhmm")]

    path = root / "panther" / version / "PAINT_Annotations"
    for child in path.iterdir():
        if child.suffix == ".json":
            name = f"panther/{version}/{child.parent.name}/{child.name}"
            members.append((child, name))

    path = root / "panther" / version / "Tree_MSF"
    members.append((path, f"panther/{version}/{path.name}"))
    return members


def iter_pfam(root: Path, version: str) -> list[tuple[Path, str]]:
    members = []
    for member in ["pfam_a.hmm", "pfam_clans", "pfam_a.dat", "pfam_a.seed"]:
        path = root / "pfam" / version / member
        members.append((path, f"pfam/{version}/{path.name}"))

    return members


def iter_pirsf(root: Path, version: str) -> list[tuple[Path, str]]:
    return [
        (root / "pirsf" / version / "pirsf.dat", f"pirsf/{version}/pirsf.dat"),
        (root / "pirsf" / version / "sf_hmm_all", f"pirsf/{version}/sf_hmm_all"),
    ]


def iter_pirsr(root: Path, version: str) -> list[tuple[Path, str]]:
    return [
        (root / "pirsr" / version / "sr_hmm_all", f"pirsr/{version}/sr_hmm_all"),
        (root / "pirsr" / version / "sr_uru.json", f"pirsr/{version}/sr_uru.json"),
    ]


def iter_prints(root: Path, version: str) -> list[tuple[Path, str]]:
    members = []
    for src, dst in [
        ("FingerPRINTShierarchy21Feb2012", "FingerPRINTShierarchy.db"),
        ("prints42_0.pval_blos62", "prints.pval"),
    ]:
        path = root / "prints" / version / src
        members.append((path, f"prints/{version}/{dst}"))

    return members


def iter_prosite(root: Path, version: str) -> list[tuple[Path, str]]:
    members = []
    for member in [
        "evaluator.dat",
        "prosite_patterns.dat",
        "prosite_profiles",
        "skip_flagged_profiles.txt",
    ]:
        path = root / "prosite" / version / member
        members.append((path, f"prosite/{version}/{path.name}"))

    return members


def iter_sfld(root: Path, version: str) -> list[tuple[Path, str]]:
    members = []
    for member in ["sfld.hmm", "sfld_sites.annot", "sfld_hierarchy.txt"]:
        path = root / "sfld" / version / member
        members.append((path, f"sfld/{version}/{path.name}"))

    return members


def iter_smart(root: Path, version: str) -> list[tuple[Path, str]]:
    members = []
    for member in ["smart.HMMs", "smart.HMMs.bin"]:
        path = root / "smart" / version / member
        members.append((path, f"smart/{version}/{path.name}"))

    return members


def iter_superfamily(root: Path, version: str) -> list[tuple[Path, str]]:
    files = [
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
    members = []
    for member in files:
        path = root / "superfamily" / version / member
        members.append((path, f"superfamily/{version}/{path.name}"))

    return members
