import json
import tarfile
from pathlib import Path
from tempfile import mkstemp

import oracledb

from interpro7dw import uniprot
from ..oracle.entries import REPR_DOM_TYPES, REPR_DOM_DATABASES
from interpro7dw.utils import logger


def package_data(ipr_uri: str, goa_uri: str, data_dir: str, output: str):
    logger.info("Exporting JSON files")
    con = oracledb.connect(ipr_uri)
    cur = con.cursor()
    pathways_file, entry2pathways_file = _export_pathways(cur)
    terms_file, entry2terms_file = _export_go_terms(cur, goa_uri)
    entries_file = _export_entries(cur)

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
    with tarfile.open(output, "w:gz") as tar:
        logger.info("Archiving JSON files")
        for file, name in [
            (pathways_file, "pathways.json"),
            (entry2pathways_file, "pathways.ipr.json"),
            (terms_file, "goterms.json"),
            (entry2terms_file, "goterms.ipr.json"),
            (entries_file, "entries.json")
        ]:
            tar.add(file, arcname=f"xrefs/{name}")
            file.unlink()

        logger.info("Archiving AntiFam")
        pkg_antifam(data_dir, versions["antifam"], tar)

        logger.info("Archiving CATH (Gene3D + FunFam)")
        pkg_cath(data_dir, versions["cathgene3d"], tar)

        logger.info("Archiving CDD")
        pkg_cdd(data_dir, versions["cdd"], tar)

        logger.info("Archiving HAMAP")
        pkg_hamap(data_dir, versions["hamap"], tar)

        logger.info("Archiving NCBIfam")
        pkg_ncbifam(data_dir, versions["ncbifam"], tar)

        logger.info("Archiving PANTHER")
        pkg_panther(data_dir, versions["panther"], tar)

        logger.info("Archiving Pfam")
        pkg_pfam(data_dir, versions["pfam"], tar)

        logger.info("Archiving PIRSF")
        pkg_pirsf(data_dir, versions["pirsf"], tar)

        logger.info("Archiving PIRSR")
        pkg_pirsr(data_dir, versions["pirsr"], tar)

        logger.info("Archiving PRINTS")
        pkg_prints(data_dir, versions["prints"], tar)

        logger.info("Archiving PROSITE (Patterns + Profiles)")
        pkg_prosite(data_dir, versions["prosite"], tar)

        logger.info("Archiving SFLD")
        pkg_sfld(data_dir, versions["sfld"], tar)

        logger.info("Archiving SMART")
        pkg_smart(data_dir, versions["smart"], tar)

        logger.info("Archiving SUPERFAMILY")
        pkg_superfamily(data_dir, versions["ssf"], tar)

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
        json.dump({
            "version": version,
            "terms": terms
        }, fh)

    fd, entry2terms_file = mkstemp()
    with open(fd, "wt") as fh:
        json.dump(interpro2go, fh)

    return Path(terms_file), Path(entry2terms_file)


def _export_entries(cur: oracledb.Cursor) -> Path:
    cur.execute("SELECT CODE, ABBREV FROM INTERPRO.CV_ENTRY_TYPE")
    types = dict(cur.fetchall())

    cur.execute(
        """
        SELECT D.DBCODE, LOWER(D.DBSHORT), D.DBNAME, V.VERSION
        FROM INTERPRO.CV_DATABASE D
        INNER JOIN INTERPRO.DB_VERSION V ON D.DBCODE = V.DBCODE
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
        WHERE NOT REGEXP_LIKE(M.METHOD_AC, 'PTHR\d+:SF\d+')
        UNION ALL
        SELECT FM.METHOD_AC, FM.NAME, FM.NAME, 'G', FM.DBCODE, NULL
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
            "database": {
                "name": dbname,
                "version": dbversion
            }
        }

    fd, entries_file = mkstemp()
    with open(fd, "wt") as fh:
        json.dump(entries, fh)

    return Path(entries_file)


def pkg_antifam(root: Path, version: str, tar: tarfile.TarFile):
    path = root / "antifam" / version / "AntiFam.hmm"
    tar.add(path, arcname="antifam/AntiFam.hmm")


def pkg_cath(root: Path, version: str, tar: tarfile.TarFile):
    base = root / "cath-gene3d" / version

    path = base / "gene3d_main.hmm"
    tar.add(path, arcname=f"cath/gene3d/{path.name}")

    path = base / "discontinuous" / "discontinuous_regs.pkl.py3"
    tar.add(path, arcname=f"cath/gene3d/{path.name}")

    path = base / "model_to_family_map.tsv"
    tar.add(path, arcname=f"cath/gene3d/{path.name}")

    path = base / "funfam" / "models"
    for child in path.iterdir():
        tar.add(child, arcname=f"cath/funfam/{child.name}")


def pkg_cdd(root: Path, version: str, tar: tarfile.TarFile):
    path = root / "cdd" / version / "data"
    tar.add(path, arcname="cdd/data")

    path = root / "cdd" / version / "db"
    tar.add(path, arcname="cdd/db")


def pkg_hamap(root: Path, version: str, tar: tarfile.TarFile):
    path = root / "hamap" / version / "hamap.prf"
    tar.add(path, arcname=f"hamap/{path.name}")

    path = root / "hamap" / version / "hamap.hmm.lib"
    tar.add(path, arcname=f"hamap/{path.name}")

    path = root / "hamap" / version / "profiles"
    tar.add(path, arcname=f"hamap/{path.name}")


def pkg_ncbifam(root: Path, version: str, tar: tarfile.TarFile):
    path = root / "ncbifam" / version / "ncbifam.hmm"
    tar.add(path, arcname=f"ncbifam/{path.name}")


def pkg_panther(root: Path, version: str, tar: tarfile.TarFile):
    path = root / "panther" / version / "famhmm"
    tar.add(path, arcname=f"panther/{path.name}")

    path = root / "panther" / version / "PAINT_Annotations"
    for child in path.iterdir():
        if child.suffix == ".json":
            tar.add(child, arcname=f"panther/{child.parent.name}/{child.name}")

    path = root / "panther" / version / "Tree_MSF"
    tar.add(path, arcname=f"panther/{path.name}")


def pkg_pfam(root: Path, version: str, tar: tarfile.TarFile):
    path = root / "pfam" / version / "pfam_a.hmm"
    tar.add(path, arcname=f"pfam/{path.name}")

    path = root / "pfam" / version / "pfam_clans"
    tar.add(path, arcname=f"pfam/{path.name}")

    path = root / "pfam" / version / "pfam_a.dat"
    tar.add(path, arcname=f"pfam/{path.name}")

    path = root / "pfam" / version / "pfam_a.seed"
    tar.add(path, arcname=f"pfam/{path.name}")


def pkg_pirsf(root: Path, version: str, tar: tarfile.TarFile):
    path = root / "pirsf" / version / "pirsf.dat"
    tar.add(path, arcname=f"pirsf/{path.name}")

    path = root / "pirsf" / version / "sf_hmm_all"
    tar.add(path, arcname=f"pirsf/{path.name}")


def pkg_pirsr(root: Path, version: str, tar: tarfile.TarFile):
    path = root / "pirsr" / version / "sr_hmm_all"
    tar.add(path, arcname=f"pirsr/{path.name}")

    path = root / "pirsr" / version / "sr_uru.json"
    tar.add(path, arcname=f"pirsr/{path.name}")


def pkg_prints(root: Path, version: str, tar: tarfile.TarFile):
    members = [
        ("FingerPRINTShierarchy21Feb2012", "FingerPRINTShierarchy.db"),
        ("prints42_0.pval_blos62", "prints.pval"),
        ("prints42_0.kdat", "prints.kdat"),
    ]

    for src, dst in members:
        path = root / "prints" / version / src
        tar.add(path, arcname=f"prints/{dst}")


def pkg_prosite(root: Path, version: str, tar: tarfile.TarFile):
    path = root / "prosite" / version / "evaluator.dat"
    tar.add(path, arcname=f"prosite/{path.name}")

    path = root / "prosite" / version / "prosite_patterns.dat"
    tar.add(path, arcname=f"prosite/{path.name}")

    path = root / "prosite" / version / "prosite_profiles"
    tar.add(path, arcname=f"prosite/{path.name}")

    path = root / "prosite" / version / "skip_flagged_profiles.txt"
    tar.add(path, arcname=f"prosite/{path.name}")


def pkg_sfld(root: Path, version: str, tar: tarfile.TarFile):
    members = [
        "sfld.hmm",
        "sfld_sites.annot",
        "sfld_hierarchy_flat.txt",
    ]

    for member in members:
        path = root / "sfld" / version / member
        tar.add(path, arcname=f"sfld/{path.name}")


def pkg_smart(root: Path, version: str, tar: tarfile.TarFile):
    members = [
        "smart.HMMs",
        "smart.HMMs.bin",
    ]

    for member in members:
        path = root / "smart" / version / member
        tar.add(path, arcname=f"smart/{path.name}")


def pkg_superfamily(root: Path, version: str, tar: tarfile.TarFile):
    members = [
        "hmmlib_1.75",
        "self_hits.tab",
        "dir.cla.scop.txt_1.75",
        "dir.des.scop.txt_1.75",
        "model.tab",
        "pdbj95d",
        "LICENSE",
    ]

    for member in members:
        path = root / "superfamily" / version / member
        tar.add(path, arcname=f"superfamily/{path.name}")
