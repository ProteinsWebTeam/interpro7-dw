import json
import math
import os
import pickle
import shutil
from typing import Optional
from xml.sax.saxutils import escape

from interpro7dw.utils import logger
from interpro7dw.interpro.oracle.entries import Entry
from interpro7dw.utils.store import Directory, BasicStore


def _init_fields(entry: Entry, clan_acc: Optional[str],
                 integrates: dict[str, list[str]],
                 relationships: list[str]) -> tuple[list, list]:
    fields = [
        {
            "name": "id",
            "value": entry.accession
        },
        {
            "name": "name",
            "value": entry.name or entry.accession
        },
        {
            "name": "short_name",
            "value": entry.short_name or entry.accession
        },
        {
            "name": "type",
            "value": entry.type.lower()
        },
        {
            "name": "creation_date",
            "value": entry.creation_date.strftime("%Y-%m-%d")
        },
        {
            "name": "description",
            "value": escape(' '.join(entry.descriptions))
        },
        {
            "name": "source_database",
            "value": entry.database
        }
    ]
    xrefs = []

    if clan_acc:
        # Add clan
        fields.append({
            "name": "set",
            "value": clan_acc
        })

    if entry.database.lower() == "interpro":
        # InterPro entry

        for database, signatures in integrates.items():
            # Add contributing database to fields
            fields.append({
                "name": "contributing_database",
                "value": database
            })

            for signature_acc in signatures:
                xrefs.append({
                    "dbname": database.upper(),
                    "dbkey": signature_acc
                })

        for database, references in entry.cross_references.items():
            for ref_id in references:
                xrefs.append({
                    "dbname": database.upper(),
                    "dbkey": ref_id
                })

        for citation in entry.literature.values():
            try:
                pmid = citation["PMID"]
            except KeyError:
                continue
            else:
                xrefs.append({
                    "dbname": "PUBMED",
                    "dbkey": pmid
                })

        for term in entry.go_terms:
            xrefs.append({
                "dbname": "GO",
                "dbkey": term["identifier"]
            })

        for rel_acc in relationships:
            xrefs.append({
                "dbname": "INTERPRO",
                "dbkey": rel_acc
            })
    else:
        # Member database signature
        if entry.integrated_in:
            xrefs.append({
                "dbname": "INTERPRO",
                "dbkey": entry.integrated_in
            })

        for citation in entry.literature.values():
            try:
                pmid = citation["PMID"]
            except KeyError:
                continue
            else:
                xrefs.append({
                    "dbname": "PUBMED",
                    "dbkey": pmid
                })

    return fields, xrefs


def export(clans_file: str, databases_file: str, entries_file: str,
           taxa_file: str, entry2xrefs_file: str, outdir: str,
           fields_per_file: int = 1000000):
    """Creates JSON files containing entries (InterPro + signatures) and
    cross-references to be ingested by EBISearch

    :param clans_file: File of clans information
    :param databases_file: File of databases information
    :param entries_file: File of InterPro entries
    and member DB signatures
    :param taxa_file: File of taxonomic information
    :param entry2xrefs_file: File of entries cross-references
    :param outdir: Output directory
    :param fields_per_file: Maximum number of fields in a JSON file
    """
    logger.info("loading clan members")
    entry2clan = {}
    with open(clans_file, "rb") as fh:
        for clan_acc, clan in pickle.load(fh).items():
            for entry_acc, _, _, _, _ in clan["members"]:
                entry2clan[entry_acc] = clan_acc

    logger.info("loading databases information")
    release_version = release_date = None
    with open(databases_file, "rb") as fh:
        for db in pickle.load(fh).values():
            if db["name"].lower() == "interpro":
                release_version = db["release"]["version"]
                release_date = db["release"]["date"].strftime("%Y-%m-%d")
                break

    if release_version is None:
        raise RuntimeError("missing release version/date for InterPro")

    logger.info("loading taxonomic information")
    taxon_names = {}
    with open(taxa_file, "rb") as fh:
        for taxon_id, taxon in pickle.load(fh).items():
            taxon_names[taxon_id] = taxon["sci_name"]

    with open(entries_file, "rb") as fh:
        entries: dict[str, Entry] = pickle.load(fh)

    integrates = {}  # InterPro accession > member database > sig. accession
    relationships = {}  # InterPro accession > InterPro parent + children
    for entry in entries.values():
        if entry.parent:
            if entry.accession in relationships:
                relationships[entry.accession].append(entry.parent)
            else:
                relationships[entry.accession] = [entry.parent]

            if entry.parent in relationships:
                relationships[entry.parent].append(entry.accession)
            else:
                relationships[entry.parent] = [entry.accession]
        elif entry.integrated_in:
            try:
                mem_dbs = integrates[entry.integrated_in]
            except KeyError:
                mem_dbs = integrates[entry.integrated_in] = {}

            try:
                mem_dbs[entry.database].append(entry.accession)
            except KeyError:
                mem_dbs[entry.database] = [entry.accession]

    try:
        shutil.rmtree(outdir)
    except FileNotFoundError:
        pass

    logger.info("starting")
    i = 0
    step = milestone = math.ceil(0.1 * len(entries))
    types = {}
    num_fields_by_type = {}
    with BasicStore(entry2xrefs_file, mode="r") as store:
        for entry_acc, entry_xrefs in store:
            entry = entries.pop(entry_acc)
            fields, xrefs = _init_fields(entry, entry2clan.get(entry_acc),
                                         integrates.get(entry_acc, {}),
                                         relationships.get(entry_acc, []))

            proteins = entry_xrefs["proteins"]
            for uniprot_acc, uniprot_id, in_alphaphold in proteins:
                xrefs.append({
                    "dbname": "UNIPROT",
                    "dbkey": uniprot_acc
                })

                fields.append({
                    "name": "uniprot_id",
                    "value": uniprot_id
                })

                if in_alphaphold:
                    xrefs.append({
                        "dbname": "ALPHAFOLD",
                        "dbkey": uniprot_acc
                    })

            for gene_name in entry_xrefs["genes"]:
                fields.append({
                    "name": "uniprot_gene",
                    "value": gene_name
                })

            for taxon_id in entry_xrefs["taxa"]["all"]:
                xrefs.append({
                    "dbname": "TAXONOMY",
                    "dbkey": taxon_id
                })

                fields.append({
                    "name": "taxonomy_name",
                    "value": taxon_names[taxon_id]
                })

            for upid in entry_xrefs["proteomes"]:
                xrefs.append({
                    "dbname": "PROTEOMES",
                    "dbkey": upid
                })

            for pdbe_id in entry_xrefs["structures"]:
                xrefs.append({
                    "dbname": "PDB",
                    "dbkey": pdbe_id
                })

            if entry_xrefs["enzymes"]:
                for ecno in entry_xrefs["enzymes"]:
                    xrefs.append({
                        "dbname": "EC",
                        "dbkey": ecno
                    })

            for key in ["metacyc", "reactome"]:
                if entry_xrefs[key]:
                    for pathway_id, pathway_name in entry_xrefs[key]:
                        xrefs.append({
                            "dbname": key.upper(),
                            "dbkey": pathway_id
                        })

            entry_type = entry.type.lower()
            try:
                directory, items = types[entry_type]
            except KeyError:
                directory = Directory(root=os.path.join(outdir, entry_type))
                items = []
                types[entry_type] = (directory, items)
                num_fields_by_type[entry_type] = 0

            num_fields = len(fields) + len(xrefs)
            if num_fields_by_type[entry_type] + num_fields >= fields_per_file:
                # Too many fields in memory for this type already
                path = directory.mktemp(suffix=".json")
                with open(path, "wt") as fh:
                    json.dump({
                        "name": "InterPro",
                        "release": release_version,
                        "release_date": release_date,
                        "entry_count": len(items),
                        "entries": items
                    }, fh, indent=4)

                items.clear()
                num_fields_by_type[entry_type] = 0

            items.append({
                "fields": fields,
                "cross_references": xrefs
            })
            num_fields_by_type[entry_type] += num_fields

            i += 1
            if i == milestone:
                logger.info(f"{i:>15,}")
                milestone += step

    for entry in entries.values():
        if not entry.public:
            continue

        fields, xrefs = _init_fields(entry, entry2clan.get(entry_acc),
                                     integrates.get(entry_acc, {}),
                                     relationships.get(entry_acc, []))

        entry_type = entry.type.lower()
        try:
            directory, items = types[entry_type]
        except KeyError:
            directory = Directory(root=os.path.join(outdir, entry_type))
            items = []
            types[entry_type] = (directory, items)

        items.append({
            "fields": fields,
            "cross_references": xrefs
        })
        i += 1

    logger.info(f"{i:>15,}")

    for entry_type, (directory, items) in types.items():
        if items:
            path = directory.mktemp(suffix=".json")
            with open(path, "wt") as fh:
                json.dump({
                    "name": "InterPro",
                    "release": release_version,
                    "release_date": release_date,
                    "entry_count": len(items),
                    "entries": items
                }, fh, indent=4)

    logger.info("complete")


def publish(src: str, dst: str):
    try:
        shutil.rmtree(dst)
    except FileNotFoundError:
        pass
    finally:
        shutil.copytree(src, dst)
