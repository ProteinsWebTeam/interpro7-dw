# -*- coding: utf-8 -*-

import MySQLdb

from interpro7dw import logger
from interpro7dw.ebi import pfam
from interpro7dw.ebi.interpro.utils import Table
from interpro7dw.utils import DumpFile
from interpro7dw.utils import loadobj, url2dict
from .utils import jsonify, reduce


def insert_entries(pfam_url: str, stg_url: str, p_entries: str,
                   p_entry2xrefs: str, p_uniprot_models: str):
    logger.info("fetching Wikipedia data for Pfam entries")
    wiki = pfam.get_wiki(pfam_url)

    logger.info("loading Pfam curation/family details")
    pfam_details = pfam.get_details(pfam_url)

    logger.info("populating webfront_entry")
    entries = loadobj(p_entries)

    # Load UniProt entries having a structural model
    uniprot_models = set()
    if p_uniprot_models:
        with open(p_uniprot_models, "rt") as fh:
            for uniprot_acc in map(str.rstrip, fh):
                uniprot_models.add(uniprot_acc)

    con = MySQLdb.connect(**url2dict(stg_url), charset="utf8mb4")
    cur = con.cursor()
    cur.execute("DROP TABLE IF EXISTS webfront_entry")
    cur.execute(
        """
        CREATE TABLE webfront_entry
        (
            entry_id VARCHAR(10) DEFAULT NULL,
            accession VARCHAR(25) PRIMARY KEY NOT NULL,
            type VARCHAR(50) NOT NULL,
            name LONGTEXT,
            short_name VARCHAR(100),
            source_database VARCHAR(10) NOT NULL,
            member_databases LONGTEXT,
            integrated_id VARCHAR(25),
            go_terms LONGTEXT,
            description LONGTEXT,
            wikipedia LONGTEXT,
            details LONGTEXT,
            literature LONGTEXT,
            hierarchy LONGTEXT,
            cross_references LONGTEXT,
            interactions LONGTEXT,
            pathways LONGTEXT,
            overlaps_with LONGTEXT,
            is_featured TINYINT NOT NULL,
            is_alive TINYINT NOT NULL,
            history LONGTEXT,
            entry_date DATETIME NOT NULL,
            deletion_date DATETIME,
            counts LONGTEXT NOT NULL
        ) CHARSET=utf8mb4 DEFAULT COLLATE=utf8mb4_unicode_ci
        """
    )

    """
    Count number of structural models per entry
    (Right now we have only tRosetta)
    """
    cur.execute(
        """
        SELECT algorithm, accession, COUNT(*)
        FROM webfront_structuralmodel
        GROUP BY algorithm, accession
        """
    )
    struct_models_algorithms = {}
    for algorithm, accession, count in cur:
        try:
            struct_models_algorithms[algorithm][accession] = count
        except KeyError:
            struct_models_algorithms[algorithm] = {accession: count}

    cur.close()

    sql = """
        INSERT INTO webfront_entry
        VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s,
          %s, %s, %s, %s, %s, %s, %s, %s, %s)
    """

    with Table(con, sql) as table:
        with DumpFile(p_entry2xrefs) as df:
            for accession, xrefs in df:
                num_struct_models = {}
                for algorithm, counts in struct_models_algorithms.items():
                    num_struct_models[algorithm] = counts.get(accession, 0)

                num_struct_models["full_length"] = 0
                for uniprot_acc, _ in xrefs["proteins"]:
                    if uniprot_acc in uniprot_acc:
                        num_struct_models["full_length"] += 1

                entry = entries[accession]
                counts = reduce(xrefs)
                counts.update({
                    "interactions": len(entry.ppi),
                    "pathways": sum([len(v) for v in entry.pathways.values()]),
                    "sets": 1 if entry.clan else 0,
                    "structural_models": num_struct_models
                })

                table.insert((
                    None,
                    accession,
                    entry.type.lower(),
                    entry.name,
                    entry.short_name,
                    entry.database,
                    jsonify(entry.integrates),
                    entry.integrated_in,
                    jsonify(entry.go_terms),
                    jsonify(entry.description),
                    jsonify(wiki.get(accession)),
                    jsonify(pfam_details.get(accession)),
                    jsonify(entry.literature),
                    jsonify(entry.hierarchy),
                    jsonify(entry.cross_references),
                    jsonify(entry.ppi),
                    jsonify(entry.pathways),
                    jsonify(entry.overlaps_with),
                    0,
                    0 if entry.is_deleted else 1,
                    jsonify(entry.history),
                    entry.creation_date,
                    entry.deletion_date,
                    jsonify(counts)
                ))

    con.commit()
    con.close()
    logger.info("complete")
