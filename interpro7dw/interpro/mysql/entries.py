import pickle
from typing import Optional

import MySQLdb

from interpro7dw import pfam
from interpro7dw.utils import logger
from interpro7dw.utils.mysql import uri2dict
from interpro7dw.utils.store import BasicStore
from .utils import jsonify


def populate_annotations(uri: str, hmms_file: str, pfam_alignments: str):
    logger.info("creating webfront_entryannotation")
    con = MySQLdb.connect(**uri2dict(uri))
    cur = con.cursor()
    cur.execute("DROP TABLE IF EXISTS webfront_entryannotation")
    cur.execute(
        """
        CREATE TABLE webfront_entryannotation
        (
            annotation_id INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
            accession VARCHAR(25) NOT NULL,
            type VARCHAR(20) NOT NULL,
            value LONGBLOB NOT NULL,
            mime_type VARCHAR(32) NOT NULL,
            num_sequences INT
        ) CHARSET=utf8mb4 DEFAULT COLLATE=utf8mb4_unicode_ci
        """
    )

    for file in [hmms_file, pfam_alignments]:
        with BasicStore(file, mode="r") as store:
            for accession, anno_type, anno_value, count in store:
                if anno_type == "logo":
                    mime_type = "application/json"
                else:
                    mime_type = "application/gzip"

                cur.execute(
                    """
                    INSERT INTO webfront_entryannotation (
                        accession, type, value, mime_type, num_sequences
                    ) VALUES (%s, %s, %s, %s, %s)
                    """, (accession, anno_type, anno_value, mime_type, count)
                )

    con.commit()

    logger.info("indexing")
    cur.execute(
        """
        CREATE INDEX i_entryannotation 
        ON webfront_entryannotation (accession)
        """
    )

    cur.close()
    con.close()

    logger.info("done")





# def _create_record(entry: Entry, wiki_entry: Optional[dict],
#                    details: Optional[dict],
#                    taxa_tree: Optional[dict]) -> tuple:
#     return (
#         None,
#         entry.accession,
#         entry.type.lower(),
#         entry.name,
#         entry.short_name,
#         entry.database,
#         jsonify(entry.integrates, nullable=True),
#         entry.integrated_in,
#         jsonify(entry.go_terms, nullable=True),
#         jsonify(entry.descriptions, nullable=True),
#         jsonify(wiki_entry, nullable=True),
#         jsonify(details, nullable=True),
#         jsonify(entry.literature, nullable=True),
#         jsonify(entry.hierarchy, nullable=True),
#         jsonify(entry.xrefs, nullable=True),
#         jsonify(entry.ppi, nullable=True),
#         jsonify(entry.pathways, nullable=True),
#         jsonify(entry.overlaps_with, nullable=True),
#         jsonify(taxa_tree, nullable=True),
#         0,
#         1 if entry.is_public else 0,
#         jsonify(entry.history, nullable=True),
#         entry.creation_date,
#         entry.deletion_date,
#         jsonify(entry.counts, nullable=False)
#     )
#
#
# def insert_entries(ipr_url: str, pfam_url: str, entries_file: str,
#                    entry2xrefs_file: str):
#     logger.info("fetching Wikipedia data for Pfam entries")
#     wiki = pfam.get_wiki(pfam_url)
#
#     logger.info("loading Pfam curation/family details")
#     pfam_details = pfam.get_details(pfam_url)
#
#     logger.info("populating webfront_entry")
#     entries = loadobj(entries_file)
#
#     con = MySQLdb.connect(**url2dict(ipr_url), charset="utf8mb4")
#     cur = con.cursor()
#     cur.execute("DROP TABLE IF EXISTS webfront_entry")
#     cur.execute(
#         """
#         CREATE TABLE webfront_entry
#         (
#             entry_id VARCHAR(10) DEFAULT NULL,
#             accession VARCHAR(25) PRIMARY KEY NOT NULL,
#             type VARCHAR(50) NOT NULL,
#             name LONGTEXT,
#             short_name VARCHAR(100),
#             source_database VARCHAR(10) NOT NULL,
#             member_databases LONGTEXT,
#             integrated_id VARCHAR(25),
#             go_terms LONGTEXT,
#             description LONGTEXT,
#             wikipedia LONGTEXT,
#             details LONGTEXT,
#             literature LONGTEXT,
#             hierarchy LONGTEXT,
#             cross_references LONGTEXT,
#             interactions LONGTEXT,
#             pathways LONGTEXT,
#             overlaps_with LONGTEXT,
#             taxa LONGTEXT,
#             is_featured TINYINT NOT NULL,
#             is_alive TINYINT NOT NULL,
#             history LONGTEXT,
#             entry_date DATETIME NOT NULL,
#             deletion_date DATETIME,
#             counts LONGTEXT NOT NULL
#         ) CHARSET=utf8mb4 DEFAULT COLLATE=utf8mb4_unicode_ci
#         """
#     )
#
#     query = """
#         INSERT INTO webfront_entry
#         VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s,
#           %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
#     """
#
#     with SimpleStore(entry2xrefs_file) as store:
#         for accession, xrefs in store:
#             entry = entries.pop(accession)
#             rec = _create_record(entry, wiki_entry=wiki.get(entry.accession),
#                                  details=pfam_details.get(entry.accession),
#                                  taxa_tree=xrefs["taxa"]["tree"])
#             cur.execute(query, rec)
#
#     # Add entries without cross-references
#     for entry in entries.values():
#         rec = _create_record(entry, wiki_entry=wiki.get(entry.accession),
#                              details=pfam_details.get(entry.accession),
#                              taxa_tree=None)
#         cur.execute(query, rec)
#
#     con.commit()
#     cur.close()
#     con.close()
#
#     logger.info("done")


# def insert_release_notes(stg_url: str, rel_url: str, entries_file: str,
#                          proteomeinfo_file: str, structures_file: str,
#                          taxa_file: str, proteins_file: str, matches_file: str,
#                          proteomes_file: str, **kwargs):
#     logger.info("loading entries")
#     entries = loadobj(entries_file)
#
#     logger.info("loading structures")
#     structures = []
#     protein2structures = {}
#     for pdbe_id, entry in loadobj(structures_file).items():
#         structures.append(entry)
#         for protein_acc in entry["proteins"]:
#             try:
#                 protein2structures[protein_acc].add(pdbe_id)
#             except KeyError:
#                 protein2structures[protein_acc] = {pdbe_id}
#
#     logger.info("loading sequence databases")
#     seq_databases = {}
#     con = MySQLdb.connect(**url2dict(stg_url))
#     cur = con.cursor()
#     cur.execute(
#         """
#         SELECT name, name_long, version
#         FROM webfront_database
#         WHERE name in ('reviewed', 'unreviewed', 'uniprot')
#         """
#     )
#
#     for name, name_long, version in cur:
#         seq_databases[name] = {
#             "name": name_long,
#             "version": version,
#             "total": 0,         # total number of proteins
#             "hit": 0,           # number of proteins with at least one hit
#             "integrated": 0     # number of integrated proteins
#         }
#
#     cur.close()
#     con.close()
#
#     # Number of proteins with GO terms from InterPro
#     uniprot2go = 0
#
#     # Entities found in InterPro
#     integrated_proteomes = set()
#     integrated_structures = set()
#     integrated_taxa = set()
#
#     proteins = Store(proteins_file, "r")
#     matches = Store(matches_file, "r")
#     proteomes = Store(proteomes_file, "r")
#     i = 0
#     for i, (protein_acc, protein) in enumerate(proteins.items()):
#         if (i + 1) % 10e6 == 0:
#             logger.info(f"{i+1:>15,}")
#
#         if protein["reviewed"]:
#             database = seq_databases["reviewed"]
#         else:
#             database = seq_databases["unreviewed"]
#
#         database["total"] += 1
#
#         try:
#             protein_matches = matches[protein_acc]
#         except KeyError:
#             continue  # No matches
#
#         # Protein matched by at least one signature
#         database["hit"] += 1
#
#         is_integrated = False
#         for entry_acc in protein_matches:
#             entry = entries[entry_acc]
#             if entry.database == "interpro":
#                 """
#                 Protein matched by at least one InterPro entry,
#                 i.e. at least one integrated signature
#                 """
#                 is_integrated = True
#
#                 if entry.go_terms:
#                     uniprot2go += 1
#                     break
#
#         if is_integrated:
#             database["integrated"] += 1
#
#             try:
#                 proteome_id = proteomes[protein_acc]
#             except KeyError:
#                 pass
#             else:
#                 integrated_proteomes.add(proteome_id)
#
#             try:
#                 protein_structures = protein2structures[protein_acc]
#             except KeyError:
#                 pass
#             else:
#                 integrated_structures |= protein_structures
#
#             integrated_taxa.add(protein["taxid"])
#
#     logger.info(f"{i + 1:>15,}")
#
#     proteins.close()
#     matches.close()
#     proteomes.close()
#
#     # Sums Swiss-Prot and TrEMBL counts
#     for key in ("total", "hit", "integrated"):
#         seq_databases["uniprot"][key] = (seq_databases["reviewed"][key]
#                                          + seq_databases["unreviewed"][key])
#
#     logger.info("tracking changes since last releases")
#     con = MySQLdb.connect(**url2dict(rel_url), charset="utf8mb4")
#     cur = con.cursor()
#     cur.execute(
#         """
#         SELECT accession, source_database, integrated_id
#         FROM webfront_entry
#         WHERE is_alive = 1
#         """
#     )
#     public_entries = set()
#     public_integrated = set()
#     for entry_acc, database, integrated_in in cur:
#         if database == "interpro":
#             public_entries.add(entry_acc)
#         elif integrated_in:
#             # Signature already integrated in the previous release
#             public_integrated.add(entry_acc)
#
#     cur.execute(
#         """
#         SELECT name, version
#         FROM webfront_database
#         WHERE type = 'entry'
#         """
#     )
#     public_databases = dict(cur.fetchall())
#
#     cur.execute("SELECT * FROM webfront_release_note")
#     prev_releases = cur.fetchall()
#     cur.close()
#     con.close()
#
#     con = MySQLdb.connect(**url2dict(stg_url), charset="utf8mb4")
#     cur = con.cursor()
#     cur.execute(
#         """
#         SELECT name, name_long, version, release_date
#         FROM webfront_database
#         WHERE type = 'entry'
#         """
#     )
#     staging_databases = {row[0]: (row[1], row[2], row[3]) for row in cur}
#
#     interpro_new = []
#     interpro_types = {}
#     member_databases = {}
#     pubmed_citations = set()
#     interpro2go = 0
#     latest_entry = None
#
#     for entry in sorted(entries.values(), key=lambda e: e.creation_date):
#         if not entry.is_public:
#             continue
#
#         if entry.database == "interpro":
#             for pub in entry.literature.values():
#                 if pub["PMID"] is not None:
#                     pubmed_citations.add(pub["PMID"])
#
#             try:
#                 interpro_types[entry.type.lower()] += 1
#             except KeyError:
#                 interpro_types[entry.type.lower()] = 1
#
#             if entry.accession not in public_entries:
#                 interpro_new.append(entry.accession)
#
#             interpro2go += len(entry.go_terms)
#             latest_entry = entry.accession
#         else:
#             try:
#                 obj = member_databases[entry.database]
#             except KeyError:
#                 database, version, _ = staging_databases[entry.database]
#
#                 is_new = is_updated = False
#                 if entry.database not in public_databases:
#                     is_new = True
#                 elif version != public_databases[entry.database]:
#                     is_updated = True
#
#                 obj = member_databases[entry.database] = {
#                     "name": database,
#                     "version": version,
#                     "signatures": 0,
#                     "integrated_signatures": 0,
#                     "recently_integrated": [],
#                     "is_new": is_new,
#                     "is_updated": is_updated,
#                     "sets": set()
#                 }
#
#             obj["signatures"] += 1
#             if entry.integrated_in:
#                 obj["integrated_signatures"] += 1
#
#                 if entry.accession not in public_integrated:
#                     # Recent integration
#                     obj["recently_integrated"].append(entry.accession)
#
#             if entry.clan:
#                 obj["sets"].add(entry.clan["accession"])
#
#     # Transforms sets of clans to counts:
#     for obj in member_databases.values():
#         obj["sets"] = len(obj["sets"])
#
#     # Checks that "integrated" proteomes and taxa are valid
#     proteomes = set(loadobj(proteomeinfo_file).keys())
#     errors = integrated_proteomes - proteomes
#     if errors:
#         raise RuntimeError(f"invalid proteomes: {', '.join(errors)}")
#
#     taxa = set(loadobj(taxa_file).keys())
#     errors = integrated_taxa - taxa
#     if errors:
#         raise RuntimeError(f"invalid taxa: {', '.join(errors)}")
#
#     logger.info("creating webfront_release_note")
#     cur.execute("DROP TABLE IF EXISTS webfront_release_note")
#     cur.execute(
#         """
#         CREATE TABLE webfront_release_note
#         (
#             version VARCHAR(20) PRIMARY KEY NOT NULL,
#             release_date DATETIME NOT NULL,
#             content LONGTEXT NOT NULL
#         ) CHARSET=utf8mb4 DEFAULT COLLATE=utf8mb4_unicode_ci
#         """
#     )
#
#     # Adds previous release notes
#     cur.executemany(
#         """
#         INSERT INTO webfront_release_note
#         VALUES (%s, %s, %s)
#         """, prev_releases
#     )
#     con.commit()
#
#     content = {
#         "notes": kwargs.get("notes", []),
#         "interpro": {
#             "entries": sum(interpro_types.values()),
#             "new_entries": interpro_new,
#             "latest_entry": latest_entry,
#             "types": interpro_types,
#             "go_terms": interpro2go,
#             "uniprot2go": uniprot2go
#         },
#         "member_databases": member_databases,
#         "proteins": seq_databases,
#         "structures": {
#             "total": len(structures),
#             "integrated": len(integrated_structures),
#             "version": max(entry["date"]
#                            for entry in structures).strftime("%Y-%m-%d")
#         },
#         "proteomes": {
#             "total": len(proteomes),
#             "integrated": len(integrated_proteomes),
#             "version": seq_databases["uniprot"]["version"]
#         },
#         "taxonomy": {
#             "total": len(taxa),
#             "integrated": len(integrated_taxa),
#             "version": seq_databases["uniprot"]["version"]
#         },
#         "citations": len(pubmed_citations)
#     }
#
#     _, version, date = staging_databases["interpro"]
#     cur.execute(
#         """
#         SELECT COUNT(*)
#         FROM webfront_release_note
#         WHERE version = %s
#         """, (version,)
#     )
#     n_rows, = cur.fetchone()
#
#     if n_rows:
#         # Release notes already in the table: update it
#         cur.execute(
#             """
#             UPDATE webfront_release_note
#             SET content = %s
#             WHERE version = %s
#             """, (jsonify(content), version)
#         )
#     else:
#         # Adds new release notes
#         cur.execute(
#             """
#             INSERT INTO webfront_release_note
#             VALUES (%s, %s, %s)
#             """, (version, date, jsonify(content))
#         )
#
#     con.commit()
#     cur.close()
#     con.close()
#
#     logger.info("done")
