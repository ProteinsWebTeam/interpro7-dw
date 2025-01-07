# interpro7-dw

Command line utilities for building InterPro7 data warehouse.

## Overview

This repository contains the code used to build the data warehouse for InterPro's new website (*InterPro 7*, in our vernacular). The data warehouse is based on MySQL and Elasticsearch. Data is exported from InterPro's Oracle production database to binary, compressed, indexed files in order to perform out-of-core operations, by using the EBI cluster instead of relying exclusively on Oracle.

## Table of Contents

1. [Getting Started](#getting-started)
2. [Configuration](#configuration)
3. [Workflow Description](#workflow-description)
4. [Usage](#usage)

## Getting Started

### Prerequisites

Requirements:
- Python 3.11+
- Python packages `oracledb`, `elasticsearch`, and `mysqlclient`
- Python package [`mundone`](https://github.com/matthiasblum/mundone/)

### Installing

```bash
git clone https://github.com/ProteinsWebTeam/interpro7-dw.git
cd interpro7-dw
pip install .
```

## Configuration

Copy or edit `config.toml` to set the options described below.

### release

| Option        | Description                                                     |
|---------------|-----------------------------------------------------------------|
| version       | InterPro release version, e.g. 80.0                             |
| date          | InterPro release date (expected format: YYYY-MM-DD)             |
| update        | Should the Oracle database be updated with `version` and `date` |

### data

| Option    | Description                                             |
|-----------|---------------------------------------------------------|
| alphafold | Mapping between UniProtKB accessions and AlpaFold IDs   |
| metacyc   | MetaCyc data file (expects a `.tar.gz` archive)         |
| path      | Directory used to store data files                      |
| src       | Directory of member database files                      |
| tmp       | Directory used for temporary files created in each task |

### databases

Expected format: `user/password@host:port/schema`.

| Option              | Description                                |
|---------------------|--------------------------------------------|
| interpro.production | Oracle production database (interpro user) |
| iprscan.production  | Oracle production database (iprscan user)  |
| interpro.staging    | InterPro release/staging MySQL database    |
| interpro.release    | InterPro release/offsite MySQL database    |
| goa                 | GOA Oracle database                        |
| pdbe                | PDBe Oracle database                       |
| uniprot             | UniProtKB/Swiss-Prot Oracle database       |

### elasticsearch

For each Elasticsearch cluster, the following options need to be provided:

* `nodes`: list of Elasticsearch nodes (format: `host:port`)
* `user`: Elasticsearch user
* `password`: Password for the user
* `fingerprint`: SSL certificate fingerprint

Assuming we have two clusters (`test` and `prod`), the configuration would be something like:

e.g.
```toml
[elasticsearch]
test.nodes = [ "es-test-node1:9200", "es-test-node2:9200" ]
test.user = "elastic"
test.password = "..."
test.fingerprint = "..."
prod.nodes = [ "es-prod-node1:9200", "es-prod-node2:9200" ]
prod.user = "elastic"
prod.password = "..."
prod.fingerprint = "..."
```

### exchange

| Option    | Description                                                                                                                |
|-----------|----------------------------------------------------------------------------------------------------------------------------|
| ebisearch | Directory monitored by EBI Search to index cross-references                                                                |
| goa       | Directory for mappings required by the GOA team                                                                            |
| interpro  | Directory for archived FTP files (should not finish with the release number, as `release.version` is appended at run time) |

### workflow

| Option    | Description                                                                         |
|-----------|-------------------------------------------------------------------------------------|
| path      | Directory for job input/output files                                                |
| scheduler | Scheduler and queue (format: `scheduler[:queue]`, e.g. `lsf:production` or `slurm`) |

## Workflow Description

**Exporting data from Oracle**

| Task name                  | Description                                                                  |
|----------------------------|------------------------------------------------------------------------------|
| export-clans               | Export clan information, including profile-profile alignments                |
| export-databases           | Export database information                                                  |
| export-entries             | Export InterPro entries and member database signatures                       |
| export-isoforms            | Export Varsplic matches                                                      |
| export-features            | Export sequence feature matches (e.g. Coils, MobiDB-Lite)                    |
| export-matches             | Export protein matches from member databases                                 |
| export-proteins            | Export protein information such a taxon ID, length, UniProt identifier, etc. |
| export-residues            | Export residue annotations (site matches)                                    |
| export-uniparc             | Export all matches found in UniParc                                          |
| export-pdbe-matches        | Export matches against sequences in PDBe                                     |
| export-taxa                | Export taxonomic data                                                        |
| export-structures          | Export PDBe structures                                                       |
| export-structure-chains    | Export the UniProt-PDBe mapping                                              |
| export-pfam-alignments     | Export sequences alignments of Pfam families                                 |
| export-alphafold           | Export the UniProt-AlphaFold mapping                                         |
| export-hmms                | Export models for HMMER3-based member databases                              |
| export-sequences           | Export protein sequences                                                     |
| export-reference-proteomes | Export references proteomes                                                  |
| export-evidences           | Export UniProt evidences/genes                                               |
| export-functions           | Export Swiss-Prot function comments                                          |
| export-names               | Export UniProt descriptions/names                                            |
| export-proteomes           | Export UniProt-proteome mapping                                              |

**Tracking cross-references between clans, entries, proteomes, proteins, structures, and taxa**

| Task name              | Description                                                                                          |
|------------------------|------------------------------------------------------------------------------------------------------|
| export-dom-orgs        | Track the domain architecture/organisation for each UniProt protein                                  |
| export-sim-entries     | Track the similar entries based on how much they overlap in all the proteins they match              |
| export-clan2xrefs      | Clans × (domain organisations, proteomes, proteins, structures, taxa)                                |  
| export-entry2xrefs     | Entries × (AlphaFold models, domain organisations, EC numbers, pathways, proteins, structures, taxa) |
| export-proteome2xrefs  | Proteomes × (domain organisations, clans, entries, proteins, structures, taxa)                       |
| export-structure2xrefs | Structures × (domain organisations, clans, entries, proteomes, proteins, taxa)                       |
| export-taxon2xrefs     | Taxa × (entries, proteomes, proteins, structures, taxa)                                              |

**Creating/populating MySQL tables**

| Task name            | Description                                                                                         |
|----------------------|-----------------------------------------------------------------------------------------------------|
| insert-annotations   | Insert Pfam sequence alignments, profile HMMs, and logos from profile HMMs                          |
| insert-clans         | Insert clans                                                                                        |
| insert-databases     | Insert database information                                                                         |
| insert-entries       | Insert InterPro entries and member database signatures                                              |
| insert-entries-taxa  | Insert the taxonomic distribution of InterPro entries and member database signatures                |
| insert-features      | Insert sequence feature matches                                                                     |
| insert-isoforms      | Insert alternatively spliced isoforms                                                               |
| insert-residues      | Insert site annotations                                                                             |
| insert-proteins      | Insert UniProt proteins with enriched information (e.g. structural features, sequences, etc.)       |
| insert-proteomes     | Insert UniProt reference proteomes                                                                  |
| insert-structures    | Insert PDBe structures with enriched information (e.g. secondary structures, literature references) |
| insert-struct-models | Insert structural models                                                                            |
| insert-taxa          | Insert taxonomic data                                                                               |
| insert-release-notes | Insert (or update) release notes (number of entries, proteins, recent integrations, etc.)           |

**Indexing (some) MySQL tables**

From some tables, the indexes creation is in a different task not to have to re-populate the table from scratch should an error occur.

| Task name          |
|--------------------|
| index-annotations  |
| index-entries      |
| index-features     |
| index-residues     |
| index-proteins     |
| index-taxa         |

**Creating Elasticsearch clusters**

In the following tasks, *id* represents the cluster identifier, as defined in the [config file](#elasticsearch)

| Task name            | Description                                                         |
|----------------------|---------------------------------------------------------------------|
| es-export            | Export documents for Elasticsearch (IDA and relationship documents) |
| es-init-*id*         | Create the staging indexes on cluster *id*                          |
| es-index-*id*        | Index documents on cluster *id*                                     |
| es-publish-*id*      | Make staging indexes live on cluster *id*                           |

**Exporting files for the public FTP**

| Task name             | Description                                                                                                    |
|-----------------------|----------------------------------------------------------------------------------------------------------------|
| export-features-xml   | Export an XML file of sequence feature matches (MobiDB-Lite, TMHMM, Phobius, Coils)                            |
| export-flat-files     | Export flat files (list of entries, InterPro-GO mapping, protein matches, etc.)                                |
| export-interpro-xml   | Export an XML file of InterPro entries and their annotations (e.g. abstract, member database signatures, etc.) |
| export-matches-xml    | Export an XML file of member databases protein matches                                                         |
| export-release-notes  | Export a text file containing the release notes                                                                |
| export-uniparc-xml    | Export a `tar.gz` archive of all UniParc matches                                                               |

**Exporting files for internal use (other EMBL-EBI groups)**

| Task name         | Description                                                                                      |
|-------------------|--------------------------------------------------------------------------------------------------|
| export-ebisearch  | Export JSON files of InterPro entries and member database signatures, and their cross-references |
| publish-ebisearch | Move JSON files created in `ebisearch` to a directory monitored by EBI Search                    |
| export-goa        | Export mappings between PDBe, InterPro, GO, and UniProt                                          |
| publish-goa       | Move files to the directory monitored by the GOA team                                            |

**Tasks related to the match lookup service**

| Task name      | Description                                              |
|----------------|----------------------------------------------------------|
| lookup-md5     | Build Oracle table of protein MD5 hashes (legacy lookup) |
| lookup-matches | Build Oracle table of matches (legacy lookup)            |
| lookup-sites   | Build Oracle table of site matches (legacy lookup)       |
| lookup         | Build lookup database                                    |

## Usage

**interprodw-build**

```
usage: interprodw-build [-h] [-t [TASK [TASK ...]]] [--dry-run] [--detach] [-v] config.ini

Build InterPro7 data warehouse

positional arguments:
  config.ini            configuration file

optional arguments:
  -h, --help            show this help message and exit
  -t [TASK [TASK ...]], --tasks [TASK [TASK ...]]
                        tasks to run
  --dry-run             list tasks to run and exit
  --detach              enqueue tasks to run and exit
  -v, --version         show the version and exit
```
