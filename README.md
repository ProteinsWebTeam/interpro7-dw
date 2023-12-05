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
python setup.py install
```

## Configuration

Copy or edit `config.ini` to set the options described below.

### release

| Option        | Description                                                     |
|---------------|-----------------------------------------------------------------|
| version       | InterPro release version, e.g. 80.0                             |
| date          | InterPro release date (expected format: YYYY-MM-DD)             |
| update        | Should the Oracle database be updated with `version` and `date` |

> **Note**: `update` accepts the following values: `yes`, `true`, `no`, and `false` (case-insensitive).

### data

| Option     | Description                                              |
|------------|----------------------------------------------------------|
| alphafold  | Mapping between UniProtKB accessions and AlpaFold IDs    |
| metacyc    | MetaCyc data file (expects a `.tar.gz` archive)          |
| path       | Directory used to store data files                       |
| tmp        | Directory used for temporary files created in each task  |

### databases

Expected format: `user/password@host:port/schema`.

| Option              | Description                                |
|---------------------|--------------------------------------------|
| interpro_production | Oracle production database (interpro user) |
| iprscan_production  | Oracle production database (iprscan user)  |
| interpro_staging    | InterPro release/staging MySQL database    |
| interpro_release    | InterPro release/offsite MySQL database    |
| interpro_fallback   | InterPro release/fallback MySQL database   |
| goa                 | GOA Oracle database                        |
| intact              | IntAct Oracle database                     |
| pdbe                | PDBe Oracle database                       |
| pfam                | Pfam release MySQL database                |
| swissprot           | UniProtKB/Swiss-Prot Oracle database       |

### elasticsearch

Each key/value pair in this section corresponds to an Elastic cluster identifier (e.g. `dev`) and its nodes separated by commas (expected format: `host[:port]`).

e.g.
```
[elasticsearch]
release =  interpro-rl-01:9200,interpro-rl-02:9200,interpro-rl-03:9200
fallback = interpro-fb-01:9200,interpro-fb-02:9200
```

### exchange

| Option    | Description                                                                                                                |
|-----------|----------------------------------------------------------------------------------------------------------------------------|
| ebisearch | Directory monitored by EBI Search to index cross-references                                                                |
| goa       | Directory for mappings required by the GOA team                                                                            |
| interpro  | Directory for archived FTP files (should not finish with the release number, as `release.version` is appended at run time) |
| pdbe      | Directory for mappings required by the PDBe team                                                                           |

### email

Use to send emails to people/groups. As of May 2021, only used during the `notify-curators` steps, to inform curators they can resume using the production database.

| Option | Description                           |
|--------|---------------------------------------|
| server | SMTP sever (format: `host:port`       |
| from   | Sender address                        |
| to     | Comma-separated addressees' addresses |

### workflow

| Option    | Description                                                            |
|-----------|------------------------------------------------------------------------|
| path      | Directory for job input/output files                                   |
| scheduler | Scheduler and queue (format: `scheduler:queue`, e.g. `lsf:production`) |

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
| export-uniparc             | Export all member database matches again UniParc                             |
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
| export-structures-xml | Export an XML file of structural matches (PDBe, CATH, SCOP)                                                    |
| export-uniparc-xml    | Export a `tar.gz` archive of all UniParc matches                                                               |

**Exporting files for internal use (other EMBL-EBI groups)**

| Task name         | Description                                                                                      |
|-------------------|--------------------------------------------------------------------------------------------------|
| export-ebisearch  | Export JSON files of InterPro entries and member database signatures, and their cross-references |
| publish-ebisearch | Move JSON files created in `ebisearch` to a directory monitored by EBI Search                    |
| export-goa        | Export mappings between PDBe, InterPro, GO, and UniProt                                          |
| publish-goa       | Move files to the directory monitored by the GOA team                                            |
| export-pdbe       | Export mappings between PDBe and InterPro                                                        |
| publish-pdbe      | Move files to the directory monitored by the PDBe team                                           |

**Building Oracle tables for the match look-up sercice**

| Task name                 | Description                                           |
|---------------------------|-------------------------------------------------------|
| build-upi-md5             | Build table of protein MD5 checksums                  |
| build-lookup-tab          | Build table of protein matches                        |
| build-lookup-tab-idx      | Index table of protein matches                        |
| build-site-lookup-tab     | Build table of site annotations                       |
| build-site-lookup-tab_idx | Index table of site annotations                       |

**Others**

| Task name            | Description                                                                                                     |
|----------------------|-----------------------------------------------------------------------------------------------------------------|
| notify-curators      | Informing curators that tasks relying on the production database are done, so they can resume curating entries  |

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

**interprodw-dropdb**

```
usage: interprodw-dropdb [-h] config.ini {release,fallback}

Drop release/fallback MySQL database

positional arguments:
  config.ini          configuration file
  {release,fallback}

optional arguments:
  -h, --help          show this help message and exit
```
