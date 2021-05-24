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
- Python 3.7+
- Python packages `cx_Oracle`, `elasticsearch`, and `mysqlclient`
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
| path       | Directory used to store data files                       |
| tmp        | Directory used for temporary files created in each task  |

### databases

Expected format: `user/password@host:port/schema`.

| Option     | Description                                 |
| -----------|---------------------------------------------|
| production | InterPro Oracle production database         |
| staging    | InterPro release/staging MySQL database     |
| release    | InterPro release/offsite MySQL database     |
| fallback   | InterPro release/fallback MySQL database    |
| pfam       | Pfam release MySQL database                 |

### elasticsearch

Each key/value pair in this section corresponds to an Elastic cluster identifier (e.g. `dev`) and its nodes separated by commas (expected format: `host[:port]`).

e.g.
```
[elasticsearch]
release =  interpro-rl-01:9200,interpro-rl-02:9200,interpro-rl-03:9200
fallback = interpro-fb-01:9200,interpro-fb-02:9200
```

### exchange

| Option   | Description                                                                                                                |
| ---------|----------------------------------------------------------------------------------------------------------------------------|
| ebisearch| Directory monitored by EBI Search to index cross-references                                                                |
| goa      | Directory for mappings required by the GOA team                                                                            |
| interpro | Directory for archived FTP files (should not finish with the release number, as `release.version` is appended at run time) |
| pdbe     | Directory for mappings required by the PDBe team                                                                           |

### metacyc

| Option   | Description                                                                                                                |
| ---------|----------------------------------------------------------------------------------------------------------------------------|
| path     | Path to the MetaCyc data file (expects a `.tar.gz` archive)                                                                |

### email

Use to send emails to people/groups. As of August 2020, only use during the `notify-curators` steps, to inform curators they can resume using the production database.

| Option   | Description                           |
| ---------|---------------------------------------|
| server   | SMTP host                             |
| port     | SMTP port number                      |
| from     | Sender address                        |
| to       | Comma-separated addressees' addresses |

### workflow

| Option             | Description                           |
| -------------------|---------------------------------------|
| path               | Directory for job input/output files  |
| lsf_queue          | Name of the queue to submit jobs to   |

## Workflow Description

**Exporting data from Oracle**

| Task name         | Description                                                                   |
|-------------------|-------------------------------------------------------------------------------|
| export-clans      | Export clan information, including profile-profile alignments                 |
| export-proteomes  | Export proteomes data                                                         |
| export-structures | Export structures from the PDBe database                                      |
| export-taxonomy   | Export taxonomic data                                                         |
| init-export       | Split proteins into chunks to avoid having to load all proteins into memory   |
| export-proteins   | Export protein information such a taxon ID, length, UniProt identifier, etc.  |
| uniprot2comments  | Export Swiss-Prot function comments                                           |
| uniprot2evidence  | Export UniProt evidences and genes                                            |
| uniprot2features  | Export sequence feature matches (MobiDB-Lite, TMHMM, Phobius, Coils)          |
| uniprot2matches   | Export protein matches from member databases                                  |
| uniprot2name      | Export UniProt descriptions/names                                             |
| uniprot2proteome  | Export UniProt-proteome mapping                                               |
| uniprot2sequence  | Export protein sequences from UniParc                                         |
| export-entries    | Export InterPro entries and member database signatures                        |

**Creating/populating MySQL tables**

| Task name           | Description                                                                                            |
|---------------------|--------------------------------------------------------------------------------------------------------|
| insert-isoforms     | Insert alternatively spliced isoforms                                                                  |
| insert-databases    | Update InterPro version in Oracle, and insert database information in MySQL                            |
| insert-annotations  | Insert Pfam sequence alignments, profile HMMs, and logos from profile HMMs                             |
| insert-entries      | Insert InterPro entries, and member database signatures                                                |
| insert-clans        | Insert clans                                                                                           |
| insert-proteins     | Insert UniProt proteins with enriched information (e.g. residue annotations, structural features)      |
| insert-proteomes    | Insert UniProt proteomes                                                                               |
| insert-structures   | Insert PDBe structures with enriched information (e.g. secondary structures, literature references)    |
| insert-struct-models| Insert structural models from Baker's group                                                            |
| insert-taxonomy     | Insert taxonomic data                                                                                  |
| insert-release-notes| Generate and insert release notes (number of entries, proteins, recent integrations, etc.)             |
| insert-features     | Insert sequence feature matches                                                                        |
| insert-residues     | Insert site matches                                                                                    |

**Creating Elasticsearch clusters**

In the following tasks, *<id>* represents the cluster identifier, as defined in the [config file](#elasticsearch)

| Task name            | Description                                                         |
|----------------------|---------------------------------------------------------------------|
| es-export            | Export documents for Elasticsearch (IDA and relationship documents) |
| es-init-<id>         | Create the staging indexes on cluster <id>                          |
| es-index-<id>        | Index documents on cluster <id>                                     |
| es-publish-<id>      | Make staging indexes live on cluster <id>                           |

**Exporting files for the public FTP**

| Task name            | Description                                                                                                    |
|----------------------|----------------------------------------------------------------------------------------------------------------|
| export-features-xml  | Export an XML file of sequence feature matches (MobiDB-Lite, TMHMM, Phobius, Coils)                            |
| export-flat-files    | Export flat files (list of entries, InterPro-GO mapping, protein matches, etc.)                                |
| export-interpro-xml  | Export an XML file of InterPro entries and their annotations (e.g. abstract, member database signatures, etc.) |
| export-matches-xml   | Export an XML file of member databases protein matches                                                         |
| export-release-notes | Export a text file containing the release notes                                                                |
| export-structures-xml| Export an XML file of structural matches (PDBe, CATH, SCOP)                                                    |
| export-uniparc-xml   | Export a `tar.gz` archive of all UniParc matches                                                               |

**Exporting files for internal use (other EMBL-EBI groups)**

| Task name         | Description                                                                                      |
|-------------------|--------------------------------------------------------------------------------------------------|
| export-ebisearch  | Export JSON files of InterPro entries and member database signatures, and their cross-references |
| publish-ebisearch | Move JSON files created in `ebisearch` to a directory monitored by EBI Search                    |
| export-goa        | Export mappings between PDBe, InterPro, GO, and UniProt                                          |
| publish-goa       | Move files to the directory monitored by the GOA team                                            |
| export-pdbe       | Export mappings between PDBe and InterPro                                                        |
| publish-pdbe      | Move files to the directory monitored by the PDBe team                                           |

**Others**

| Task name            | Description                                                                                                     |
|----------------------|-----------------------------------------------------------------------------------------------------------------|
| notify-curators      | Informing curators that tasks relying on the production database are done, so they can resume curating entries  |

## Usage

**interpro7dw-build**

```
usage: interpro7dw-build [-h] [-t [TASK [TASK ...]]] [--dry-run] [--detach] [-v] config.ini

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

**interpro7dw-dblinks**

```
usage: interpro7dw-dblinks [-h] config.ini

Test Oracle public database links

positional arguments:
  config.ini          configuration file

optional arguments:
  -h, --help          show this help message and exit
```

**interpro7dw-dropdb**

```
usage: interpro7dw-dropdb [-h] config.ini {release,fallback}

Drop release/fallback MySQL database

positional arguments:
  config.ini          configuration file
  {release,fallback}

optional arguments:
  -h, --help          show this help message and exit
```
