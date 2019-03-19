# i7dw

Command line utilities for building InterPro7 data warehouse.

## Overview

This repository contains the code used to build the data warehouse for InterPro's new website (_InterPro 7_, in our vernacular). The data warehouse is based on MySQL and Elasticsearch. Data is exported from InterPro's Oracle production database to binary, compressed, indexed files in order to perform out-of-core operations, by using the EBI cluster instead of relying exclusively on Oracle.

## Getting started

### Prerequisites

Requirements:
- Python 3.6+
- Python packages `cx_Oracle`, `elasticsearch`, and `mysqlclient`
- Python package [`mundone`](https://github.com/matthiasblum/mundone/)

### Installing

```bash
git clone https://github.com/matthiasblum/i7dw.git
cd i7dw/
python setup.py install
```

## Configuration

Copy or edit `config.ini` to set the options described below.

### meta

| Option        | Description                         | Notes                       |
|---------------|-------------------------------------|-----------------------------|
| name          | Project name for EBI Search         | Should always be _InterPro_ |
| release       | InterPro release version, e.g. 70.0 |                             |
| release_date  | InterPro release date               | Expected format: YYYY-MM-DD |

### databases

For connection strings, the expected format is: `driver:user/password@host:port/schema`. With `driver` being `mysql` or `oracle`.

| Option             | Description                                              | Notes                       |
| -------------------|----------------------------------------------------------|-----------------------------|
| interpro_oracle    | Connection string to InterPro Oracle database            |                             |
| interpro_mysql_stg | Connection string to InterPro __staging__ MySQL database |                             |
| interpro_mysql_rel | Connection string to InterPro __release__ MySQL database |                             |
| pdbe_oracle        | Connection string to PDBe Oracle database                |                             |
| pfam_mysql         | Connection string to Pfam release MySQL database         |                             |

### jaccard

| Option             | Description                                              | Notes                         |
| -------------------|----------------------------------------------------------|-------------------------------|
| threshold          | Threshold for Jaccard similarity/containment index       | _0.75_ was decided by curators|

### export

| Option             | Description                                                | Notes                         |
| -------------------|------------------------------------------------------------|-------------------------------|
| dir                | Output directory for binary files containing exported data |                               |

### elastic

| Option             | Description                                                | Notes                         |
| -------------------|------------------------------------------------------------|-------------------------------|
| clusters           | Semi-colon separated Elastic clusters. For each cluster, at least one node (`host[:port]`) must be provided. Multiple nodes from the same cluster can be provided (separated by comma) | e.g. _ES1-node1,ES1-node2;ES2-node1_ |
| type               | Document type                                              | Should always be _relationship_ |
| body               | Path to the JSON file defining type mappings, and indices settings (e.g. analyzers) | |
| indices            | Path to the JSON file containing the custom number of shards for some Elastic indices | We create one index per member database, with a fixed number of shards by default; but larger databases may profit from a greater number of shards |
| shards             | Default number of shards | _5_ as of January 2019 |
| dir            | Output directory for JSON files containing documents to index | JSON files will be created in the `documents` sub-directory. For each cluster, a `cluster-N` sub-direcotry will be created. |
| alias            | Index alias, i.e. synonym for all the member database indices | _current_ as of January 2019 |

### ebisearch

| Option             | Description                                                 | Notes                         |
| -------------------|-------------------------------------------------------------|-------------------------------|
| dir                | Output directory for JSON files containing cross-references | Any file in the directory will be __deleted__ (but not the directory itself) |

### workflow

| Option             | Description                                                 | Notes                                              |
| -------------------|-------------------------------------------------------------|----------------------------------------------------|
| queue              | LSF queue name                                              |                                                    |
| dir                | Directory for job input/output files                        | You shouldn't need to change this between releases |

## Usage

### Steps

#### Exporting data from Oracle

| Task name         | Description                                                                                 |
|-------------------|---------------------------------------------------------------------------------------------|
| chunk-proteins    | Split proteins into chunks to avoid having to load all proteins into memory                 |
| export-features   | Dump sequence feature matches, e.g. TMHMM, Phobius, Coils, etc. (__not__ MobiDB-Lite)       |
| export-matches    | Dump protein matches, and MobiDB-Lite sequence features |
| export-residues   | Dump site matches, i.e. residue annotations |
| export-proteins   | Dump protein information such a taxon ID, length, UniProt identifier, etc. |
| export-sequences  | Dump protein sequences from UniParc |
| export-comments   | Dump Swiss-Prot function comments |
| export-names      | Dump UniProt descriptions |
| export-misc       | Dump UniProt evidences and genes |
| export-proteomes  | Dump UniProt proteomes |

#### Creating/populating MySQL tables

| Task name         | Description                                                                                 |
|-------------------|---------------------------------------------------------------------------------------------|
| init-tables       | Drop existing tables and recreate them                |
| insert-taxa       | Load taxonomy data       |
| insert-proteomes  | Load UniProt proteomes |
| insert-databases  | Load database information such as short/long name, version, previous version, description, etc. |
| insert-entries    | Load InterPro entries, and member database signatures |
| insert-annotations| Load Pfam signature annotations (HMM logo) |
| insert-structures | Load PDBe structures |
| insert-sets       | Load sets (e.g. Pfam clans, CDD superfamilies) and profile-profile alignments |
| insert-proteins   | Load proteins with enriched information (e.g. residue annotations, structural features/predictions) |
| release-notes     | Generate release notes, i.e. compare the number of entries in this release and in the previous release |

#### Cross data and update MySQL tables

#### EBI Search

#### Elasticsearch
