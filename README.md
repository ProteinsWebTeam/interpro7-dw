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

### release

| Option        | Description                         | Notes                       |
|---------------|-------------------------------------|-----------------------------|
| version       | InterPro release version, e.g. 70.0 |                             |
| date          | InterPro release date               | Expected format: YYYY-MM-DD |

### databases

For connection strings, the expected format is: `user/password@[host:port/]schema`.

| Option              | Description                                                | Notes                       |
| --------------------|------------------------------------------------------------|-----------------------------|
| interpro_production | Connection string to InterPro production database          | Oracle database             |
| interpro_staging    | Connection string to InterPro release/staging database     | MySQL database              |
| interpro_offsite    | Connection string to InterPro release/offsite database     | MySQL database              |
| interpro_fallback   | Connection string to InterPro release/fallback database    | MySQL database              |
| pdbe                | Connection string to PDBe production database              | Oracle database             |
| pfam                | Connection string to Pfam release database                 | MySQL database              |

### ebisearch

| Option    | Description                                          | Notes                                                           |
| ----------|------------------------------------------------------|-----------------------------------------------------------------|
| path-stg  | Directory where to write cross-references JSON files | All files and sub-directories in `path-stg` will be __deleted__ |
| path-rel  | Directory used by EBI Search                         |                                                                 |

### elasticsearch

| Option   | Description                                                                                                  | Notes                                                                                                                   |
| ---------|--------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------|
| nodes    | Elasticsearch nodes grouped in clusters. Clusters are separated by semicolons, nodes are separated by commas | e.g. _Cluster1-node1,Cluster1-node2;Cluster1-node1_                                                                     |
| path     | Directory where to write JSON documents to index                                                             | A `documents` sub-directory will be created. For each cluster, an additional `cluster-N` sub-direcotry will be created. |

### goa

| Option  | Description                                                | Notes                         |
| --------|------------------------------------------------------------|-------------------------------|
| path    | Output directory for mappings required by the GOA team     |                               |

### stores

| Option  | Description                                                | Notes                         |
| --------|------------------------------------------------------------|-------------------------------|
| path    | Output directory for binary files containing exported data |                               |


### workflow

| Option             | Description                           | Notes                                              |
| -------------------|---------------------------------------|----------------------------------------------------|
| path               | Directory for job input/output files  | You shouldn't need to change this between releases |
| queue              | LSF queue name                        |                                                    |

## Usage

### Steps

#### Exporting data from Oracle

| Task name         | Description                                                                                 |
|-------------------|---------------------------------------------------------------------------------------------|
| chunk-proteins    | Split proteins into chunks to avoid having to load all proteins into memory                 |
| export-features   | Dump sequence feature matches, e.g. TMHMM, Phobius, Coils, etc. (__not__ MobiDB-Lite)       |
| export-matches    | Dump protein matches, and MobiDB-Lite sequence features                                     |
| export-residues   | Dump site matches, i.e. residue annotations                                                 |
| export-proteins   | Dump protein information such a taxon ID, length, UniProt identifier, etc.                  |
| export-sequences  | Dump protein sequences from UniParc                                                         |
| export-comments   | Dump Swiss-Prot function comments                                                           |
| export-names      | Dump UniProt descriptions                                                                   |
| export-misc       | Dump UniProt evidences and genes                                                            |
| export-proteomes  | Dump UniProt proteomes                                                                      |

#### Creating/populating small MySQL tables

| Task name         | Description                                                                                     |
|-------------------|-------------------------------------------------------------------------------------------------|
| init-tables       | Drop existing tables and recreate them                                                          |
| insert-taxa       | Load taxonomy data                                                                              |
| insert-proteomes  | Load UniProt proteomes                                                                          |
| insert-databases  | Load database information such as short/long name, version, previous version, description, etc. |
| insert-entries    | Load InterPro entries, and member database signatures                                           |
| insert-annotations| Load Pfam signature annotations (HMM logo)                                                      |
| insert-structures | Load PDBe structures                                                                            |
| insert-sets       | Load sets (e.g. Pfam clans, CDD superfamilies) and profile-profile alignments                   |

#### Preparing additional proteins-related data

| Task name            | Description                                                                                                     |
|----------------------|-----------------------------------------------------------------------------------------------------------------|
| export-ida           | Calculate domain architectures so we can know how many proteins share the same domain architecture              |
| release-notes        | Generate release notes, i.e. compare the number of entries in this release and in the previous release          |
| overlapping-families | Compare how InterPro entries overlap, and define overlapping (super)familiy relationships                       |
| export-xrefs         | Evaluate the cross-references or relationships between entries, proteomes, taxa, structures, proteins, and sets |


| release-notes     | Generate release notes, i.e. compare the number of entries in this release and in the previous release |

#### Populating the protein MySQL table

| Task name            | Description                                                                                                     |
|----------------------|-----------------------------------------------------------------------------------------------------------------|
| insert-proteins      | Load proteins with enriched information (e.g. residue annotations, structural features/predictions)             |

#### Updating MySQL tables

| Task name            | Description                                                                                                     |
|----------------------|-----------------------------------------------------------------------------------------------------------------|
| update-counts        | Update MySQL tables with the number of cross-references each entry, taxon, proteome, structure, and set has     |

#### EBI Search

| Task name            | Description                                                                                                                       |
|----------------------|-----------------------------------------------------------------------------------------------------------------------------------|
| ebi-search           | Generate JSON files of InterPro entries and member database signatures so they can be indexed into EBI Search and made searchable |

#### Elasticsearch

| Task name            | Description                                                                                                                                                                     |
|----------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| init-elastic         | Create the output directory if needed, delete existing Elasticsearch documents if any                                                                                           |
| create-documents     | Generate Elastic documents and store them in JSON files                                                                                                                         |
| index-*N*            | Create indices and index documents by loading JSON files to an Elasticsearch cluster. N is an integer representing the target Elasticsearch cluster                             |
| complete-index-*N*   | Index documents that failed to be indexed during the previous step (often due to network errors). Loop until all documents are indexed or a specific number of tries is reached |
| update-alias-*N*     | Update the alias to use the latest indices. Indices that previously used the alias are deleted                                                                                  |
