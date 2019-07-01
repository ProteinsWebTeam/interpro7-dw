from ... import dbms


def init(uri):
    con, cur = dbms.connect(uri)

    cur.execute('DROP TABLE IF EXISTS webfront_structure')
    cur.execute('DROP TABLE IF EXISTS webfront_alignment')
    cur.execute('DROP TABLE IF EXISTS webfront_set')
    cur.execute('DROP TABLE IF EXISTS webfront_entryannotation')
    cur.execute('DROP TABLE IF EXISTS webfront_entry')
    cur.execute('DROP TABLE IF EXISTS webfront_varsplic')
    cur.execute('DROP TABLE IF EXISTS webfront_protein')
    cur.execute('DROP TABLE IF EXISTS webfront_database')
    cur.execute('DROP TABLE IF EXISTS webfront_proteome')
    cur.execute('DROP TABLE IF EXISTS webfront_taxonomy')

    cur.execute(
        """
        CREATE TABLE webfront_database
        (
            name VARCHAR(10) NOT NULL PRIMARY KEY,
            name_long VARCHAR(25) NOT NULL,
            description LONGTEXT,
            type ENUM('protein', 'entry', 'other') NOT NULL,
            version VARCHAR(20),
            release_date DATETIME,
            prev_version VARCHAR(20),
            prev_release_date DATETIME
        ) CHARSET=utf8 DEFAULT COLLATE=utf8_unicode_ci
        """
    )

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
            literature LONGTEXT,
            hierarchy LONGTEXT,
            cross_references LONGTEXT,
            overlaps_with LONGTEXT DEFAULT NULL,
            is_featured TINYINT NOT NULL DEFAULT 0,
            is_alive TINYINT NOT NULL DEFAULT 1,
            entry_date DATETIME NOT NULL,
            deletion_date DATETIME,
            counts LONGTEXT DEFAULT NULL,
            CONSTRAINT fk_entry_entry
              FOREIGN KEY (integrated_id)
              REFERENCES webfront_entry (accession),
            CONSTRAINT fk_entry_database
              FOREIGN KEY (source_database)
              REFERENCES webfront_database (name)
        ) CHARSET=utf8 DEFAULT COLLATE=utf8_unicode_ci
        """
    )

    cur.execute(
        """
        CREATE TABLE webfront_entryannotation
        (
            annotation_id VARCHAR(255) PRIMARY KEY NOT NULL,
            accession_id VARCHAR(25) NOT NULL,
            type VARCHAR(32) NOT NULL,
            value LONGBLOB NOT NULL,
            mime_type VARCHAR(32) NOT NULL,
            CONSTRAINT fk_entryannotation_entry
              FOREIGN KEY (accession_id)
              REFERENCES webfront_entry (accession)
        ) CHARSET=utf8 DEFAULT COLLATE=utf8_unicode_ci
        """
    )

    cur.execute(
        """
        CREATE TABLE webfront_taxonomy
        (
            accession VARCHAR(20) PRIMARY KEY NOT NULL,
            scientific_name VARCHAR(255) NOT NULL,
            full_name VARCHAR(512) NOT NULL,
            lineage LONGTEXT NOT NULL,
            parent_id VARCHAR(20),
            rank VARCHAR(20) NOT NULL,
            children LONGTEXT NOT NULL,
            counts LONGTEXT DEFAULT NULL,
            CONSTRAINT fk_taxonomy_taxonomy
              FOREIGN KEY (parent_id)
              REFERENCES webfront_taxonomy (accession)
        ) CHARSET=utf8 DEFAULT COLLATE=utf8_unicode_ci
        """
    )

    cur.execute(
        """
        CREATE TABLE webfront_proteome
        (
            accession VARCHAR(20) PRIMARY KEY NOT NULL,
            name VARCHAR(215) NOT NULL,
            is_reference TINYINT NOT NULL,
            strain VARCHAR(512),
            assembly VARCHAR(512),
            taxonomy_id VARCHAR(20) NOT NULL,
            counts LONGTEXT DEFAULT NULL,
            CONSTRAINT fk_proteome_taxonomy
              FOREIGN KEY (taxonomy_id)
              REFERENCES webfront_taxonomy (accession)
        ) CHARSET=utf8 DEFAULT COLLATE=utf8_unicode_ci
        """
    )

    cur.execute(
        """
        CREATE TABLE webfront_protein
        (
            accession VARCHAR(15) PRIMARY KEY NOT NULL,
            identifier VARCHAR(16) NOT NULL,
            organism LONGTEXT NOT NULL,
            name VARCHAR(255) NOT NULL,
            other_names LONGTEXT NOT NULL,
            description LONGTEXT NOT NULL,
            sequence LONGTEXT NOT NULL,
            length INT(11) NOT NULL,
            proteome VARCHAR(20),
            gene VARCHAR(70),
            go_terms LONGTEXT NOT NULL,
            evidence_code INT(11) NOT NULL,
            source_database VARCHAR(10) NOT NULL,
            residues LONGTEXT NOT NULL,
            is_fragment TINYINT NOT NULL,
            structure LONGTEXT NOT NULL,
            tax_id VARCHAR(20) NOT NULL,
            extra_features LONGTEXT NOT NULL,
            ida_id VARCHAR(40) DEFAULT NULL,
            ida TEXT DEFAULT NULL,
            counts LONGTEXT DEFAULT NULL,
            CONSTRAINT fk_protein_taxonomy
              FOREIGN KEY (tax_id)
              REFERENCES webfront_taxonomy (accession),
            CONSTRAINT fk_protein_database
              FOREIGN KEY (source_database)
              REFERENCES webfront_database (name),
            CONSTRAINT fk_protein_proteome
              FOREIGN KEY (proteome)
              REFERENCES webfront_proteome (accession)
        ) CHARSET=utf8 DEFAULT COLLATE=utf8_unicode_ci
        """
    )

    cur.execute(
        """
        CREATE TABLE webfront_varsplic
        (
            accession VARCHAR(20) PRIMARY KEY NOT NULL,
            protein_acc VARCHAR(15) NOT NULL,
            length INT(11) NOT NULL,
            sequence LONGTEXT NOT NULL,
            features LONGTEXT NOT NULL
        ) CHARSET=utf8 DEFAULT COLLATE=utf8_unicode_ci
        """
    )

    cur.execute(
        """
        CREATE TABLE webfront_structure
        (
            accession VARCHAR(4) PRIMARY KEY NOT NULL,
            name VARCHAR(512) NOT NULL,
            short_name VARCHAR(32) DEFAULT NULL,
            other_names LONGTEXT DEFAULT NULL,
            source_database VARCHAR(10) NOT NULL,
            experiment_type VARCHAR(16) NOT NULL,
            release_date DATETIME NOT NULL,
            resolution FLOAT DEFAULT NULL,
            literature LONGTEXT NOT NULL,
            chains LONGTEXT NOT NULL,
            proteins LONGTEXT NOT NULL,
            secondary_structures LONGTEXT NOT NULL,
            counts LONGTEXT DEFAULT NULL,
            CONSTRAINT fk_structure_database
              FOREIGN KEY (source_database)
              REFERENCES webfront_database (name)
        ) CHARSET=utf8 DEFAULT COLLATE=utf8_unicode_ci
        """
    )

    cur.execute(
        """
        CREATE TABLE webfront_set
        (
            accession VARCHAR(20) PRIMARY KEY NOT NULL,
            name VARCHAR(100),
            description TEXT,
            source_database VARCHAR(10) NOT NULL,
            is_set TINYINT NOT NULL,
            relationships LONGTEXT NOT NULL,
            integrated LONGTEXT DEFAULT NULL,
            counts LONGTEXT DEFAULT NULL,
            CONSTRAINT fk_set_database
              FOREIGN KEY (source_database)
              REFERENCES webfront_database (name)
        ) CHARSET=utf8 DEFAULT COLLATE=utf8_unicode_ci
        """
    )

    cur.execute(
        """
        CREATE TABLE webfront_alignment
        (
            id INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
            set_acc VARCHAR(20) NOT NULL,
            entry_acc VARCHAR(25) NOT NULL,
            target_acc VARCHAR(25) NOT NULL,
            target_set_acc VARCHAR(20),
            score DOUBLE NOT NULL,
            seq_length MEDIUMINT NOT NULL,
            domains TEXT NOT NULL,
            CONSTRAINT fk_alignment_set
              FOREIGN KEY (set_acc)
              REFERENCES webfront_set (accession)
        ) CHARSET=utf8 DEFAULT COLLATE=utf8_unicode_ci
        """
    )

    cur.execute(
        """
        CREATE TABLE IF NOT EXISTS webfront_release_note
        (
            version VARCHAR(20) PRIMARY KEY NOT NULL,
            release_date DATETIME NOT NULL,
            content LONGTEXT NOT NULL
        ) CHARSET=utf8 DEFAULT COLLATE=utf8_unicode_ci
        """
    )

    cur.close()
    con.close()


def reduce(src: dict) -> dict:
    dst = {}
    for k, v in src.items():
        if isinstance(v, dict):
            dst[k] = reduce(v)
        elif isinstance(v, set):
            dst[k] = len(v)
        else:
            dst[k] = v

    return dst


def update_counts(uri: str, src_entries: str, src_proteomes: str,
                  src_structures: str, src_taxa: str, tmpdir=None):
    taxonomy.update_counts(uri, src_taxa, tmpdir)
    proteome.update_counts(uri, src_proteomes)
    structure.update_counts(uri, src_structures)
    entry.update_counts(uri, src_entries, tmpdir)


from . import database, entry, protein, proteome, relnote, structure, taxonomy
