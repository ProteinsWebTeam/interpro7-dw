IDA_INDEX = "ida"
REL_DEFAULT_INDEX = "others"  # for docs without an entry (no mem DB)

EXTENSION = ".dat"  # Extension given to files after they have been written
LOAD_SUFFIX = ".load"  # For sentinel file while docs are written
DONE_SUFFIX = ".done"  # For sentinel file once all files have been written

IDA_ALIAS = "ida"
REL_ALIAS = "rel"
STAGING_ALIAS_SUFFIX = "_staging"
LIVE_ALIAS_SUFFIX = "_current"
PREVIOUS_ALIAS_SUFFIX = "_previous"

NUM_SHARDS = {
    "antifam": 1,
    "cathgene3d": 10,
    "hamap": 1,
    "ida": 1,
    "interpro": 10,
    "ncbifam": 3,
    "panther": 5,
    "pfam": 10,
    "pirsf": 1,
    "prints": 1,
    "profile": 3,
    "prosite": 2,
    "sfld": 1,
    "smart": 3,
    REL_DEFAULT_INDEX: 1
}
DEFAULT_SHARDS = 5

IDA_BODY = {
    "mappings": {
        "properties": {
            "ida_id": {"type": "keyword"},
            "ida": {"type": "keyword"},
            "representative": {
                "type": "object",
                "enabled": False
            },
            "counts": {"type": "integer"}
        }
    }
}
REL_BODY = {
    "settings": {
        "analysis": {
            "analyzer": {
                "autocomplete": {
                    "tokenizer": "autocomplete",
                    "filter": [
                        "lowercase"
                    ]
                }
            },
            "tokenizer": {
                "autocomplete": {
                    "type": "edge_ngram",
                    "min_gram": 2,
                    "max_gram": 20,
                    "token_chars": [
                        "letter",
                        "digit"
                    ]
                }
            }
        }
    },
    "mappings": {
        "properties": {
            # Protein
            "protein_acc": {"type": "keyword"},
            "protein_length": {"type": "long"},
            "protein_is_fragment": {"type": "keyword"},
            "protein_af_score": {"type": "float"},
            "protein_db": {"type": "keyword"},
            "protein_gene": {"type": "keyword"},
            "text_protein": {"type": "text", "analyzer": "autocomplete"},

            # Domain architecture
            "ida_id": {"type": "keyword"},
            "ida": {"type": "keyword"},

            # Taxonomy
            "tax_id": {"type": "long"},
            "tax_name": {"type": "keyword"},
            "tax_lineage": {"type": "keyword"},
            "tax_rank": {"type": "keyword"},
            "text_taxonomy": {"type": "text", "analyzer": "autocomplete"},

            # Proteome
            "proteome_acc": {"type": "keyword"},
            "proteome_name": {"type": "keyword"},
            "proteome_is_reference": {"type": "keyword"},
            "text_proteome": {"type": "text", "analyzer": "autocomplete"},

            # Structure
            "structure_acc": {"type": "keyword"},
            "structure_resolution": {"type": "float"},
            "structure_date": {"type": "date"},
            "structure_evidence": {"type": "keyword"},
            "protein_structure": {"type": "object", "enabled": False},
            "text_structure": {"type": "text", "analyzer": "autocomplete"},

            # Chain
            "structure_chain_acc": {"type": "text", "analyzer": "keyword"},
            "structure_chain": {"type": "text", "analyzer": "keyword", "fielddata": True},
            "structure_protein_length": {"type": "long"},
            "structure_protein_locations": {"type": "object", "enabled": False},

            # Entry
            "entry_acc": {"type": "keyword"},
            "entry_db": {"type": "keyword"},
            "entry_type": {"type": "keyword"},
            "entry_date": {"type": "date"},
            "entry_protein_locations": {"type": "object", "enabled": False},
            "entry_structure_locations": {"type": "object", "enabled": False},
            "entry_go_terms": {"type": "keyword"},
            "entry_integrated": {"type": "keyword"},
            "text_entry": {"type": "text", "analyzer": "autocomplete"},

            # Clan/set
            "set_acc": {"type": "keyword"},
            "set_db": {"type": "keyword"},
            "set_integrated": {"type": "keyword"},         # todo: remove?
            "text_set": {"type": "text", "analyzer": "autocomplete"},
        }
    }
}
