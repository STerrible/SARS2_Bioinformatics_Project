DEFAULT_REFERENCE_ACCESSION = "NC_045512"

METADATA_SHEET = "Metadata_counts"
CDS_SHEET = "CDS_features"
SEQUENCES_SHEET = "Sequences"
MUTATIONS_SHEET = "Mutations"
QC_SUMMARY_SHEET = "QC_summary"

MUTATION_COLUMNS = [
    "sample_accession",
    "reference_accession",
    "position",
    "ref_base",
    "sample_base",
    "mutation",
    "gene",
    "product",
    "protein_id",
    "is_in_cds",
    "requires_alignment",
    "comparison_note",
]
