DEFAULT_REFERENCE_ACCESSION = "NC_045512"

DEFAULT_INPUT_FASTA = "data/raw/sequence.fasta"
DEFAULT_NEXTCLADE_EXE = r"C:\Games\Nextclade\nextclade.exe"
DEFAULT_NEXTCLADE_DATASET = r"C:\Games\Nextclade\sars-cov-2"
DEFAULT_NEXTCLADE_TSV = "results/nextclade.tsv"

METADATA_SHEET = "Metadata_counts"
CDS_SHEET = "CDS_features"
SEQUENCES_SHEET = "Sequences"
MUTATIONS_SHEET = "Mutations"
QC_SUMMARY_SHEET = "QC_summary"
NEXTCLADE_SHEET = "Nextclade_results"

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
