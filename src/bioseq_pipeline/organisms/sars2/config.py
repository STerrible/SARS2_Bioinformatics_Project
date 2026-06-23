import os


DEFAULT_INPUT_FASTA = "data/raw/covid_data/sequence.fasta"
DEFAULT_NEXTCLADE_EXE = os.environ.get("NEXTCLADE_EXE", r"C:\Games\Nextclade\nextclade.exe")
DEFAULT_NEXTCLADE_DATASET = os.environ.get("NEXTCLADE_DATASET", r"C:\Games\Nextclade\sars-cov-2")
DEFAULT_NEXTCLADE_TSV = "results/covid_data/nextclade.tsv"
DEFAULT_NEXTCLADE_ALIGNED_FASTA = "results/covid_data/aligned.fasta"

METADATA_SHEET = "Metadata_counts"
CDS_SHEET = "CDS_features"
SEQUENCES_SHEET = "Sequences"
QC_SUMMARY_SHEET = "QC_summary"
NEXTCLADE_SHEET = "Nextclade_results"
NEXTCLADE_QC_SHEET = "Nextclade_QC"
NEXTCLADE_MUTATIONS_SHEET = "Nextclade_Mutations"
NEXTCLADE_SUMMARY_SHEET = "Nextclade_Summary"
NEXTCLADE_GENE_SUMMARY_SHEET = "Nextclade_Gene_Summary"
NEXTCLADE_TOP_MUTATIONS_SHEET = "Nextclade_Top_Mutations"
COUNTRY_SUMMARY_SHEET = "Country_Summary"
COUNTRY_MUTATIONS_SHEET = "Country_Mutations"
AMINO_ACID_CHANGES_SHEET = "Amino_Acid_Changes"
AMINO_ACID_CHANGES_BY_GENE_SHEET = "Amino_Acid_Changes_By_Gene"
ALIGNED_FASTA_SHEET = "aligned_fasta"
RUN_METADATA_SHEET = "Run_Metadata"
