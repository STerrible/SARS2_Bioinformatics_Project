from pathlib import Path


PROJECT_DIR = Path(__file__).resolve().parents[4]
DEFAULT_ASSEMBLY_DIR = PROJECT_DIR / "data" / "raw" / "tuberculosis_data" / "ncbi_assemblies"
DEFAULT_REFERENCE_DIR = PROJECT_DIR / "data" / "raw" / "tuberculosis_data" / "reference"
DEFAULT_REFERENCE_ACCESSION = "NC_000962.3"
DEFAULT_OUTPUT_XLSX = PROJECT_DIR / "results" / "tuberculosis_data" / "tb_metadata.xlsx"
DEFAULT_RESULTS_DIR = PROJECT_DIR / "results" / "tuberculosis_data"
DEFAULT_INPUT_MANIFEST_TSV = DEFAULT_RESULTS_DIR / "tb_input_manifest.tsv"
DEFAULT_INPUT_MANIFEST_XLSX = DEFAULT_RESULTS_DIR / "tb_input_manifest.xlsx"
DEFAULT_INPUT_VALIDATION_REPORT = DEFAULT_RESULTS_DIR / "tb_input_validation.md"
DEFAULT_PREPARED_FASTA_DIR = DEFAULT_RESULTS_DIR / "inputs" / "fasta"

TB_METADATA_SHEET = "TB_Metadata_counts"
TB_ASSEMBLY_METADATA_SHEET = "TB_Assembly_metadata"
TB_INPUT_MANIFEST_SHEET = "TB_Input_manifest"
TB_INPUT_VALIDATION_SHEET = "TB_Input_validation"
TB_RUN_METADATA_SHEET = "TB_Run_Metadata"
