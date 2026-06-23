from pathlib import Path


PROJECT_DIR = Path(__file__).resolve().parents[4]
DEFAULT_ASSEMBLY_DIR = PROJECT_DIR / "data" / "raw" / "tuberculosis_data" / "ncbi_assemblies"
DEFAULT_OUTPUT_XLSX = PROJECT_DIR / "results" / "tuberculosis_data" / "tb_metadata.xlsx"

TB_METADATA_SHEET = "TB_Metadata_counts"
TB_ASSEMBLY_METADATA_SHEET = "TB_Assembly_metadata"
TB_RUN_METADATA_SHEET = "TB_Run_Metadata"
