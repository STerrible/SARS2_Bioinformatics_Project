from argparse import ArgumentParser
from pathlib import Path

from sars2_pipeline.config import (
    DEFAULT_INPUT_FASTA,
    DEFAULT_NEXTCLADE_DATASET,
    DEFAULT_NEXTCLADE_EXE,
    DEFAULT_NEXTCLADE_TSV,
)
from sars2_pipeline.excel_export import write_excel
from sars2_pipeline.genbank_parser import extract_genbank_tables, read_genbank_records
from sars2_pipeline.nextclade import read_nextclade_tsv, run_nextclade
from sars2_pipeline.nextclade_summary import build_nextclade_summary_tables
from sars2_pipeline.qc import build_qc_summary


def parse_genbank_to_excel(
    input_gb,
    output_xlsx,
    input_fasta=DEFAULT_INPUT_FASTA,
    nextclade_output_tsv=DEFAULT_NEXTCLADE_TSV,
    nextclade_exe=DEFAULT_NEXTCLADE_EXE,
    nextclade_dataset=DEFAULT_NEXTCLADE_DATASET,
):
    input_gb = Path(input_gb)
    output_xlsx = Path(output_xlsx)

    records = read_genbank_records(input_gb)
    metadata_rows, cds_rows, sequence_rows = extract_genbank_tables(records)
    qc_summary_rows = build_qc_summary(records, cds_rows)

    nextclade_output_tsv = run_nextclade(
        input_fasta,
        nextclade_output_tsv,
        nextclade_exe,
        nextclade_dataset,
    )
    nextclade_df = read_nextclade_tsv(nextclade_output_tsv)
    nextclade_qc_df, nextclade_mutations_df, nextclade_summary_df = build_nextclade_summary_tables(
        nextclade_df,
    )

    write_excel(
        output_xlsx,
        metadata_rows,
        cds_rows,
        sequence_rows,
        qc_summary_rows,
        nextclade_df,
        nextclade_qc_df,
        nextclade_mutations_df,
        nextclade_summary_df,
    )

    print(f"Done: {output_xlsx}")
    print(f"Records processed: {len(records)}")
    print(f"CDS processed: {len(cds_rows)}")
    print(f"Nextclade results: {nextclade_output_tsv}")


def parse_args():
    project_dir = Path(__file__).resolve().parents[2]
    parser = ArgumentParser(description="Convert GenBank records to an Excel workbook.")
    parser.add_argument(
        "--input",
        type=Path,
        default=project_dir / "data" / "raw" / "sequence.gb",
        help="Input GenBank file.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=project_dir / "results" / "genbank_table.xlsx",
        help="Output Excel workbook.",
    )
    parser.add_argument(
        "--fasta",
        type=Path,
        default=project_dir / DEFAULT_INPUT_FASTA,
        help="Input FASTA file for Nextclade.",
    )
    parser.add_argument(
        "--nextclade-exe",
        type=Path,
        default=Path(DEFAULT_NEXTCLADE_EXE),
        help="Path to nextclade.exe.",
    )
    parser.add_argument(
        "--nextclade-dataset",
        type=Path,
        default=Path(DEFAULT_NEXTCLADE_DATASET),
        help="Path to the local Nextclade dataset.",
    )
    parser.add_argument(
        "--nextclade-output",
        type=Path,
        default=project_dir / DEFAULT_NEXTCLADE_TSV,
        help="Output TSV file for Nextclade results.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    try:
        parse_genbank_to_excel(
            args.input,
            args.output,
            args.fasta,
            args.nextclade_output,
            args.nextclade_exe,
            args.nextclade_dataset,
        )
    except (PermissionError, ValueError) as exc:
        raise SystemExit(f"Error: {exc}") from None
