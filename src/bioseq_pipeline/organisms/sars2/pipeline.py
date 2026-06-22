from argparse import ArgumentParser
from datetime import datetime
from pathlib import Path

from bioseq_pipeline.core.genbank import extract_genbank_tables, read_genbank_records
from bioseq_pipeline.core.qc import build_qc_summary
from bioseq_pipeline.organisms.sars2.config import (
    DEFAULT_INPUT_FASTA,
    DEFAULT_NEXTCLADE_ALIGNED_FASTA,
    DEFAULT_NEXTCLADE_DATASET,
    DEFAULT_NEXTCLADE_EXE,
    DEFAULT_NEXTCLADE_TSV,
)
from bioseq_pipeline.organisms.sars2.excel_export import write_excel
from bioseq_pipeline.organisms.sars2.nextclade_summary import (
    build_amino_acid_changes,
    build_amino_acid_changes_by_gene,
    build_country_mutations,
    build_country_summary,
    build_nextclade_gene_summary,
    build_nextclade_mutations,
    build_nextclade_qc,
    build_nextclade_summary,
    build_nextclade_top_mutations,
    build_run_metadata,
    enrich_mutations_with_metadata,
    enrich_nextclade_with_metadata,
)
from bioseq_pipeline.tools.nextclade import (
    get_nextclade_version,
    read_aligned_fasta,
    read_nextclade_dataset_info,
    read_nextclade_tsv,
    run_nextclade,
)


def project_dir():
    return Path(__file__).resolve().parents[4]


def parse_genbank_to_excel(
    input_gb,
    output_xlsx,
    input_fasta=DEFAULT_INPUT_FASTA,
    nextclade_output_tsv=DEFAULT_NEXTCLADE_TSV,
    nextclade_aligned_fasta=DEFAULT_NEXTCLADE_ALIGNED_FASTA,
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
        nextclade_aligned_fasta,
    )
    nextclade_df = read_nextclade_tsv(nextclade_output_tsv)
    aligned_fasta_df = read_aligned_fasta(nextclade_aligned_fasta)
    nextclade_dataset_info = read_nextclade_dataset_info(nextclade_dataset)
    nextclade_mutations_df = build_nextclade_mutations(nextclade_df)
    nextclade_summary_df = build_nextclade_summary(nextclade_df)
    enriched_nextclade_df = enrich_nextclade_with_metadata(nextclade_df, metadata_rows)
    enriched_mutations_df = enrich_mutations_with_metadata(nextclade_mutations_df, metadata_rows)
    nextclade_qc_df = build_nextclade_qc(enriched_nextclade_df)
    nextclade_mutations_df = enriched_mutations_df
    nextclade_gene_summary_df = build_nextclade_gene_summary(enriched_mutations_df)
    nextclade_top_mutations_df = build_nextclade_top_mutations(enriched_mutations_df, len(nextclade_df))
    country_summary_df = build_country_summary(enriched_nextclade_df)
    country_mutations_df = build_country_mutations(enriched_mutations_df, enriched_nextclade_df)
    amino_acid_changes_df = build_amino_acid_changes(enriched_mutations_df)
    amino_acid_changes_by_gene_df = build_amino_acid_changes_by_gene(enriched_mutations_df)
    run_metadata_df = build_run_metadata(
        nextclade_df,
        metadata_rows,
        cds_rows,
        input_gb,
        input_fasta,
        output_xlsx,
        nextclade_output_tsv,
        nextclade_aligned_fasta,
        nextclade_exe,
        nextclade_dataset,
        get_nextclade_version(nextclade_exe),
        datetime.now().astimezone().isoformat(timespec="seconds"),
        nextclade_dataset_info,
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
        nextclade_gene_summary_df,
        nextclade_top_mutations_df,
        country_summary_df,
        country_mutations_df,
        amino_acid_changes_df,
        amino_acid_changes_by_gene_df,
        aligned_fasta_df,
        run_metadata_df,
    )

    print(f"Done: {output_xlsx}")
    print(f"Records processed: {len(records)}")
    print(f"CDS processed: {len(cds_rows)}")
    print(f"Nextclade results: {nextclade_output_tsv}")
    print(f"Aligned FASTA: {nextclade_aligned_fasta}")


def add_sars2_arguments(parser):
    root_dir = project_dir()
    parser.add_argument(
        "--input",
        type=Path,
        default=root_dir / "data" / "raw" / "sequence.gb",
        help="Input GenBank file.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=root_dir / "results" / "genbank_table.xlsx",
        help="Output Excel workbook.",
    )
    parser.add_argument(
        "--fasta",
        type=Path,
        default=root_dir / DEFAULT_INPUT_FASTA,
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
        default=root_dir / DEFAULT_NEXTCLADE_TSV,
        help="Output TSV file for Nextclade results.",
    )
    parser.add_argument(
        "--nextclade-aligned-fasta",
        type=Path,
        default=root_dir / DEFAULT_NEXTCLADE_ALIGNED_FASTA,
        help="Output aligned FASTA file from Nextclade.",
    )


def parse_args():
    parser = ArgumentParser(description="Run the SARS-CoV-2 GenBank/Nextclade pipeline.")
    add_sars2_arguments(parser)
    return parser.parse_args()


def run_from_args(args):
    parse_genbank_to_excel(
        args.input,
        args.output,
        args.fasta,
        args.nextclade_output,
        args.nextclade_aligned_fasta,
        args.nextclade_exe,
        args.nextclade_dataset,
    )


def main():
    try:
        run_from_args(parse_args())
    except (PermissionError, ValueError) as exc:
        raise SystemExit(f"Error: {exc}") from None
