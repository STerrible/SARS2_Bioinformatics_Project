from argparse import ArgumentParser
from pathlib import Path

from sars2_pipeline.config import DEFAULT_REFERENCE_ACCESSION
from sars2_pipeline.excel_export import write_excel
from sars2_pipeline.genbank_parser import extract_genbank_tables, read_genbank_records
from sars2_pipeline.mutation_analysis import build_mutation_rows, require_reference_record
from sars2_pipeline.qc import build_qc_summary


def parse_genbank_to_excel(input_gb, output_xlsx, reference_accession=DEFAULT_REFERENCE_ACCESSION):
    input_gb = Path(input_gb)
    output_xlsx = Path(output_xlsx)

    records = read_genbank_records(input_gb)
    reference_record = require_reference_record(records, reference_accession, input_gb)
    metadata_rows, cds_rows, sequence_rows = extract_genbank_tables(records)

    mutation_rows, length_mismatch_records, not_comparable_records, mutation_count = build_mutation_rows(
        records,
        reference_record,
        reference_accession,
    )
    qc_summary_rows = build_qc_summary(
        records,
        cds_rows,
        reference_record,
        reference_accession,
        length_mismatch_records,
        not_comparable_records,
        mutation_count,
    )

    write_excel(output_xlsx, metadata_rows, cds_rows, sequence_rows, mutation_rows, qc_summary_rows)

    print(f"Done: {output_xlsx}")
    print(f"Records processed: {len(records)}")
    print(f"CDS processed: {len(cds_rows)}")
    print(f"Mutations found: {mutation_count}")
    print(f"Records requiring alignment: {len(not_comparable_records)}")


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
        "--reference-accession",
        default=DEFAULT_REFERENCE_ACCESSION,
        help="Reference accession for nucleotide mutation analysis.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    try:
        parse_genbank_to_excel(args.input, args.output, args.reference_accession)
    except ValueError as exc:
        raise SystemExit(f"Error: {exc}") from None
