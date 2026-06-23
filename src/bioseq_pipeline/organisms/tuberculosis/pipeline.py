import argparse
import gzip
from datetime import datetime
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from bioseq_pipeline.core.excel import write_dataframes_to_excel
from bioseq_pipeline.core.genbank import (
    count_nucleotides,
    extract_country,
    get_accession,
    get_source_feature,
    get_source_qualifier,
)
from bioseq_pipeline.organisms.tuberculosis.config import (
    DEFAULT_ASSEMBLY_DIR,
    DEFAULT_OUTPUT_XLSX,
    TB_ASSEMBLY_METADATA_SHEET,
    TB_METADATA_SHEET,
    TB_RUN_METADATA_SHEET,
)


def read_assembly_metadata(assembly_dir):
    metadata_path = Path(assembly_dir) / "assembly_metadata.tsv"
    if not metadata_path.exists():
        return {}

    metadata = pd.read_csv(metadata_path, sep="\t", keep_default_na=False, dtype=str)
    if "accession" not in metadata.columns:
        return {}

    return {
        str(row["accession"]): row.to_dict()
        for _, row in metadata.iterrows()
    }


def find_gbff_files(assembly_dir):
    assembly_dir = Path(assembly_dir)
    if not assembly_dir.exists():
        raise ValueError(f"Assembly directory was not found: {assembly_dir}")

    files = sorted(assembly_dir.glob("*/*_genomic.gbff.gz"))
    if not files:
        raise ValueError(f"No *_genomic.gbff.gz files were found in: {assembly_dir}")
    return files


def read_genbank_records_gzip(gbff_path):
    with gzip.open(gbff_path, "rt", encoding="utf-8") as handle:
        yield from SeqIO.parse(handle, "genbank")


def source_db_xrefs(record):
    source = get_source_feature(record)
    if source is None:
        return ""
    return ";".join(source.qualifiers.get("db_xref", []))


def metadata_row_from_record(record, assembly_accession, gbff_path, assembly_metadata):
    country, geo_loc_name = extract_country(record)
    counts = count_nucleotides(record.seq)

    return {
        "assembly_accession": assembly_accession,
        "assembly_name": assembly_metadata.get("assembly_name", ""),
        "assembly_status": assembly_metadata.get("assembly_status", ""),
        "refseq_category": assembly_metadata.get("refseq_category", ""),
        "biosample": assembly_metadata.get("biosample", ""),
        "bioproject": assembly_metadata.get("bioproject", ""),
        "assembly_uid": assembly_metadata.get("uid", ""),
        "taxid": assembly_metadata.get("taxid", ""),
        "species_name": assembly_metadata.get("species_name", ""),
        "assembly_organism": assembly_metadata.get("organism", ""),
        "gbff_path": str(gbff_path),
        "accession": get_accession(record),
        "accession_version": record.id,
        "description": record.description,
        "organism": record.annotations.get("organism", ""),
        "molecule_type": record.annotations.get("molecule_type", ""),
        "country": country,
        "geo_loc_name": geo_loc_name,
        "collection_date": get_source_qualifier(record, "collection_date"),
        "host": get_source_qualifier(record, "host"),
        "isolate": get_source_qualifier(record, "isolate"),
        "strain": get_source_qualifier(record, "strain"),
        "isolation_source": get_source_qualifier(record, "isolation_source"),
        "lat_lon": get_source_qualifier(record, "lat_lon"),
        "source_db_xref": source_db_xrefs(record),
        **counts,
    }


def build_tb_metadata_counts(assembly_dir):
    assembly_metadata_by_accession = read_assembly_metadata(assembly_dir)
    rows = []
    files = find_gbff_files(assembly_dir)

    for gbff_path in files:
        assembly_accession = gbff_path.parent.name
        assembly_metadata = assembly_metadata_by_accession.get(assembly_accession, {})
        record_count = 0
        for record in read_genbank_records_gzip(gbff_path):
            record_count += 1
            rows.append(metadata_row_from_record(record, assembly_accession, gbff_path, assembly_metadata))
        if record_count == 0:
            rows.append(
                {
                    "assembly_accession": assembly_accession,
                    "gbff_path": str(gbff_path),
                    "accession": "",
                    "accession_version": "",
                    "description": "No GenBank records found",
                }
            )

    return pd.DataFrame(rows)


def build_run_metadata(assembly_dir, output_xlsx, metadata_df):
    return pd.DataFrame(
        [
            ("generated_at", datetime.now().astimezone().isoformat(timespec="seconds")),
            ("assembly_dir", str(Path(assembly_dir))),
            ("output_xlsx", str(Path(output_xlsx))),
            ("rows", len(metadata_df)),
            ("assemblies", metadata_df["assembly_accession"].nunique() if "assembly_accession" in metadata_df else 0),
        ],
        columns=["metric", "value"],
    )


def write_tb_excel(assembly_dir, output_xlsx):
    output_xlsx = Path(output_xlsx)
    metadata_df = build_tb_metadata_counts(assembly_dir)
    assembly_metadata_path = Path(assembly_dir) / "assembly_metadata.tsv"
    assembly_metadata_df = (
        pd.read_csv(assembly_metadata_path, sep="\t", keep_default_na=False, dtype=str)
        if assembly_metadata_path.exists()
        else pd.DataFrame()
    )
    run_metadata_df = build_run_metadata(assembly_dir, output_xlsx, metadata_df)

    write_dataframes_to_excel(
        output_xlsx,
        [
            (TB_METADATA_SHEET, metadata_df),
            (TB_ASSEMBLY_METADATA_SHEET, assembly_metadata_df),
            (TB_RUN_METADATA_SHEET, run_metadata_df),
        ],
    )
    return {
        "output": output_xlsx,
        "rows": len(metadata_df),
        "assemblies": metadata_df["assembly_accession"].nunique() if "assembly_accession" in metadata_df else 0,
    }


def add_tuberculosis_arguments(parser):
    parser.add_argument(
        "--assembly-dir",
        type=Path,
        default=DEFAULT_ASSEMBLY_DIR,
        help="Directory containing downloaded NCBI assembly folders with *_genomic.gbff.gz files.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT_XLSX,
        help="Output Excel workbook.",
    )


def parse_args():
    parser = argparse.ArgumentParser(description="Build tuberculosis metadata/counts workbook.")
    add_tuberculosis_arguments(parser)
    return parser.parse_args()


def run_from_args(args):
    result = write_tb_excel(args.assembly_dir, args.output)
    print(f"Done: {result['output']}")
    print(f"Assemblies processed: {result['assemblies']}")
    print(f"Metadata rows: {result['rows']}")


def main():
    try:
        run_from_args(parse_args())
    except (PermissionError, ValueError) as exc:
        raise SystemExit(f"Error: {exc}") from None


if __name__ == "__main__":
    main()
