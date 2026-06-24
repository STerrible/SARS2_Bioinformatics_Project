import argparse
import gzip
import re
import shutil
from datetime import datetime
from pathlib import Path
from urllib.error import URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

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
    DEFAULT_INPUT_MANIFEST_TSV,
    DEFAULT_INPUT_MANIFEST_XLSX,
    DEFAULT_INPUT_VALIDATION_REPORT,
    DEFAULT_OUTPUT_XLSX,
    DEFAULT_PREPARED_FASTA_DIR,
    DEFAULT_REFERENCE_ACCESSION,
    DEFAULT_REFERENCE_DIR,
    TB_INPUT_MANIFEST_SHEET,
    TB_INPUT_VALIDATION_SHEET,
    TB_ASSEMBLY_METADATA_SHEET,
    TB_METADATA_SHEET,
    TB_RUN_METADATA_SHEET,
)


USER_AGENT = "SARS2_Bioinformatics_Project/1.0"
MIN_EXPECTED_TB_GENOME_LENGTH = 4_000_000
MAX_EXPECTED_TB_GENOME_LENGTH = 4_800_000


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


def find_assembly_dirs(assembly_dir):
    assembly_dir = Path(assembly_dir)
    if not assembly_dir.exists():
        raise ValueError(f"Assembly directory was not found: {assembly_dir}")

    assembly_dirs = sorted(path for path in assembly_dir.iterdir() if path.is_dir() and path.name.startswith(("GCF_", "GCA_")))
    if not assembly_dirs:
        raise ValueError(f"No assembly directories were found in: {assembly_dir}")
    return assembly_dirs


def single_file(directory, pattern):
    files = sorted(Path(directory).glob(pattern))
    return files[0] if files else None


def sample_id_from_accession(accession):
    sample_id = str(accession).replace(".", "_")
    sample_id = re.sub(r"[^A-Za-z0-9_.-]+", "_", sample_id).strip("_")
    return sample_id or "sample"


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


def metadata_by_assembly(metadata_df):
    if metadata_df.empty or "assembly_accession" not in metadata_df.columns:
        return {}
    return {
        str(row["assembly_accession"]): row.to_dict()
        for _, row in metadata_df.iterrows()
    }


def gc_percent(row):
    g = pd.to_numeric(pd.Series([row.get("G_%", "")]), errors="coerce").iloc[0]
    c = pd.to_numeric(pd.Series([row.get("C_%", "")]), errors="coerce").iloc[0]
    if pd.isna(g) or pd.isna(c):
        return ""
    return round(float(g) + float(c), 2)


def unpack_fasta_gz(source_path, destination_path, overwrite=False):
    source_path = Path(source_path)
    destination_path = Path(destination_path)
    destination_path.parent.mkdir(parents=True, exist_ok=True)

    if destination_path.exists() and destination_path.stat().st_size > 0 and not overwrite:
        return "skipped"

    partial_path = destination_path.with_name(destination_path.name + ".part")
    with gzip.open(source_path, "rb") as source, partial_path.open("wb") as destination:
        shutil.copyfileobj(source, destination)
    partial_path.replace(destination_path)
    return "written"


def efetch_reference(accession, rettype, output_path):
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    query = urlencode(
        {
            "db": "nuccore",
            "id": accession,
            "rettype": rettype,
            "retmode": "text",
        }
    )
    request = Request(
        f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?{query}",
        headers={"User-Agent": USER_AGENT},
    )
    with urlopen(request, timeout=90) as response:
        data = response.read()
    output_path.write_bytes(data)


def prepare_reference(reference_dir, reference_accession, download_reference=True):
    reference_dir = Path(reference_dir)
    reference_dir.mkdir(parents=True, exist_ok=True)

    fasta_path = reference_dir / f"{reference_accession}.fasta"
    genbank_path = reference_dir / f"{reference_accession}.gb"
    statuses = []

    for path, rettype in [(fasta_path, "fasta"), (genbank_path, "gbwithparts")]:
        if path.exists() and path.stat().st_size > 0:
            statuses.append((path, "exists"))
            continue
        if not download_reference:
            statuses.append((path, "missing"))
            continue
        try:
            efetch_reference(reference_accession, rettype, path)
            statuses.append((path, "downloaded"))
        except URLError as exc:
            statuses.append((path, f"download_failed: {exc.reason}"))

    return {
        "reference_accession": reference_accession,
        "reference_dir": reference_dir,
        "reference_fasta": fasta_path,
        "reference_genbank": genbank_path,
        "statuses": statuses,
    }


def build_manifest_row(assembly_dir, metadata_row, prepared_fasta_dir, unpack_fasta=True, overwrite=False):
    assembly_dir = Path(assembly_dir)
    assembly_accession = assembly_dir.name
    sample_id = sample_id_from_accession(assembly_accession)
    fna_path = single_file(assembly_dir, "*_genomic.fna.gz")
    gbff_path = single_file(assembly_dir, "*_genomic.gbff.gz")
    report_path = single_file(assembly_dir, "*_assembly_report.txt")
    prepared_fasta_path = Path(prepared_fasta_dir) / f"{sample_id}.fasta"

    warnings = []
    if fna_path is None:
        warnings.append("missing_fna")
    if gbff_path is None:
        warnings.append("missing_gbff")
    if report_path is None:
        warnings.append("missing_assembly_report")

    fasta_status = "not_requested"
    if unpack_fasta and fna_path is not None:
        fasta_status = unpack_fasta_gz(fna_path, prepared_fasta_path, overwrite=overwrite)

    organism = str(metadata_row.get("organism", ""))
    species_name = str(metadata_row.get("species_name", ""))
    assembly_organism = str(metadata_row.get("assembly_organism", ""))
    if "tuberculosis" not in f"{organism} {species_name} {assembly_organism}".lower():
        warnings.append("organism_not_confirmed_as_tuberculosis")

    length = pd.to_numeric(pd.Series([metadata_row.get("length", "")]), errors="coerce").iloc[0]
    if pd.isna(length):
        warnings.append("missing_length")
    elif length < MIN_EXPECTED_TB_GENOME_LENGTH or length > MAX_EXPECTED_TB_GENOME_LENGTH:
        warnings.append("unexpected_genome_length")

    return {
        "sample_id": sample_id,
        "assembly_accession": assembly_accession,
        "input_type": "RefSeq assembly FASTA",
        "raw_reads_used": "no",
        "fasta_gz_path": str(fna_path or ""),
        "prepared_fasta_path": str(prepared_fasta_path if unpack_fasta and fna_path is not None else ""),
        "fasta_status": fasta_status,
        "gbff_path": str(gbff_path or ""),
        "assembly_report_path": str(report_path or ""),
        "assembly_name": metadata_row.get("assembly_name", ""),
        "assembly_status": metadata_row.get("assembly_status", ""),
        "refseq_category": metadata_row.get("refseq_category", ""),
        "biosample": metadata_row.get("biosample", ""),
        "bioproject": metadata_row.get("bioproject", ""),
        "species_name": metadata_row.get("species_name", ""),
        "organism": metadata_row.get("organism", ""),
        "country": metadata_row.get("country", ""),
        "geo_loc_name": metadata_row.get("geo_loc_name", ""),
        "collection_date": metadata_row.get("collection_date", ""),
        "host": metadata_row.get("host", ""),
        "isolate": metadata_row.get("isolate", ""),
        "strain": metadata_row.get("strain", ""),
        "length": metadata_row.get("length", ""),
        "N_count": metadata_row.get("N_count", ""),
        "GC_percent": gc_percent(metadata_row),
        "warnings": ";".join(warnings),
    }


def build_input_manifest(
    assembly_dir,
    prepared_fasta_dir,
    unpack_fasta=True,
    overwrite=False,
):
    metadata_df = build_tb_metadata_counts(assembly_dir)
    metadata_lookup = metadata_by_assembly(metadata_df)
    rows = []
    for directory in find_assembly_dirs(assembly_dir):
        metadata_row = metadata_lookup.get(directory.name, {})
        rows.append(
            build_manifest_row(
                directory,
                metadata_row,
                prepared_fasta_dir,
                unpack_fasta=unpack_fasta,
                overwrite=overwrite,
            )
        )
    return pd.DataFrame(rows)


def validation_rows(manifest_df, reference_info):
    rows = [
        ("input_type", "RefSeq assemblies"),
        ("raw_reads_used", "no"),
        ("raw_reads_note", "Raw FASTQ reads are intentionally not used to keep the training workflow reproducible on a laptop."),
        ("assemblies", len(manifest_df)),
        ("missing_fna", int(manifest_df["fasta_gz_path"].eq("").sum())),
        ("missing_gbff", int(manifest_df["gbff_path"].eq("").sum())),
        ("missing_assembly_report", int(manifest_df["assembly_report_path"].eq("").sum())),
        ("prepared_fasta_files", int(manifest_df["prepared_fasta_path"].ne("").sum())),
        ("samples_with_warnings", int(manifest_df["warnings"].astype(str).str.strip().ne("").sum())),
        ("reference_accession", reference_info["reference_accession"]),
        ("reference_fasta_exists", reference_info["reference_fasta"].exists()),
        ("reference_genbank_exists", reference_info["reference_genbank"].exists()),
    ]
    for path, status in reference_info["statuses"]:
        rows.append((f"reference_{path.name}", status))
    return pd.DataFrame(rows, columns=["metric", "value"])


def write_validation_report(output_path, validation_df, manifest_df, reference_info):
    warning_rows = manifest_df[manifest_df["warnings"].astype(str).str.strip().ne("")]
    lines = [
        "# TB input validation",
        "",
        "This input layer uses RefSeq genome assemblies, not raw FASTQ reads.",
        "",
        "## Summary",
        "",
        "| metric | value |",
        "| --- | --- |",
    ]
    for _, row in validation_df.iterrows():
        lines.append(f"| {row['metric']} | {row['value']} |")

    lines.extend(
        [
            "",
            "## Reference",
            "",
            f"- accession: `{reference_info['reference_accession']}`",
            f"- FASTA: `{reference_info['reference_fasta']}`",
            f"- GenBank: `{reference_info['reference_genbank']}`",
            "",
            "## Warnings",
            "",
        ]
    )
    if warning_rows.empty:
        lines.append("- None")
    else:
        for _, row in warning_rows.iterrows():
            lines.append(f"- `{row['sample_id']}`: {row['warnings']}")

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def prepare_inputs(
    assembly_dir,
    manifest_tsv,
    manifest_xlsx,
    validation_report,
    prepared_fasta_dir,
    reference_dir,
    reference_accession=DEFAULT_REFERENCE_ACCESSION,
    unpack_fasta=True,
    overwrite=False,
    download_reference=True,
):
    manifest_tsv = Path(manifest_tsv)
    manifest_xlsx = Path(manifest_xlsx)
    validation_report = Path(validation_report)

    reference_info = prepare_reference(reference_dir, reference_accession, download_reference=download_reference)
    manifest_df = build_input_manifest(
        assembly_dir,
        prepared_fasta_dir,
        unpack_fasta=unpack_fasta,
        overwrite=overwrite,
    )
    validation_df = validation_rows(manifest_df, reference_info)

    manifest_tsv.parent.mkdir(parents=True, exist_ok=True)
    manifest_df.to_csv(manifest_tsv, sep="\t", index=False)
    write_dataframes_to_excel(
        manifest_xlsx,
        [
            (TB_INPUT_MANIFEST_SHEET, manifest_df),
            (TB_INPUT_VALIDATION_SHEET, validation_df),
        ],
    )
    write_validation_report(validation_report, validation_df, manifest_df, reference_info)

    return {
        "manifest_tsv": manifest_tsv,
        "manifest_xlsx": manifest_xlsx,
        "validation_report": validation_report,
        "prepared_fasta_dir": Path(prepared_fasta_dir),
        "reference_fasta": reference_info["reference_fasta"],
        "reference_genbank": reference_info["reference_genbank"],
        "assemblies": len(manifest_df),
        "warnings": int(manifest_df["warnings"].astype(str).str.strip().ne("").sum()),
    }


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
    parser.set_defaults(tb_command="metadata")
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

    subparsers = parser.add_subparsers(dest="tb_command")
    metadata_parser = subparsers.add_parser(
        "metadata",
        help="Build TB metadata/counts workbook.",
    )
    metadata_parser.add_argument(
        "--assembly-dir",
        type=Path,
        default=DEFAULT_ASSEMBLY_DIR,
        help="Directory containing downloaded NCBI assembly folders with *_genomic.gbff.gz files.",
    )
    metadata_parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT_XLSX,
        help="Output Excel workbook.",
    )

    prepare_parser = subparsers.add_parser(
        "prepare-inputs",
        help="Prepare TB manifest, FASTA inputs, validation report, and H37Rv reference.",
    )
    prepare_parser.add_argument(
        "--assembly-dir",
        type=Path,
        default=DEFAULT_ASSEMBLY_DIR,
        help="Directory containing downloaded NCBI assembly folders.",
    )
    prepare_parser.add_argument(
        "--manifest-tsv",
        type=Path,
        default=DEFAULT_INPUT_MANIFEST_TSV,
        help="Output TSV manifest path.",
    )
    prepare_parser.add_argument(
        "--manifest-xlsx",
        type=Path,
        default=DEFAULT_INPUT_MANIFEST_XLSX,
        help="Output Excel manifest path.",
    )
    prepare_parser.add_argument(
        "--validation-report",
        type=Path,
        default=DEFAULT_INPUT_VALIDATION_REPORT,
        help="Output Markdown validation report path.",
    )
    prepare_parser.add_argument(
        "--prepared-fasta-dir",
        type=Path,
        default=DEFAULT_PREPARED_FASTA_DIR,
        help="Directory for unpacked FASTA files with stable sample names.",
    )
    prepare_parser.add_argument(
        "--reference-dir",
        type=Path,
        default=DEFAULT_REFERENCE_DIR,
        help="Directory for H37Rv reference FASTA and GenBank files.",
    )
    prepare_parser.add_argument(
        "--reference-accession",
        default=DEFAULT_REFERENCE_ACCESSION,
        help="Reference accession to download/check.",
    )
    prepare_parser.add_argument(
        "--no-unpack-fasta",
        action="store_true",
        help="Do not unpack *_genomic.fna.gz files into the prepared FASTA directory.",
    )
    prepare_parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite prepared FASTA files if they already exist.",
    )
    prepare_parser.add_argument(
        "--skip-reference-download",
        action="store_true",
        help="Only check reference files; do not download missing files from NCBI.",
    )


def parse_args():
    parser = argparse.ArgumentParser(description="Build tuberculosis metadata/counts workbook.")
    add_tuberculosis_arguments(parser)
    return parser.parse_args()


def run_from_args(args):
    if getattr(args, "tb_command", "metadata") == "prepare-inputs":
        result = prepare_inputs(
            args.assembly_dir,
            args.manifest_tsv,
            args.manifest_xlsx,
            args.validation_report,
            args.prepared_fasta_dir,
            args.reference_dir,
            reference_accession=args.reference_accession,
            unpack_fasta=not args.no_unpack_fasta,
            overwrite=args.overwrite,
            download_reference=not args.skip_reference_download,
        )
        print(f"Manifest TSV: {result['manifest_tsv']}")
        print(f"Manifest Excel: {result['manifest_xlsx']}")
        print(f"Validation report: {result['validation_report']}")
        print(f"Prepared FASTA directory: {result['prepared_fasta_dir']}")
        print(f"Reference FASTA: {result['reference_fasta']}")
        print(f"Reference GenBank: {result['reference_genbank']}")
        print(f"Assemblies prepared: {result['assemblies']}")
        print(f"Samples with warnings: {result['warnings']}")
        return

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
