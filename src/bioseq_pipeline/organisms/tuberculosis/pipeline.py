import argparse
import gzip
import re
import shlex
import shutil
import subprocess
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
    DEFAULT_FASTTREE_LOG,
    DEFAULT_FASTTREE_README,
    DEFAULT_FASTTREE_TREE,
    DEFAULT_GENE_SUMMARY_TSV,
    DEFAULT_ITOL_COUNTRY_STRIP,
    DEFAULT_ITOL_README,
    DEFAULT_ITOL_TARGET_HEATMAP,
    DEFAULT_ITOL_VARIANT_BARS,
    DEFAULT_MUTATION_REPORT,
    DEFAULT_MUTATION_SUMMARY_XLSX,
    DEFAULT_PHYLOGENETICS_DIR,
    DEFAULT_REFERENCE_ACCESSION,
    DEFAULT_REFERENCE_DIR,
    DEFAULT_SAMPLE_SUMMARY_TSV,
    DEFAULT_SNIPPY_COMMANDS,
    DEFAULT_SNIPPY_CORE_ALN,
    DEFAULT_SNIPPY_CORE_NO_REFERENCE_ALN,
    DEFAULT_SNIPPY_CORE_TAB,
    DEFAULT_SNIPPY_CORE_TXT,
    DEFAULT_SNIPPY_CORE_VCF,
    DEFAULT_SNIPPY_DIR,
    DEFAULT_SNIPPY_MANIFEST,
    DEFAULT_SNIPPY_README,
    DEFAULT_SNIPPY_REFERENCE,
    DEFAULT_SNIPPY_RUNS_DIR,
    DEFAULT_SNIPPY_SUMMARY_REPORT,
    DEFAULT_SNIPPY_SUMMARY_XLSX,
    DEFAULT_TARGET_MATRIX_TSV,
    DEFAULT_TB_REPORT,
    DEFAULT_VARIANTS_TSV,
    TB_INPUT_MANIFEST_SHEET,
    TB_INPUT_VALIDATION_SHEET,
    TB_ASSEMBLY_METADATA_SHEET,
    TB_METADATA_SHEET,
    TB_RUN_METADATA_SHEET,
    TB_SNIPPY_CORE_STATS_SHEET,
    TB_SNIPPY_RUN_STATUS_SHEET,
    TB_SNIPPY_SUMMARY_SHEET,
)
from bioseq_pipeline.organisms.tuberculosis.variants import (
    export_itol_annotations,
    write_final_report,
    write_mutation_summary,
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


def path_to_wsl(path):
    text = str(path).strip()
    if not text:
        return ""

    normalized = text.replace("\\", "/")
    drive_match = re.match(r"^([A-Za-z]):/(.*)$", normalized)
    if drive_match:
        drive = drive_match.group(1).lower()
        rest = drive_match.group(2)
        return f"/mnt/{drive}/{rest}"
    return normalized


def local_path_from_possible_wsl(path):
    text = str(path).strip().replace("\\", "/")
    mount_match = re.match(r"^/mnt/([A-Za-z])/(.*)$", text)
    if mount_match:
        drive = mount_match.group(1).upper()
        rest = mount_match.group(2).replace("/", "\\")
        return Path(f"{drive}:\\{rest}")
    return Path(path)


def shell_quote(value):
    return shlex.quote(str(value))


def write_text_lf(path, text):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(text)


def read_input_manifest(manifest_path):
    manifest_path = Path(manifest_path)
    if not manifest_path.exists():
        raise ValueError(f"Input manifest was not found: {manifest_path}")

    manifest_df = pd.read_csv(manifest_path, sep="\t", keep_default_na=False, dtype=str)
    required_columns = {"sample_id", "prepared_fasta_path"}
    missing_columns = sorted(required_columns - set(manifest_df.columns))
    if missing_columns:
        raise ValueError(f"Input manifest is missing required columns: {', '.join(missing_columns)}")
    if "warnings" not in manifest_df.columns:
        manifest_df["warnings"] = ""
    return manifest_df


def build_snippy_plan_rows(manifest_df, include_warnings=False):
    rows = []
    skipped = []
    seen_sample_ids = set()

    for _, row in manifest_df.iterrows():
        sample_id = str(row.get("sample_id", "")).strip()
        prepared_fasta_path = str(row.get("prepared_fasta_path", "")).strip()
        warnings = str(row.get("warnings", "")).strip()
        skip_reasons = []

        if not sample_id:
            skip_reasons.append("missing_sample_id")
        elif sample_id in seen_sample_ids:
            skip_reasons.append("duplicate_sample_id")

        if not prepared_fasta_path:
            skip_reasons.append("missing_prepared_fasta_path")
        elif not local_path_from_possible_wsl(prepared_fasta_path).exists():
            skip_reasons.append("prepared_fasta_not_found")

        if warnings and not include_warnings:
            skip_reasons.append(f"manifest_warnings: {warnings}")

        if skip_reasons:
            skipped.append(
                {
                    "sample_id": sample_id or "(missing)",
                    "prepared_fasta_path": prepared_fasta_path,
                    "reason": "; ".join(skip_reasons),
                }
            )
            continue

        seen_sample_ids.add(sample_id)
        rows.append(
            {
                "sample_id": sample_id,
                "prepared_fasta_path": prepared_fasta_path,
                "prepared_fasta_wsl": path_to_wsl(prepared_fasta_path),
                "warnings": warnings,
            }
        )

    if not rows:
        raise ValueError("No valid samples were available for Snippy planning.")
    return rows, skipped


def write_snippy_multi(rows, output_path):
    output_path = Path(output_path)
    lines = [f"{row['sample_id']}\t{row['prepared_fasta_wsl']}" for row in rows]
    write_text_lf(output_path, "\n".join(lines) + "\n")


def write_snippy_commands(
    rows,
    commands_path,
    reference_path,
    snippy_dir,
    runs_dir,
    snippy_manifest,
    cpus=4,
    ram=8,
    force=False,
    cleanup=False,
):
    commands_path = Path(commands_path)

    extra_flags = []
    if force:
        extra_flags.append("--force")
    if cleanup:
        extra_flags.append("--cleanup")
    extra_options = f" {' '.join(extra_flags)}" if extra_flags else ""

    lines = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        "",
        "# Generated by: python main.py tuberculosis snippy-plan",
        "# Run inside WSL after activating the conda environment with Snippy.",
        f"PLAN_DIR={shell_quote(path_to_wsl(snippy_dir))}",
        f"RUNS_DIR={shell_quote(path_to_wsl(runs_dir))}",
        f"MULTI={shell_quote(path_to_wsl(snippy_manifest))}",
        f"REF={shell_quote(path_to_wsl(reference_path))}",
        f"SNIPPY_CPUS={int(cpus)}",
        f"SNIPPY_RAM={int(ram)}",
        "",
        'if ! command -v snippy >/dev/null 2>&1; then',
        '  echo "snippy was not found. Activate the conda environment first: conda activate snippy" >&2',
        "  exit 1",
        "fi",
        "",
        'mkdir -p "$PLAN_DIR" "$RUNS_DIR"',
        "snippy --check",
        "",
        "# These calls use --ctgs because this project currently uses assembly FASTA files, not raw FASTQ reads.",
    ]

    for row in rows:
        outdir = Path(runs_dir) / row["sample_id"]
        lines.append(
            " ".join(
                [
                    "snippy",
                    '--cpus "$SNIPPY_CPUS"',
                    '--ram "$SNIPPY_RAM"',
                    f"--outdir {shell_quote(path_to_wsl(outdir))}",
                    '--ref "$REF"',
                    f"--ctgs {shell_quote(row['prepared_fasta_wsl'])}{extra_options}",
                ]
            )
        )

    lines.extend(
        [
            "",
            "# Build a core SNP alignment after all sample runs finish.",
            'snippy-core --ref "$REF" --prefix "$PLAN_DIR/core" "$RUNS_DIR"/*',
            "",
            'echo "Done. Snippy outputs are in: $PLAN_DIR"',
        ]
    )
    write_text_lf(commands_path, "\n".join(lines) + "\n")


def write_snippy_readme(
    readme_path,
    rows,
    skipped,
    commands_path,
    snippy_manifest,
    reference_path,
    snippy_dir,
    runs_dir,
    cpus,
    ram,
):
    readme_path = Path(readme_path)
    skipped_lines = ["- None"] if not skipped else [
        f"- `{row['sample_id']}`: {row['reason']}"
        for row in skipped
    ]

    lines = [
        "# TB Snippy plan",
        "",
        "This directory contains a planned Snippy run for the tuberculosis assembly dataset.",
        "It does not mean Snippy has already been executed.",
        "",
        "## Generated files",
        "",
        f"- Snippy multi manifest: `{snippy_manifest}`",
        f"- Shell command script: `{commands_path}`",
        f"- Planned Snippy output directory: `{runs_dir}`",
        f"- Reference GenBank: `{reference_path}`",
        "",
        "## Planned run",
        "",
        f"- Samples: {len(rows)}",
        f"- CPU threads per sample: {cpus}",
        f"- RAM setting per sample: {ram} GB",
        "- Input mode: assembly FASTA via `snippy --ctgs`",
        "- Raw FASTQ reads: not used in this training workflow",
        "",
        "## How to run in WSL",
        "",
        "```bash",
        "source ~/miniforge3/etc/profile.d/conda.sh",
        "conda activate snippy",
        f"bash {path_to_wsl(commands_path)}",
        "```",
        "",
        "## Main outputs after the script finishes",
        "",
        f"- Per-sample Snippy folders: `{runs_dir}`",
        f"- Core SNP alignment prefix: `{Path(snippy_dir) / 'core'}`",
        "",
        "## Skipped samples",
        "",
        *skipped_lines,
    ]
    write_text_lf(readme_path, "\n".join(lines) + "\n")


def create_snippy_plan(
    manifest_path=DEFAULT_INPUT_MANIFEST_TSV,
    reference_path=DEFAULT_SNIPPY_REFERENCE,
    snippy_dir=DEFAULT_SNIPPY_DIR,
    snippy_manifest=None,
    commands_path=None,
    readme_path=None,
    runs_dir=None,
    cpus=4,
    ram=8,
    force=False,
    cleanup=False,
    include_warnings=False,
):
    if int(cpus) < 1:
        raise ValueError("--cpus must be at least 1")
    if int(ram) < 1:
        raise ValueError("--ram must be at least 1")

    snippy_dir = Path(snippy_dir)
    snippy_manifest = Path(snippy_manifest) if snippy_manifest else snippy_dir / DEFAULT_SNIPPY_MANIFEST.name
    commands_path = Path(commands_path) if commands_path else snippy_dir / DEFAULT_SNIPPY_COMMANDS.name
    readme_path = Path(readme_path) if readme_path else snippy_dir / DEFAULT_SNIPPY_README.name
    runs_dir = Path(runs_dir) if runs_dir else snippy_dir / DEFAULT_SNIPPY_RUNS_DIR.name
    reference_path = Path(reference_path)

    if not local_path_from_possible_wsl(reference_path).exists():
        raise ValueError(f"Reference GenBank was not found: {reference_path}")

    manifest_df = read_input_manifest(manifest_path)
    rows, skipped = build_snippy_plan_rows(manifest_df, include_warnings=include_warnings)

    snippy_dir.mkdir(parents=True, exist_ok=True)
    write_snippy_multi(rows, snippy_manifest)
    write_snippy_commands(
        rows,
        commands_path,
        reference_path,
        snippy_dir,
        runs_dir,
        snippy_manifest,
        cpus=cpus,
        ram=ram,
        force=force,
        cleanup=cleanup,
    )
    write_snippy_readme(
        readme_path,
        rows,
        skipped,
        commands_path,
        snippy_manifest,
        reference_path,
        snippy_dir,
        runs_dir,
        cpus,
        ram,
    )

    return {
        "samples": len(rows),
        "skipped": len(skipped),
        "snippy_manifest": snippy_manifest,
        "commands": commands_path,
        "readme": readme_path,
        "runs_dir": runs_dir,
        "reference": reference_path,
    }


def count_data_lines(path, comment_prefix=None, has_header=False):
    path = Path(path)
    if not path.exists():
        return 0

    count = 0
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if comment_prefix and line.startswith(comment_prefix):
                continue
            if line.strip():
                count += 1
    if has_header and count > 0:
        count -= 1
    return count


def fasta_count_and_length(path):
    path = Path(path)
    if not path.exists():
        return 0, 0

    record_count = 0
    first_length = 0
    current_length = 0
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if record_count == 1 and first_length == 0:
                    first_length = current_length
                record_count += 1
                current_length = 0
            else:
                current_length += len(line)
    if record_count == 1 and first_length == 0:
        first_length = current_length
    return record_count, first_length


def write_fasta_excluding_ids(input_path, output_path, excluded_ids):
    input_path = Path(input_path)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    excluded_ids = set(excluded_ids)

    total = 0
    written = 0
    excluded = 0
    write_current = False
    with input_path.open("r", encoding="utf-8", errors="replace") as source, output_path.open(
        "w",
        encoding="utf-8",
        newline="\n",
    ) as destination:
        for line in source:
            if line.startswith(">"):
                total += 1
                record_id = line[1:].strip().split()[0]
                write_current = record_id not in excluded_ids
                if write_current:
                    written += 1
                    destination.write(line)
                else:
                    excluded += 1
                continue
            if write_current:
                destination.write(line)

    return {
        "input": input_path,
        "output": output_path,
        "total_records": total,
        "written_records": written,
        "excluded_records": excluded,
    }


def read_snippy_core_stats(core_txt):
    core_txt = Path(core_txt)
    if not core_txt.exists():
        raise ValueError(f"Snippy core summary was not found: {core_txt}")

    stats_df = pd.read_csv(core_txt, sep="\t", keep_default_na=False, dtype=str)
    numeric_columns = ["LENGTH", "ALIGNED", "UNALIGNED", "VARIANT", "HET", "MASKED", "LOWCOV"]
    for column in numeric_columns:
        if column in stats_df.columns:
            stats_df[column] = pd.to_numeric(stats_df[column], errors="coerce").fillna(0).astype("int64")

    if {"ALIGNED", "LENGTH"}.issubset(stats_df.columns):
        stats_df["ALIGNED_%"] = stats_df.apply(
            lambda row: round(row["ALIGNED"] / row["LENGTH"] * 100, 2) if row["LENGTH"] else 0,
            axis=1,
        )
    if {"UNALIGNED", "LENGTH"}.issubset(stats_df.columns):
        stats_df["UNALIGNED_%"] = stats_df.apply(
            lambda row: round(row["UNALIGNED"] / row["LENGTH"] * 100, 2) if row["LENGTH"] else 0,
            axis=1,
        )
    return stats_df


def build_snippy_run_status(runs_dir):
    runs_dir = Path(runs_dir)
    if not runs_dir.exists():
        raise ValueError(f"Snippy runs directory was not found: {runs_dir}")

    rows = []
    for run_dir in sorted(path for path in runs_dir.iterdir() if path.is_dir()):
        snps_tab = run_dir / "snps.tab"
        snps_vcf = run_dir / "snps.vcf"
        snps_aligned = run_dir / "snps.aligned.fa"
        snps_consensus = run_dir / "snps.consensus.fa"
        snps_log = run_dir / "snps.log"
        required_paths = [snps_tab, snps_vcf, snps_aligned, snps_consensus, snps_log]
        missing = [path.name for path in required_paths if not path.exists() or path.stat().st_size == 0]
        rows.append(
            {
                "sample_id": run_dir.name,
                "run_dir": str(run_dir),
                "status": "complete" if not missing else "incomplete",
                "missing_outputs": ";".join(missing),
                "snps_tab_rows": count_data_lines(snps_tab, has_header=True),
                "snps_vcf_variants": count_data_lines(snps_vcf, comment_prefix="#"),
                "snps_aligned_exists": snps_aligned.exists() and snps_aligned.stat().st_size > 0,
                "snps_consensus_exists": snps_consensus.exists() and snps_consensus.stat().st_size > 0,
                "snps_log_exists": snps_log.exists() and snps_log.stat().st_size > 0,
            }
        )
    return pd.DataFrame(rows)


def build_snippy_summary_table(core_stats_df, run_status_df, core_aln, core_tab, core_vcf):
    sample_rows = core_stats_df[core_stats_df["ID"].astype(str) != "Reference"] if "ID" in core_stats_df.columns else core_stats_df
    complete_runs = int(run_status_df["status"].eq("complete").sum()) if "status" in run_status_df.columns else 0
    fasta_records, alignment_length = fasta_count_and_length(core_aln)
    rows = [
        ("generated_at", datetime.now().astimezone().isoformat(timespec="seconds")),
        ("snippy_runs", len(run_status_df)),
        ("complete_runs", complete_runs),
        ("incomplete_runs", len(run_status_df) - complete_runs),
        ("core_stats_samples", len(sample_rows)),
        ("core_alignment_records", fasta_records),
        ("core_alignment_length", alignment_length),
        ("core_tab_variant_rows", count_data_lines(core_tab, has_header=True)),
        ("core_vcf_variant_rows", count_data_lines(core_vcf, comment_prefix="#")),
        ("median_sample_variants", float(sample_rows["VARIANT"].median()) if "VARIANT" in sample_rows.columns and not sample_rows.empty else 0),
        ("max_sample_variants", int(sample_rows["VARIANT"].max()) if "VARIANT" in sample_rows.columns and not sample_rows.empty else 0),
        ("min_sample_variants", int(sample_rows["VARIANT"].min()) if "VARIANT" in sample_rows.columns and not sample_rows.empty else 0),
    ]
    return pd.DataFrame(rows, columns=["metric", "value"])


def write_snippy_summary_report(report_path, summary_df, run_status_df, core_stats_df, core_aln, core_tab, core_vcf):
    incomplete = run_status_df[run_status_df["status"] != "complete"] if "status" in run_status_df.columns else pd.DataFrame()
    summary_lookup = {row["metric"]: row["value"] for _, row in summary_df.iterrows()}
    lines = [
        "# TB Snippy summary",
        "",
        "## Summary",
        "",
        "| metric | value |",
        "| --- | --- |",
    ]
    for _, row in summary_df.iterrows():
        lines.append(f"| {row['metric']} | {row['value']} |")

    lines.extend(
        [
            "",
            "## Key files",
            "",
            f"- Core SNP alignment: `{core_aln}`",
            f"- SNP matrix: `{core_tab}`",
            f"- Core VCF: `{core_vcf}`",
            "",
            "## Run validation",
            "",
            f"- Complete Snippy runs: {summary_lookup.get('complete_runs', 0)}",
            f"- Incomplete Snippy runs: {summary_lookup.get('incomplete_runs', 0)}",
            "",
        ]
    )

    if incomplete.empty:
        lines.append("- No incomplete runs detected.")
    else:
        for _, row in incomplete.iterrows():
            lines.append(f"- `{row['sample_id']}`: missing `{row['missing_outputs']}`")

    if "VARIANT" in core_stats_df.columns and "ID" in core_stats_df.columns:
        top_variants = core_stats_df[core_stats_df["ID"].astype(str) != "Reference"].sort_values("VARIANT", ascending=False).head(10)
        lines.extend(["", "## Highest per-sample variant counts", "", "| sample_id | variants | aligned_% |", "| --- | ---: | ---: |"])
        for _, row in top_variants.iterrows():
            lines.append(f"| {row['ID']} | {row['VARIANT']} | {row.get('ALIGNED_%', '')} |")

    write_text_lf(report_path, "\n".join(lines) + "\n")


def write_snippy_summary(
    runs_dir=DEFAULT_SNIPPY_RUNS_DIR,
    core_txt=DEFAULT_SNIPPY_CORE_TXT,
    core_aln=DEFAULT_SNIPPY_CORE_ALN,
    core_tab=DEFAULT_SNIPPY_CORE_TAB,
    core_vcf=DEFAULT_SNIPPY_CORE_VCF,
    output_xlsx=DEFAULT_SNIPPY_SUMMARY_XLSX,
    report_path=DEFAULT_SNIPPY_SUMMARY_REPORT,
):
    core_aln = Path(core_aln)
    core_tab = Path(core_tab)
    core_vcf = Path(core_vcf)
    for path, label in [(core_aln, "core alignment"), (core_tab, "core SNP table"), (core_vcf, "core VCF")]:
        if not path.exists():
            raise ValueError(f"Snippy {label} was not found: {path}")

    core_stats_df = read_snippy_core_stats(core_txt)
    run_status_df = build_snippy_run_status(runs_dir)
    summary_df = build_snippy_summary_table(core_stats_df, run_status_df, core_aln, core_tab, core_vcf)

    write_dataframes_to_excel(
        output_xlsx,
        [
            (TB_SNIPPY_SUMMARY_SHEET, summary_df),
            (TB_SNIPPY_RUN_STATUS_SHEET, run_status_df),
            (TB_SNIPPY_CORE_STATS_SHEET, core_stats_df),
        ],
    )
    write_snippy_summary_report(report_path, summary_df, run_status_df, core_stats_df, core_aln, core_tab, core_vcf)
    return {
        "output_xlsx": Path(output_xlsx),
        "report": Path(report_path),
        "runs": len(run_status_df),
        "complete_runs": int(run_status_df["status"].eq("complete").sum()) if "status" in run_status_df.columns else 0,
        "core_tab_rows": count_data_lines(core_tab, has_header=True),
        "core_vcf_rows": count_data_lines(core_vcf, comment_prefix="#"),
    }


def write_fasttree_readme(readme_path, alignment_path, tree_path, log_path, filter_info=None):
    lines = [
        "# TB phylogenetics",
        "",
        "This directory contains the FastTree result built from the Snippy core SNP alignment.",
        "",
        "## Files",
        "",
        f"- Input alignment: `{alignment_path}`",
        f"- Newick tree: `{tree_path}`",
        f"- FastTree log: `{log_path}`",
    ]
    if filter_info:
        lines.extend(
            [
                f"- Original alignment records: {filter_info['total_records']}",
                f"- Tree alignment records: {filter_info['written_records']}",
                f"- Excluded records: {filter_info['excluded_records']}",
            ]
        )
    lines.extend(
        [
            "",
            "The Newick tree can be imported into iTOL for visualization and annotation.",
        ]
    )
    write_text_lf(readme_path, "\n".join(lines) + "\n")


def run_fasttree_wsl(
    alignment_path=DEFAULT_SNIPPY_CORE_ALN,
    filtered_alignment_path=DEFAULT_SNIPPY_CORE_NO_REFERENCE_ALN,
    tree_path=DEFAULT_FASTTREE_TREE,
    log_path=DEFAULT_FASTTREE_LOG,
    readme_path=DEFAULT_FASTTREE_README,
    distro="Ubuntu",
    conda_env="snippy",
    conda_profile="~/miniforge3/etc/profile.d/conda.sh",
    fasttree_exe="FastTree",
    gtr=True,
    exclude_reference=True,
):
    alignment_path = Path(alignment_path)
    filtered_alignment_path = Path(filtered_alignment_path)
    tree_path = Path(tree_path)
    log_path = Path(log_path)
    readme_path = Path(readme_path)
    if not alignment_path.exists():
        raise ValueError(f"Core SNP alignment was not found: {alignment_path}")

    fasttree_alignment = alignment_path
    filter_info = None
    if exclude_reference:
        filter_info = write_fasta_excluding_ids(alignment_path, filtered_alignment_path, {"Reference"})
        if filter_info["written_records"] == 0:
            raise ValueError(f"No records left after excluding Reference from: {alignment_path}")
        fasttree_alignment = filtered_alignment_path

    tree_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    flags = "-nt -gtr" if gtr else "-nt"
    bash_command = " ".join(
        [
            f"source {conda_profile}",
            "&&",
            f"conda activate {shell_quote(conda_env)}",
            "&&",
            "mkdir -p",
            shell_quote(path_to_wsl(tree_path.parent)),
            "&&",
            f"{shell_quote(fasttree_exe)} {flags}",
            shell_quote(path_to_wsl(fasttree_alignment)),
            ">",
            shell_quote(path_to_wsl(tree_path)),
            "2>",
            shell_quote(path_to_wsl(log_path)),
        ]
    )
    subprocess.run(["wsl", "-d", distro, "-e", "bash", "-lc", bash_command], check=True)
    write_fasttree_readme(readme_path, fasttree_alignment, tree_path, log_path, filter_info=filter_info)
    return {
        "tree": tree_path,
        "log": log_path,
        "readme": readme_path,
        "alignment": fasttree_alignment,
        "source_alignment": alignment_path,
        "filter_info": filter_info,
        "tree_size": tree_path.stat().st_size if tree_path.exists() else 0,
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

    snippy_parser = subparsers.add_parser(
        "snippy-plan",
        help="Create WSL/Snippy command files from the prepared TB input manifest without running Snippy.",
    )
    snippy_parser.add_argument(
        "--manifest",
        type=Path,
        default=DEFAULT_INPUT_MANIFEST_TSV,
        help="Input TSV manifest from prepare-inputs.",
    )
    snippy_parser.add_argument(
        "--reference",
        type=Path,
        default=DEFAULT_SNIPPY_REFERENCE,
        help="Reference GenBank file for Snippy.",
    )
    snippy_parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_SNIPPY_DIR,
        help="Directory where Snippy plan files will be written.",
    )
    snippy_parser.add_argument(
        "--snippy-manifest",
        type=Path,
        default=None,
        help="Optional explicit output path for snippy_multi.tsv.",
    )
    snippy_parser.add_argument(
        "--commands",
        type=Path,
        default=None,
        help="Optional explicit output path for snippy_commands.sh.",
    )
    snippy_parser.add_argument(
        "--readme",
        type=Path,
        default=None,
        help="Optional explicit output path for README_snippy.md.",
    )
    snippy_parser.add_argument(
        "--runs-dir",
        type=Path,
        default=None,
        help="Optional explicit Snippy output directory. Defaults to <output-dir>/runs.",
    )
    snippy_parser.add_argument(
        "--cpus",
        type=int,
        default=4,
        help="CPU threads per Snippy sample command.",
    )
    snippy_parser.add_argument(
        "--ram",
        type=int,
        default=8,
        help="Snippy RAM setting per sample command, in GB.",
    )
    snippy_parser.add_argument(
        "--force",
        action="store_true",
        help="Add --force to generated Snippy commands.",
    )
    snippy_parser.add_argument(
        "--cleanup",
        action="store_true",
        help="Add --cleanup to generated Snippy commands.",
    )
    snippy_parser.add_argument(
        "--include-warnings",
        action="store_true",
        help="Include manifest rows with validation warnings.",
    )

    snippy_summary_parser = subparsers.add_parser(
        "snippy-summary",
        help="Validate completed Snippy outputs and write summary workbook/report.",
    )
    snippy_summary_parser.add_argument(
        "--runs-dir",
        type=Path,
        default=DEFAULT_SNIPPY_RUNS_DIR,
        help="Directory containing per-sample Snippy run folders.",
    )
    snippy_summary_parser.add_argument(
        "--core-txt",
        type=Path,
        default=DEFAULT_SNIPPY_CORE_TXT,
        help="Snippy core.txt summary file.",
    )
    snippy_summary_parser.add_argument(
        "--core-aln",
        type=Path,
        default=DEFAULT_SNIPPY_CORE_ALN,
        help="Snippy core SNP alignment.",
    )
    snippy_summary_parser.add_argument(
        "--core-tab",
        type=Path,
        default=DEFAULT_SNIPPY_CORE_TAB,
        help="Snippy core SNP table.",
    )
    snippy_summary_parser.add_argument(
        "--core-vcf",
        type=Path,
        default=DEFAULT_SNIPPY_CORE_VCF,
        help="Snippy core VCF.",
    )
    snippy_summary_parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_SNIPPY_SUMMARY_XLSX,
        help="Output Excel summary workbook.",
    )
    snippy_summary_parser.add_argument(
        "--report",
        type=Path,
        default=DEFAULT_SNIPPY_SUMMARY_REPORT,
        help="Output Markdown summary report.",
    )

    tree_parser = subparsers.add_parser(
        "build-tree",
        help="Build a FastTree Newick tree from the Snippy core SNP alignment through WSL.",
    )
    tree_parser.add_argument(
        "--alignment",
        type=Path,
        default=DEFAULT_SNIPPY_CORE_ALN,
        help="Input Snippy core SNP alignment.",
    )
    tree_parser.add_argument(
        "--filtered-alignment",
        type=Path,
        default=DEFAULT_SNIPPY_CORE_NO_REFERENCE_ALN,
        help="Output alignment used for tree building when Reference is excluded.",
    )
    tree_parser.add_argument(
        "--tree",
        type=Path,
        default=DEFAULT_FASTTREE_TREE,
        help="Output Newick tree path.",
    )
    tree_parser.add_argument(
        "--log",
        type=Path,
        default=DEFAULT_FASTTREE_LOG,
        help="Output FastTree log path.",
    )
    tree_parser.add_argument(
        "--readme",
        type=Path,
        default=DEFAULT_FASTTREE_README,
        help="Output phylogenetics README path.",
    )
    tree_parser.add_argument(
        "--wsl-distro",
        default="Ubuntu",
        help="WSL distribution name.",
    )
    tree_parser.add_argument(
        "--conda-env",
        default="snippy",
        help="Conda environment containing FastTree.",
    )
    tree_parser.add_argument(
        "--conda-profile",
        default="~/miniforge3/etc/profile.d/conda.sh",
        help="Conda shell profile path inside WSL.",
    )
    tree_parser.add_argument(
        "--fasttree-exe",
        default="FastTree",
        help="FastTree executable name or path inside WSL.",
    )
    tree_parser.add_argument(
        "--no-gtr",
        action="store_true",
        help="Use FastTree -nt without -gtr.",
    )
    tree_parser.add_argument(
        "--include-reference",
        action="store_true",
        help="Keep the Snippy Reference sequence in the FastTree alignment.",
    )

    mutation_parser = subparsers.add_parser(
        "mutation-summary",
        help="Build assembly-based mutation summaries from per-sample Snippy snps.tab files.",
    )
    mutation_parser.add_argument(
        "--runs-dir",
        type=Path,
        default=DEFAULT_SNIPPY_RUNS_DIR,
        help="Directory containing per-sample Snippy run folders.",
    )
    mutation_parser.add_argument(
        "--manifest",
        type=Path,
        default=DEFAULT_INPUT_MANIFEST_TSV,
        help="TB input manifest with metadata.",
    )
    mutation_parser.add_argument(
        "--variants-tsv",
        type=Path,
        default=DEFAULT_VARIANTS_TSV,
        help="Output combined Snippy variants TSV.",
    )
    mutation_parser.add_argument(
        "--gene-summary",
        type=Path,
        default=DEFAULT_GENE_SUMMARY_TSV,
        help="Output gene mutation summary TSV.",
    )
    mutation_parser.add_argument(
        "--sample-summary",
        type=Path,
        default=DEFAULT_SAMPLE_SUMMARY_TSV,
        help="Output sample mutation summary TSV.",
    )
    mutation_parser.add_argument(
        "--target-matrix",
        type=Path,
        default=DEFAULT_TARGET_MATRIX_TSV,
        help="Output target gene count matrix TSV.",
    )
    mutation_parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_MUTATION_SUMMARY_XLSX,
        help="Output Excel workbook.",
    )
    mutation_parser.add_argument(
        "--report",
        type=Path,
        default=DEFAULT_MUTATION_REPORT,
        help="Output Markdown report.",
    )

    itol_parser = subparsers.add_parser(
        "itol-export",
        help="Create iTOL annotation datasets from mutation summaries and metadata.",
    )
    itol_parser.add_argument(
        "--sample-summary",
        type=Path,
        default=DEFAULT_SAMPLE_SUMMARY_TSV,
        help="Sample mutation summary TSV from mutation-summary.",
    )
    itol_parser.add_argument(
        "--target-matrix",
        type=Path,
        default=DEFAULT_TARGET_MATRIX_TSV,
        help="Target gene matrix TSV from mutation-summary.",
    )
    itol_parser.add_argument(
        "--tree",
        type=Path,
        default=DEFAULT_FASTTREE_TREE,
        help="Newick tree path for README reference.",
    )
    itol_parser.add_argument(
        "--country-strip",
        type=Path,
        default=DEFAULT_ITOL_COUNTRY_STRIP,
        help="Output iTOL country color-strip file.",
    )
    itol_parser.add_argument(
        "--variant-bars",
        type=Path,
        default=DEFAULT_ITOL_VARIANT_BARS,
        help="Output iTOL variant-count bar file.",
    )
    itol_parser.add_argument(
        "--target-heatmap",
        type=Path,
        default=DEFAULT_ITOL_TARGET_HEATMAP,
        help="Output iTOL target-gene heatmap file.",
    )
    itol_parser.add_argument(
        "--readme",
        type=Path,
        default=DEFAULT_ITOL_README,
        help="Output iTOL README path.",
    )

    report_parser = subparsers.add_parser(
        "report",
        help="Write final assembly-based tuberculosis analysis report.",
    )
    report_parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_TB_REPORT,
        help="Output Markdown report path.",
    )
    report_parser.add_argument(
        "--manifest",
        type=Path,
        default=DEFAULT_INPUT_MANIFEST_TSV,
        help="TB input manifest.",
    )
    report_parser.add_argument(
        "--snippy-summary-report",
        type=Path,
        default=DEFAULT_SNIPPY_SUMMARY_REPORT,
        help="Snippy summary Markdown report.",
    )
    report_parser.add_argument(
        "--mutation-report",
        type=Path,
        default=DEFAULT_MUTATION_REPORT,
        help="Mutation summary Markdown report.",
    )
    report_parser.add_argument(
        "--sample-summary",
        type=Path,
        default=DEFAULT_SAMPLE_SUMMARY_TSV,
        help="Sample mutation summary TSV.",
    )
    report_parser.add_argument(
        "--gene-summary",
        type=Path,
        default=DEFAULT_GENE_SUMMARY_TSV,
        help="Gene mutation summary TSV.",
    )
    report_parser.add_argument(
        "--target-matrix",
        type=Path,
        default=DEFAULT_TARGET_MATRIX_TSV,
        help="Target gene matrix TSV.",
    )
    report_parser.add_argument(
        "--tree",
        type=Path,
        default=DEFAULT_FASTTREE_TREE,
        help="FastTree Newick tree path.",
    )
    report_parser.add_argument(
        "--fasttree-log",
        type=Path,
        default=DEFAULT_FASTTREE_LOG,
        help="FastTree log path.",
    )
    report_parser.add_argument(
        "--core-aln",
        type=Path,
        default=DEFAULT_SNIPPY_CORE_ALN,
        help="Snippy core alignment path.",
    )
    report_parser.add_argument(
        "--core-tab",
        type=Path,
        default=DEFAULT_SNIPPY_CORE_TAB,
        help="Snippy core SNP matrix path.",
    )
    report_parser.add_argument(
        "--core-txt",
        type=Path,
        default=DEFAULT_SNIPPY_CORE_TXT,
        help="Snippy core summary path.",
    )


def parse_args():
    parser = argparse.ArgumentParser(description="Build tuberculosis metadata/counts workbook.")
    add_tuberculosis_arguments(parser)
    return parser.parse_args()


def run_from_args(args):
    tb_command = getattr(args, "tb_command", "metadata")

    if tb_command == "snippy-plan":
        result = create_snippy_plan(
            manifest_path=args.manifest,
            reference_path=args.reference,
            snippy_dir=args.output_dir,
            snippy_manifest=args.snippy_manifest,
            commands_path=args.commands,
            readme_path=args.readme,
            runs_dir=args.runs_dir,
            cpus=args.cpus,
            ram=args.ram,
            force=args.force,
            cleanup=args.cleanup,
            include_warnings=args.include_warnings,
        )
        print(f"Snippy multi manifest: {result['snippy_manifest']}")
        print(f"Snippy command script: {result['commands']}")
        print(f"Snippy README: {result['readme']}")
        print(f"Planned Snippy run directory: {result['runs_dir']}")
        print(f"Reference GenBank: {result['reference']}")
        print(f"Samples planned: {result['samples']}")
        print(f"Samples skipped: {result['skipped']}")
        return

    if tb_command == "snippy-summary":
        result = write_snippy_summary(
            runs_dir=args.runs_dir,
            core_txt=args.core_txt,
            core_aln=args.core_aln,
            core_tab=args.core_tab,
            core_vcf=args.core_vcf,
            output_xlsx=args.output,
            report_path=args.report,
        )
        print(f"Snippy summary workbook: {result['output_xlsx']}")
        print(f"Snippy summary report: {result['report']}")
        print(f"Snippy runs: {result['runs']}")
        print(f"Complete Snippy runs: {result['complete_runs']}")
        print(f"Core SNP table rows: {result['core_tab_rows']}")
        print(f"Core VCF rows: {result['core_vcf_rows']}")
        return

    if tb_command == "build-tree":
        result = run_fasttree_wsl(
            alignment_path=args.alignment,
            filtered_alignment_path=args.filtered_alignment,
            tree_path=args.tree,
            log_path=args.log,
            readme_path=args.readme,
            distro=args.wsl_distro,
            conda_env=args.conda_env,
            conda_profile=args.conda_profile,
            fasttree_exe=args.fasttree_exe,
            gtr=not args.no_gtr,
            exclude_reference=not args.include_reference,
        )
        print(f"FastTree input alignment: {result['alignment']}")
        print(f"FastTree Newick tree: {result['tree']}")
        print(f"FastTree log: {result['log']}")
        print(f"Phylogenetics README: {result['readme']}")
        print(f"Tree size: {result['tree_size']} bytes")
        return

    if tb_command == "mutation-summary":
        result = write_mutation_summary(
            runs_dir=args.runs_dir,
            manifest_path=args.manifest,
            variants_tsv=args.variants_tsv,
            gene_summary_tsv=args.gene_summary,
            sample_summary_tsv=args.sample_summary,
            target_matrix_tsv=args.target_matrix,
            output_xlsx=args.output,
            report_path=args.report,
        )
        print(f"Mutation summary workbook: {result['output_xlsx']}")
        print(f"Mutation summary report: {result['report']}")
        print(f"Combined variants TSV: {result['variants_tsv']}")
        print(f"Gene summary TSV: {result['gene_summary_tsv']}")
        print(f"Sample summary TSV: {result['sample_summary_tsv']}")
        print(f"Target matrix TSV: {result['target_matrix_tsv']}")
        print(f"Samples: {result['samples']}")
        print(f"Variant rows: {result['variant_rows']}")
        print(f"Genes/loci with variants: {result['genes']}")
        print(f"Target genes with variants: {result['target_genes_with_variants']}")
        return

    if tb_command == "itol-export":
        result = export_itol_annotations(
            sample_summary_tsv=args.sample_summary,
            target_matrix_tsv=args.target_matrix,
            tree_path=args.tree,
            country_strip=args.country_strip,
            variant_bars=args.variant_bars,
            target_heatmap=args.target_heatmap,
            readme_path=args.readme,
        )
        print(f"iTOL country strip: {result['country_strip']}")
        print(f"iTOL variant bars: {result['variant_bars']}")
        print(f"iTOL target gene heatmap: {result['target_heatmap']}")
        print(f"iTOL README: {result['readme']}")
        return

    if tb_command == "report":
        result = write_final_report(
            report_path=args.output,
            manifest_path=args.manifest,
            snippy_summary_report=args.snippy_summary_report,
            mutation_report=args.mutation_report,
            sample_summary_tsv=args.sample_summary,
            gene_summary_tsv=args.gene_summary,
            target_matrix_tsv=args.target_matrix,
            tree_path=args.tree,
            fasttree_log=args.fasttree_log,
            core_aln=args.core_aln,
            core_tab=args.core_tab,
            core_txt=args.core_txt,
        )
        print(f"TB analysis report: {result['report']}")
        print(f"Samples: {result['samples']}")
        print(f"Genes/loci with variants: {result['genes']}")
        print(f"Target genes tracked: {result['target_genes']}")
        return

    if tb_command == "prepare-inputs":
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
