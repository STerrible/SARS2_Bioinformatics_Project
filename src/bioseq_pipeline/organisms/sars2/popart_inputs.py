import argparse
import re
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from bioseq_pipeline.organisms.sars2.phylogenetics import (
    DEFAULT_ALIGNED_FASTA,
    DEFAULT_EXCEL,
    DEFAULT_OUTPUT_DIR,
    prepare_inputs,
)


POPART_TRAITS = {
    "country": "country",
    "month": "collection_month",
    "clade": "clade",
    "lineage": "Nextclade_pango",
    "d614g": "has_S_D614G",
}


def duplicate_values(values):
    seen = set()
    duplicates = set()
    for value in values:
        if value in seen:
            duplicates.add(value)
        seen.add(value)
    return sorted(duplicates)


def nexus_token(value, prefix="item"):
    text = str(value).strip()
    if not text:
        text = "unknown"
    text = re.sub(r"[^A-Za-z0-9_.-]+", "_", text).strip("_")
    if not text:
        text = "unknown"
    if not re.match(r"^[A-Za-z_]", text):
        text = f"{prefix}_{text}"
    return text


def normalize_trait_value(value):
    text = str(value).strip()
    return text if text else "unknown"


def read_alignment(alignment_path):
    records = list(SeqIO.parse(alignment_path, "fasta"))
    if not records:
        raise ValueError(f"No FASTA records were found in: {alignment_path}")

    sequence_lengths = {len(record.seq) for record in records}
    if len(sequence_lengths) != 1:
        raise ValueError(f"PopART expects aligned sequences of equal length: {alignment_path}")

    return records, sequence_lengths.pop()


def validate_nexus_file(nexus_path, expected_taxa, expected_sequence_length, expected_trait_count):
    nexus_path = Path(nexus_path)
    messages = []
    if not nexus_path.exists():
        return [f"ERROR: NEXUS file was not created: {nexus_path}"]

    content = nexus_path.read_text(encoding="utf-8")
    upper_content = content.upper()
    if not content.startswith("#NEXUS"):
        messages.append(f"ERROR: {nexus_path.name} does not start with #NEXUS.")

    for block_name in ["TAXA", "CHARACTERS", "TRAITS"]:
        if f"BEGIN {block_name};" not in upper_content:
            messages.append(f"ERROR: {nexus_path.name} is missing BEGIN {block_name}; block.")

    ntax_match = re.search(r"DIMENSIONS\s+NTAX\s*=\s*(\d+)", content, flags=re.IGNORECASE)
    nchar_match = re.search(r"DIMENSIONS\s+NCHAR\s*=\s*(\d+)", content, flags=re.IGNORECASE)
    ntraits_match = re.search(r"DIMENSIONS\s+NTRAITS\s*=\s*(\d+)", content, flags=re.IGNORECASE)

    if not ntax_match or int(ntax_match.group(1)) != expected_taxa:
        found = ntax_match.group(1) if ntax_match else "missing"
        messages.append(f"ERROR: {nexus_path.name} NTAX={found}, expected {expected_taxa}.")
    if not nchar_match or int(nchar_match.group(1)) != expected_sequence_length:
        found = nchar_match.group(1) if nchar_match else "missing"
        messages.append(f"ERROR: {nexus_path.name} NCHAR={found}, expected {expected_sequence_length}.")
    if not ntraits_match or int(ntraits_match.group(1)) != expected_trait_count:
        found = ntraits_match.group(1) if ntraits_match else "missing"
        messages.append(f"ERROR: {nexus_path.name} NTRAITS={found}, expected {expected_trait_count}.")

    return messages


def validate_popart_inputs(alignment_path, metadata_path, nexus_files, trait_columns):
    records, sequence_length = read_alignment(alignment_path)
    metadata = pd.read_csv(metadata_path, sep="\t", keep_default_na=False)
    errors = []
    warnings = []
    details = []

    record_ids = [record.id for record in records]
    sample_tokens = [nexus_token(record.id, prefix="sample") for record in records]
    duplicate_record_ids = duplicate_values(record_ids)
    duplicate_sample_tokens = duplicate_values(sample_tokens)
    if duplicate_record_ids:
        errors.append(f"Duplicate FASTA record IDs: {', '.join(duplicate_record_ids[:20])}")
    if duplicate_sample_tokens:
        errors.append(f"Duplicate NEXUS-safe sample names: {', '.join(duplicate_sample_tokens[:20])}")

    if "short_name" not in metadata.columns:
        errors.append(f"Metadata file must contain short_name column: {metadata_path}")
        return {
            "status": "FAIL",
            "errors": errors,
            "warnings": warnings,
            "details": details,
        }

    metadata_names = metadata["short_name"].astype(str).tolist()
    duplicate_metadata_names = duplicate_values(metadata_names)
    if duplicate_metadata_names:
        errors.append(f"Duplicate metadata short_name values: {', '.join(duplicate_metadata_names[:20])}")

    record_id_set = set(record_ids)
    metadata_name_set = set(metadata_names)
    missing_metadata = sorted(record_id_set - metadata_name_set)
    extra_metadata = sorted(metadata_name_set - record_id_set)
    if missing_metadata:
        errors.append(f"Samples missing from metadata: {', '.join(missing_metadata[:20])}")
    if extra_metadata:
        warnings.append(f"Metadata rows not present in alignment: {', '.join(extra_metadata[:20])}")

    details.append(("alignment_records", len(records)))
    details.append(("alignment_sequence_length", sequence_length))
    details.append(("metadata_rows", len(metadata)))
    details.append(("metadata_columns", len(metadata.columns)))

    for output_name, trait_column in trait_columns.items():
        if trait_column not in metadata.columns:
            errors.append(f"Trait column is missing for {output_name}: {trait_column}")
            continue

        trait_by_sample = {
            str(row["short_name"]): normalize_trait_value(row[trait_column])
            for _, row in metadata.iterrows()
        }
        trait_values = [trait_by_sample.get(record_id, "unknown") for record_id in record_ids]
        unknown_count = sum(value == "unknown" for value in trait_values)
        unique_trait_values = sorted(set(trait_values))
        if unknown_count:
            warnings.append(f"{output_name}: {unknown_count} samples have unknown trait values.")

        nexus_path = nexus_files.get(output_name)
        mapping_path = nexus_path.with_name(f"{nexus_path.stem}_trait_labels.tsv") if nexus_path else None
        if mapping_path is None or not mapping_path.exists():
            errors.append(f"{output_name}: trait label mapping file was not created.")
        else:
            mapping = pd.read_csv(mapping_path, sep="\t", keep_default_na=False)
            if len(mapping) != len(unique_trait_values):
                errors.append(
                    f"{output_name}: mapping rows={len(mapping)}, expected {len(unique_trait_values)} trait values."
                )

        if nexus_path is not None:
            errors.extend(
                validate_nexus_file(
                    nexus_path,
                    len(records),
                    sequence_length,
                    len(unique_trait_values),
                )
            )

        details.append((f"{output_name}_trait_column", trait_column))
        details.append((f"{output_name}_trait_values", len(unique_trait_values)))
        details.append((f"{output_name}_unknown_values", unknown_count))

    status = "FAIL" if errors else "WARN" if warnings else "PASS"
    return {
        "status": status,
        "errors": errors,
        "warnings": warnings,
        "details": details,
    }


def write_validation_report(output_path, validation):
    lines = [
        "# PopART validation report",
        "",
        f"Status: **{validation['status']}**",
        "",
        "## Details",
        "",
        "| metric | value |",
        "| --- | --- |",
    ]
    for metric, value in validation["details"]:
        lines.append(f"| {metric} | {value} |")

    lines.extend(["", "## Errors", ""])
    if validation["errors"]:
        lines.extend(f"- {message}" for message in validation["errors"])
    else:
        lines.append("- None")

    lines.extend(["", "## Warnings", ""])
    if validation["warnings"]:
        lines.extend(f"- {message}" for message in validation["warnings"])
    else:
        lines.append("- None")

    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_popart_nexus(alignment_path, metadata_path, trait_column, output_path):
    records, sequence_length = read_alignment(alignment_path)
    metadata = pd.read_csv(metadata_path, sep="\t", keep_default_na=False)
    if "short_name" not in metadata.columns:
        raise ValueError(f"Metadata file must contain a short_name column: {metadata_path}")
    if trait_column not in metadata.columns:
        raise ValueError(f"Metadata file does not contain trait column: {trait_column}")

    trait_by_sample = {
        str(row["short_name"]): normalize_trait_value(row[trait_column])
        for _, row in metadata.iterrows()
    }
    raw_trait_labels = sorted({trait_by_sample.get(record.id, "unknown") for record in records})
    trait_labels = []
    trait_label_mapping = []
    used_trait_tokens = set()
    for raw_trait in raw_trait_labels:
        base_token = nexus_token(raw_trait, prefix="trait")
        token = base_token
        index = 2
        while token in used_trait_tokens:
            token = f"{base_token}_{index}"
            index += 1
        used_trait_tokens.add(token)
        trait_labels.append(token)
        trait_label_mapping.append((token, raw_trait))

    trait_index = {
        raw_trait: index
        for index, raw_trait in enumerate(raw_trait_labels)
    }

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("#NEXUS\n\n")

        handle.write("BEGIN TAXA;\n")
        handle.write(f"  DIMENSIONS NTAX={len(records)};\n")
        handle.write("  TAXLABELS\n")
        for record in records:
            handle.write(f"    {nexus_token(record.id, prefix='sample')}\n")
        handle.write("  ;\n")
        handle.write("END;\n\n")

        handle.write("BEGIN CHARACTERS;\n")
        handle.write(f"  DIMENSIONS NCHAR={sequence_length};\n")
        handle.write("  FORMAT DATATYPE=DNA MISSING=? GAP=-;\n")
        handle.write("  MATRIX\n")
        for record in records:
            handle.write(f"    {nexus_token(record.id, prefix='sample')} {str(record.seq).upper()}\n")
        handle.write("  ;\n")
        handle.write("END;\n\n")

        handle.write("BEGIN TRAITS;\n")
        handle.write(f"  DIMENSIONS NTRAITS={len(trait_labels)};\n")
        handle.write("  FORMAT LABELS=YES MISSING=? SEPARATOR=Comma;\n")
        handle.write("  TRAITLABELS\n")
        for trait in trait_labels:
            handle.write(f"    {trait}\n")
        handle.write("  ;\n")
        handle.write("  MATRIX\n")
        for record in records:
            values = ["0"] * len(trait_labels)
            values[trait_index[trait_by_sample.get(record.id, "unknown")]] = "1"
            handle.write(f"    {nexus_token(record.id, prefix='sample')} {','.join(values)}\n")
        handle.write("  ;\n")
        handle.write("END;\n")

    mapping_path = output_path.with_name(f"{output_path.stem}_trait_labels.tsv")
    with mapping_path.open("w", encoding="utf-8", newline="") as handle:
        handle.write("popart_label\toriginal_label\n")
        for token, raw_trait in trait_label_mapping:
            handle.write(f"{token}\t{raw_trait}\n")

    return output_path


def write_popart_readme(output_path, nexus_files):
    lines = [
        "# Входные файлы для PopART",
        "",
        "Эти NEXUS-файлы подготовлены для анализа гаплотипической сети в PopART.",
        "",
        "Рекомендуемый порядок работы в PopART:",
        "",
        "1. Откройте один из файлов `popart_*.nex` в PopART.",
        "2. Выберите метод построения сети, обычно Median Joining или TCS.",
        "3. Используйте легенду признаков для окраски по стране, месяцу, clade, lineage или D614G.",
        "4. Экспортируйте итоговую фигуру из PopART в SVG/PNG/PDF.",
        "",
        "Воспроизводимость:",
        "",
        "- Лучше пересоздавать эти файлы через Docker workflow, описанный в основном README проекта.",
        "- Docker-образ фиксирует версию Nextclade CLI и tag SARS-CoV-2 dataset.",
        "",
        "Сгенерированные файлы:",
        "",
    ]
    for name, path in nexus_files.items():
        lines.append(f"- `{path.name}`: окраска по `{name}`.")
    lines.extend(
        [
            "",
            "Для каждого NEXUS-файла создается парный файл `*_trait_labels.tsv`.",
            "Он сопоставляет PopART-safe метки, например `South_Korea`, с исходными значениями метаданных.",
            "",
            "Примечание: для 35 геномов сеть будет намного меньше, чем крупные схемы из дипломных материалов.",
            "Чтобы получить сопоставимую плотную сеть, нужен набор с большим числом геномов по разным месяцам и странам.",
            "",
        ]
    )
    output_path.write_text("\n".join(lines), encoding="utf-8")


def create_popart_inputs(aligned_fasta, excel_path, output_dir):
    prepared = prepare_inputs(aligned_fasta, excel_path, output_dir)
    output_dir = Path(output_dir)
    alignment_path = prepared["alignment"]
    metadata_path = prepared["metadata"]

    nexus_files = {}
    for output_name, trait_column in POPART_TRAITS.items():
        nexus_files[output_name] = write_popart_nexus(
            alignment_path,
            metadata_path,
            trait_column,
            output_dir / f"popart_{output_name}.nex",
        )

    readme_path = output_dir / "README_popart.md"
    write_popart_readme(readme_path, nexus_files)
    validation = validate_popart_inputs(alignment_path, metadata_path, nexus_files, POPART_TRAITS)
    validation_report_path = output_dir / "popart_validation_report.md"
    write_validation_report(validation_report_path, validation)
    if validation["errors"]:
        raise ValueError(f"PopART input validation failed. See: {validation_report_path}")

    return {
        **prepared,
        "nexus_files": nexus_files,
        "popart_readme": readme_path,
        "validation_report": validation_report_path,
        "validation_status": validation["status"],
    }


def parse_args():
    parser = argparse.ArgumentParser(
        description="Create PopART NEXUS files from aligned FASTA and pipeline metadata.",
    )
    parser.add_argument(
        "--aligned-fasta",
        type=Path,
        default=DEFAULT_ALIGNED_FASTA,
        help="Input aligned FASTA file.",
    )
    parser.add_argument(
        "--excel",
        type=Path,
        default=DEFAULT_EXCEL,
        help="Input Excel file generated by the main pipeline.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Output directory for PopART input files.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    outputs = create_popart_inputs(args.aligned_fasta, args.excel, args.output_dir)
    print(f"Done: {args.output_dir}")
    print(f"Records processed: {outputs['record_count']}")
    for name, path in outputs["nexus_files"].items():
        print(f"PopART {name}: {path}")
    print(f"Validation: {outputs['validation_status']}")
    print(f"Validation report: {outputs['validation_report']}")


if __name__ == "__main__":
    main()
