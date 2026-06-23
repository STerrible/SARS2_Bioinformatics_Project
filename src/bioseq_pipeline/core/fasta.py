import argparse
import re
from pathlib import Path

from Bio import SeqIO


PROJECT_DIR = Path(__file__).resolve().parents[3]
DEFAULT_INPUT = PROJECT_DIR / "results" / "covid_data" / "aligned.fasta"
DEFAULT_OUTPUT = PROJECT_DIR / "results" / "covid_data" / "phylogenetics" / "alignment_short.fasta"
DEFAULT_MAPPING = PROJECT_DIR / "results" / "covid_data" / "phylogenetics" / "alignment_short_mapping.tsv"


def short_name(record_id, used_names):
    name = record_id.split()[0]
    name = re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("_")
    if not name:
        name = "sequence"

    unique_name = name
    index = 2
    while unique_name in used_names:
        unique_name = f"{name}_{index}"
        index += 1

    used_names.add(unique_name)
    return unique_name


def create_short_alignment(input_fasta, output_fasta, mapping_tsv):
    input_fasta = Path(input_fasta)
    output_fasta = Path(output_fasta)
    mapping_tsv = Path(mapping_tsv)

    if not input_fasta.exists():
        raise ValueError(f"Input aligned FASTA was not found: {input_fasta}")

    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    mapping_tsv.parent.mkdir(parents=True, exist_ok=True)

    used_names = set()
    records = []
    mapping_rows = []

    for record in SeqIO.parse(input_fasta, "fasta"):
        original_description = record.description
        new_id = short_name(record.id, used_names)
        record.id = new_id
        record.name = new_id
        record.description = ""
        records.append(record)
        mapping_rows.append((new_id, original_description))

    if not records:
        raise ValueError(f"No FASTA records were found in: {input_fasta}")

    SeqIO.write(records, output_fasta, "fasta")

    with mapping_tsv.open("w", encoding="utf-8", newline="") as handle:
        handle.write("short_name\toriginal_description\n")
        for new_id, original_description in mapping_rows:
            handle.write(f"{new_id}\t{original_description}\n")

    return len(records)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Create an aligned FASTA file with short sequence names for tree-building tools.",
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help="Input aligned FASTA file.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help="Output FASTA file with short sequence names.",
    )
    parser.add_argument(
        "--mapping",
        type=Path,
        default=DEFAULT_MAPPING,
        help="TSV file mapping short names to original FASTA descriptions.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    record_count = create_short_alignment(args.input, args.output, args.mapping)
    print(f"Done: {args.output}")
    print(f"Mapping: {args.mapping}")
    print(f"Records processed: {record_count}")


if __name__ == "__main__":
    main()
