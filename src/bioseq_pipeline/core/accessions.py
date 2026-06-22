import argparse
import re
import sys
from pathlib import Path

from Bio import SeqIO


PROJECT_DIR = Path(__file__).resolve().parents[3]
DEFAULT_FASTA = PROJECT_DIR / "data" / "raw" / "sequence.fasta"

ACCESSION_PATTERN = re.compile(
    r"(?<![A-Za-z0-9_])"
    r"("
    r"(?:[A-Z]{1,4}_\d+(?:\.\d+)?)"
    r"|"
    r"(?:[A-Z]{1,6}\d{5,12}(?:\.\d+)?)"
    r")"
    r"(?![A-Za-z0-9_])"
)


def strip_accession_version(accession):
    prefix, separator, version = accession.rpartition(".")
    if separator and version.isdigit():
        return prefix
    return accession


def find_accession(text):
    match = ACCESSION_PATTERN.search(str(text))
    return match.group(1) if match else ""


def extract_fasta_accessions(input_fasta, keep_version=True):
    input_fasta = Path(input_fasta)
    if not input_fasta.exists():
        raise ValueError(f"Input FASTA file was not found: {input_fasta}")

    accessions = []
    seen = set()
    for record in SeqIO.parse(input_fasta, "fasta"):
        accession = find_accession(record.id) or find_accession(record.description)
        if not accession:
            continue
        if not keep_version:
            accession = strip_accession_version(accession)
        if accession in seen:
            continue
        seen.add(accession)
        accessions.append(accession)

    return accessions


def write_accessions(accessions, output_path=None):
    text = "\n".join(accessions)
    if text:
        text += "\n"

    if output_path is None:
        sys.stdout.write(text)
        return

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(text, encoding="utf-8")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract NCBI-style accession numbers from FASTA record headers.",
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_FASTA,
        help="Input FASTA file.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Optional output TXT file. If omitted, accessions are printed to stdout.",
    )
    parser.add_argument(
        "--strip-version",
        action="store_true",
        help="Write accessions without version suffixes, for example NC_045512 instead of NC_045512.2.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    accessions = extract_fasta_accessions(args.input, keep_version=not args.strip_version)
    write_accessions(accessions, args.output)
    if args.output is not None:
        print(f"Accessions written: {args.output}")
        print(f"Accessions found: {len(accessions)}")


if __name__ == "__main__":
    main()
