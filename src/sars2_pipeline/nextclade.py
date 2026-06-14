from pathlib import Path
import subprocess

from Bio import SeqIO
import pandas as pd


def run_nextclade(input_fasta, output_tsv, nextclade_exe, dataset_dir, output_aligned_fasta=None):
    input_fasta = Path(input_fasta)
    output_tsv = Path(output_tsv)
    nextclade_exe = Path(nextclade_exe)
    dataset_dir = Path(dataset_dir)
    output_aligned_fasta = Path(output_aligned_fasta) if output_aligned_fasta is not None else None

    if not input_fasta.exists():
        raise ValueError(f"Input FASTA file was not found: {input_fasta}")
    if not nextclade_exe.exists():
        raise ValueError(f"Nextclade executable was not found: {nextclade_exe}")
    if not dataset_dir.exists():
        raise ValueError(f"Nextclade dataset directory was not found: {dataset_dir}")

    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    if output_aligned_fasta is not None:
        output_aligned_fasta.parent.mkdir(parents=True, exist_ok=True)

    command = [
        str(nextclade_exe),
        "run",
        "--input-dataset",
        str(dataset_dir),
        "--output-tsv",
        str(output_tsv),
    ]
    if output_aligned_fasta is not None:
        command.extend(["--output-fasta", str(output_aligned_fasta)])

    command.extend([
        str(input_fasta),
    ])

    try:
        subprocess.run(command, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as exc:
        details = (exc.stderr or exc.stdout or "").strip()
        message = f"Nextclade failed with exit code {exc.returncode}."
        if details:
            message = f"{message}\n{details}"
        raise ValueError(message) from exc

    return output_tsv


def get_nextclade_version(nextclade_exe):
    nextclade_exe = Path(nextclade_exe)
    if not nextclade_exe.exists():
        raise ValueError(f"Nextclade executable was not found: {nextclade_exe}")

    try:
        result = subprocess.run(
            [str(nextclade_exe), "--version"],
            check=True,
            capture_output=True,
            text=True,
        )
    except subprocess.CalledProcessError as exc:
        details = (exc.stderr or exc.stdout or "").strip()
        message = f"Failed to read Nextclade version with exit code {exc.returncode}."
        if details:
            message = f"{message}\n{details}"
        raise ValueError(message) from exc

    return result.stdout.strip()


def read_nextclade_tsv(output_tsv):
    output_tsv = Path(output_tsv)
    if not output_tsv.exists():
        raise ValueError(f"Nextclade TSV file was not found: {output_tsv}")
    return pd.read_csv(output_tsv, sep="\t", keep_default_na=False)


def read_aligned_fasta(aligned_fasta):
    aligned_fasta = Path(aligned_fasta)
    if not aligned_fasta.exists():
        raise ValueError(f"Aligned FASTA file was not found: {aligned_fasta}")

    rows = []
    for record in SeqIO.parse(aligned_fasta, "fasta"):
        rows.append({
            "seq_id": record.id,
            "description": record.description,
            "aligned_sequence": str(record.seq).upper(),
        })
    return pd.DataFrame(rows, columns=["seq_id", "description", "aligned_sequence"])
