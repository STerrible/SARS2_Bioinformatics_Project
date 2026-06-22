from pathlib import Path
import json
import shutil
import subprocess

from Bio import SeqIO
import pandas as pd


def resolve_nextclade_executable(nextclade_exe):
    value = str(nextclade_exe)
    path = Path(value)
    if path.exists():
        return str(path)

    executable = shutil.which(value)
    if executable:
        return executable

    raise ValueError(f"Nextclade executable was not found: {nextclade_exe}")


def read_nextclade_dataset_info(dataset_dir):
    dataset_dir = Path(dataset_dir)
    pathogen_json = dataset_dir / "pathogen.json"
    info = {
        "nextclade_dataset_path": str(dataset_dir),
        "nextclade_dataset_metadata_file": str(pathogen_json),
    }

    if not pathogen_json.exists():
        info["nextclade_dataset_metadata_status"] = "missing"
        return info

    try:
        data = json.loads(pathogen_json.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        info["nextclade_dataset_metadata_status"] = f"invalid_json: {exc}"
        return info

    version = data.get("version") or {}
    compatibility = version.get("compatibility") or {}
    attributes = data.get("attributes") or {}
    shortcuts = data.get("shortcuts") or []

    info.update({
        "nextclade_dataset_metadata_status": "ok",
        "nextclade_dataset_schema_version": data.get("schemaVersion", ""),
        "nextclade_dataset_version_tag": version.get("tag", ""),
        "nextclade_dataset_updated_at": version.get("updatedAt", ""),
        "nextclade_dataset_cli_compatibility": compatibility.get("cli", ""),
        "nextclade_dataset_web_compatibility": compatibility.get("web", ""),
        "nextclade_dataset_reference_name": attributes.get("reference name", ""),
        "nextclade_dataset_reference_accession": attributes.get("reference accession", ""),
        "nextclade_dataset_shortcuts": "; ".join(str(shortcut) for shortcut in shortcuts),
    })
    return info


def run_nextclade(input_fasta, output_tsv, nextclade_exe, dataset_dir, output_aligned_fasta=None):
    input_fasta = Path(input_fasta)
    output_tsv = Path(output_tsv)
    nextclade_command = resolve_nextclade_executable(nextclade_exe)
    dataset_dir = Path(dataset_dir)
    output_aligned_fasta = Path(output_aligned_fasta) if output_aligned_fasta is not None else None

    if not input_fasta.exists():
        raise ValueError(f"Input FASTA file was not found: {input_fasta}")
    if not dataset_dir.exists():
        raise ValueError(f"Nextclade dataset directory was not found: {dataset_dir}")

    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    if output_aligned_fasta is not None:
        output_aligned_fasta.parent.mkdir(parents=True, exist_ok=True)

    command = [
        nextclade_command,
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
    nextclade_command = resolve_nextclade_executable(nextclade_exe)

    try:
        result = subprocess.run(
            [nextclade_command, "--version"],
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
    df = pd.read_csv(output_tsv, sep="\t", keep_default_na=False)
    if "seqName" in df.columns:
        df = df.sort_values("seqName", kind="stable").reset_index(drop=True)
    return df


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
