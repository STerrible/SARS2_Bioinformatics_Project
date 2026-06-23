import argparse
import csv
import json
import sys
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen


PROJECT_DIR = Path(__file__).resolve().parents[1]
DEFAULT_ACCESSIONS = PROJECT_DIR / "data" / "raw" / "tuberculosis_data" / "accessions_103.txt"
DEFAULT_OUTPUT_DIR = PROJECT_DIR / "data" / "raw" / "tuberculosis_data" / "ncbi_assemblies"
USER_AGENT = "SARS2_Bioinformatics_Project/1.0"

FILE_TYPES = {
    "fna": "_genomic.fna.gz",
    "gbff": "_genomic.gbff.gz",
    "report": "_assembly_report.txt",
}


def read_accessions(path):
    accessions = [line.strip() for line in Path(path).read_text(encoding="utf-8").splitlines()]
    accessions = [accession for accession in accessions if accession and not accession.startswith("#")]
    duplicates = sorted({accession for accession in accessions if accessions.count(accession) > 1})
    if duplicates:
        raise ValueError(f"Duplicate accessions in {path}: {', '.join(duplicates)}")
    return accessions


def request_json(url, data=None):
    encoded_data = None if data is None else urlencode(data).encode("utf-8")
    request = Request(
        url,
        data=encoded_data,
        headers={
            "User-Agent": USER_AGENT,
            "Content-Type": "application/x-www-form-urlencoded",
        },
    )
    with urlopen(request, timeout=90) as response:
        return json.load(response)


def fetch_assembly_ids(accessions):
    term = " OR ".join(f"{accession}[Assembly Accession]" for accession in accessions)
    data = {
        "db": "assembly",
        "term": term,
        "retmode": "json",
        "retmax": str(max(200, len(accessions))),
    }
    result = request_json("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi", data)
    return result.get("esearchresult", {}).get("idlist", [])


def fetch_assembly_summaries(assembly_ids):
    if not assembly_ids:
        return {}
    data = {
        "db": "assembly",
        "id": ",".join(assembly_ids),
        "retmode": "json",
    }
    result = request_json("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi", data)
    summaries = result.get("result", {})
    return {uid: summaries[uid] for uid in summaries.get("uids", [])}


def normalize_ftp_url(ftp_path):
    if not ftp_path:
        return ""
    return ftp_path.replace("ftp://ftp.ncbi.nlm.nih.gov", "https://ftp.ncbi.nlm.nih.gov")


def build_metadata(accessions):
    assembly_ids = fetch_assembly_ids(accessions)
    summaries = fetch_assembly_summaries(assembly_ids)
    by_accession = {entry.get("assemblyaccession"): entry for entry in summaries.values()}
    return [by_accession.get(accession) for accession in accessions]


def write_metadata(metadata_rows, output_path):
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "accession",
        "uid",
        "latest_accession",
        "organism",
        "species_name",
        "taxid",
        "assembly_name",
        "assembly_status",
        "refseq_category",
        "biosample",
        "bioproject",
        "ftp_path_refseq",
        "ftp_path_genbank",
    ]
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for entry in metadata_rows:
            writer.writerow(
                {
                    "accession": entry.get("assemblyaccession", ""),
                    "uid": entry.get("uid", ""),
                    "latest_accession": entry.get("latestaccession", ""),
                    "organism": entry.get("organism", ""),
                    "species_name": entry.get("speciesname", ""),
                    "taxid": entry.get("taxid", ""),
                    "assembly_name": entry.get("assemblyname", ""),
                    "assembly_status": entry.get("assemblystatus", ""),
                    "refseq_category": entry.get("refseq_category", ""),
                    "biosample": entry.get("biosampleaccn", ""),
                    "bioproject": ";".join(project.get("bioprojectaccn", "") for project in entry.get("rs_bioprojects", [])),
                    "ftp_path_refseq": entry.get("ftppath_refseq", ""),
                    "ftp_path_genbank": entry.get("ftppath_genbank", ""),
                }
            )


def download_file(url, destination, overwrite=False):
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and destination.stat().st_size > 0 and not overwrite:
        return "skipped"

    partial = destination.with_name(destination.name + ".part")
    try:
        with urlopen(Request(url, headers={"User-Agent": USER_AGENT}), timeout=180) as response:
            with partial.open("wb") as handle:
                while True:
                    chunk = response.read(1024 * 1024)
                    if not chunk:
                        break
                    handle.write(chunk)
        partial.replace(destination)
        return "downloaded"
    except HTTPError as exc:
        if partial.exists():
            partial.unlink()
        return f"http_{exc.code}"
    except URLError as exc:
        if partial.exists():
            partial.unlink()
        return f"url_error:{exc.reason}"


def build_downloads(entry, file_types):
    ftp_path = entry.get("ftppath_refseq") or entry.get("ftppath_genbank")
    base_url = normalize_ftp_url(ftp_path)
    if not base_url:
        return []
    base_name = base_url.rstrip("/").split("/")[-1]
    downloads = []
    for file_type in file_types:
        suffix = FILE_TYPES[file_type]
        filename = f"{base_name}{suffix}"
        downloads.append((file_type, f"{base_url}/{filename}", filename))
    return downloads


def write_manifest(rows, output_path):
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fields = ["accession", "file_type", "status", "path", "url"]
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def parse_args():
    parser = argparse.ArgumentParser(description="Download NCBI Assembly raw files for GCF/GCA accessions.")
    parser.add_argument("--accessions", type=Path, default=DEFAULT_ACCESSIONS)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--file-type",
        action="append",
        choices=sorted(FILE_TYPES),
        dest="file_types",
        help="File type to download. Repeat for multiple types. Default: fna, gbff, report.",
    )
    parser.add_argument("--metadata-only", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main():
    args = parse_args()
    file_types = args.file_types or ["fna", "gbff", "report"]
    accessions = read_accessions(args.accessions)
    metadata_rows = build_metadata(accessions)

    missing = [accession for accession, entry in zip(accessions, metadata_rows) if entry is None]
    if missing:
        raise SystemExit(f"Missing accessions in NCBI Assembly: {', '.join(missing)}")

    metadata_path = args.output_dir / "assembly_metadata.tsv"
    write_metadata(metadata_rows, metadata_path)
    print(f"Metadata written: {metadata_path}")

    if args.metadata_only:
        print(f"Accessions found: {len(metadata_rows)}")
        return

    manifest_rows = []
    for index, entry in enumerate(metadata_rows, start=1):
        accession = entry["assemblyaccession"]
        accession_dir = args.output_dir / accession
        for file_type, url, filename in build_downloads(entry, file_types):
            destination = accession_dir / filename
            status = download_file(url, destination, args.overwrite)
            manifest_rows.append(
                {
                    "accession": accession,
                    "file_type": file_type,
                    "status": status,
                    "path": str(destination),
                    "url": url,
                }
            )
            print(f"[{index}/{len(metadata_rows)}] {accession} {file_type}: {status}", flush=True)

    manifest_path = args.output_dir / "download_manifest.tsv"
    write_manifest(manifest_rows, manifest_path)
    failed = [row for row in manifest_rows if row["status"] not in {"downloaded", "skipped"}]
    print(f"Manifest written: {manifest_path}")
    print(f"Downloads requested: {len(manifest_rows)}")
    print(f"Downloads failed: {len(failed)}")
    if failed:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
