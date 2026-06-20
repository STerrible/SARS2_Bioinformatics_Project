import argparse
from datetime import datetime
from pathlib import Path

import pandas as pd

from create_alignment_short import create_short_alignment


PROJECT_DIR = Path(__file__).resolve().parents[1]
DEFAULT_ALIGNED_FASTA = PROJECT_DIR / "results" / "aligned.fasta"
DEFAULT_EXCEL = PROJECT_DIR / "results" / "genbank_table.xlsx"
DEFAULT_OUTPUT_DIR = PROJECT_DIR / "results" / "phylogenetics"


def read_sheet(excel_path, sheet_name):
    try:
        return pd.read_excel(excel_path, sheet_name=sheet_name, keep_default_na=False)
    except ValueError:
        return pd.DataFrame()


def accession_from_version(accession_version):
    value = str(accession_version).strip()
    if "." not in value:
        return value
    prefix, suffix = value.rsplit(".", 1)
    if suffix.isdigit():
        return prefix
    return value


def normalize_key(value):
    return str(value).strip()


def collection_month(value):
    text = str(value).strip()
    if not text:
        return ""

    for date_format in ("%d-%b-%Y", "%d-%B-%Y", "%b-%Y", "%B-%Y", "%Y-%m-%d", "%Y"):
        try:
            return datetime.strptime(text, date_format).strftime("%Y-%m")
        except ValueError:
            continue
    return ""


def d614g_accessions(nextclade_mutations):
    if nextclade_mutations.empty:
        return set()

    required = {"gene", "mutation", "accession_version"}
    if not required.issubset(nextclade_mutations.columns):
        return set()

    matches = nextclade_mutations[
        (nextclade_mutations["gene"].astype(str).str.upper() == "S")
        & (nextclade_mutations["mutation"].astype(str).str.upper() == "D614G")
    ]
    return set(matches["accession_version"].astype(str))


def select_columns(df, columns):
    available = [column for column in columns if column in df.columns]
    return df.loc[:, available].copy() if available else pd.DataFrame()


def build_metadata(mapping, excel_path):
    metadata = read_sheet(excel_path, "Metadata_counts")
    nextclade_qc = read_sheet(excel_path, "Nextclade_QC")
    nextclade_mutations = read_sheet(excel_path, "Nextclade_Mutations")
    d614g = d614g_accessions(nextclade_mutations)

    rows = mapping.copy()
    rows["accession_version"] = rows["short_name"].map(normalize_key)
    rows["accession"] = rows["accession_version"].map(accession_from_version)

    if not metadata.empty and "accession_version" in metadata.columns:
        metadata_subset = select_columns(
            metadata,
            [
                "accession_version",
                "country",
                "geo_loc_name",
                "collection_date",
                "host",
                "isolate",
                "length",
                "A_count",
                "G_count",
                "C_count",
                "T_count",
                "N_count",
                "A_%",
                "G_%",
                "C_%",
                "T_%",
                "N_%",
            ],
        )
        rows = rows.merge(metadata_subset, on="accession_version", how="left")

    if not nextclade_qc.empty and "accession_version" in nextclade_qc.columns:
        qc_subset = select_columns(
            nextclade_qc,
            [
                "accession_version",
                "clade",
                "Nextclade_pango",
                "qc.overallStatus",
                "qc.overallScore",
                "coverage",
                "totalSubstitutions",
                "totalAminoacidSubstitutions",
                "warnings",
                "errors",
            ],
        )
        rows = rows.merge(qc_subset, on="accession_version", how="left")

    rows["collection_month"] = rows.get("collection_date", "").map(collection_month)
    rows["has_S_D614G"] = rows["accession_version"].map(lambda value: "yes" if value in d614g else "no")

    preferred_order = [
        "short_name",
        "accession",
        "accession_version",
        "country",
        "geo_loc_name",
        "collection_date",
        "collection_month",
        "clade",
        "Nextclade_pango",
        "has_S_D614G",
        "qc.overallStatus",
        "qc.overallScore",
        "coverage",
        "totalSubstitutions",
        "totalAminoacidSubstitutions",
        "host",
        "isolate",
        "length",
        "A_count",
        "G_count",
        "C_count",
        "T_count",
        "N_count",
        "A_%",
        "G_%",
        "C_%",
        "T_%",
        "N_%",
        "warnings",
        "errors",
        "original_description",
    ]
    ordered = [column for column in preferred_order if column in rows.columns]
    extras = [column for column in rows.columns if column not in ordered]
    return rows.loc[:, ordered + extras].fillna("")


def write_traits(metadata, output_path):
    trait_columns = [
        "short_name",
        "country",
        "collection_month",
        "clade",
        "Nextclade_pango",
        "has_S_D614G",
        "qc.overallStatus",
    ]
    available = [column for column in trait_columns if column in metadata.columns]
    metadata.loc[:, available].to_csv(output_path, sep="\t", index=False)


def write_readme(output_path, alignment_path, metadata_path, traits_path, mapping_path):
    output_path.write_text(
        "\n".join(
            [
                "# Phylogenetic Network Inputs",
                "",
                "These files are prepared for external phylogenetic network tools.",
                "",
                "- `alignment_short.fasta`: aligned sequences with short sample names.",
                "- `network_metadata.tsv`: full metadata table for coloring and annotation.",
                "- `network_traits.tsv`: compact trait table for quick import into network tools.",
                "- `alignment_short_mapping.tsv`: mapping from short names to original FASTA descriptions.",
                "",
                "Recommended workflow:",
                "",
                "1. Use `alignment_short.fasta` as the sequence input.",
                "2. Use `network_traits.tsv` or `network_metadata.tsv` for coloring by country, month, clade, lineage, or D614G.",
                "3. Build the network in PopART, R, or another external tool.",
                "",
                "Generated paths:",
                "",
                f"- alignment: `{alignment_path}`",
                f"- metadata: `{metadata_path}`",
                f"- traits: `{traits_path}`",
                f"- mapping: `{mapping_path}`",
                "",
            ]
        ),
        encoding="utf-8",
    )


def prepare_inputs(aligned_fasta, excel_path, output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    alignment_path = output_dir / "alignment_short.fasta"
    mapping_path = output_dir / "alignment_short_mapping.tsv"
    metadata_path = output_dir / "network_metadata.tsv"
    traits_path = output_dir / "network_traits.tsv"
    readme_path = output_dir / "README_network_inputs.md"

    record_count = create_short_alignment(aligned_fasta, alignment_path, mapping_path)
    mapping = pd.read_csv(mapping_path, sep="\t", keep_default_na=False)
    metadata = build_metadata(mapping, excel_path)
    metadata.to_csv(metadata_path, sep="\t", index=False)
    write_traits(metadata, traits_path)
    write_readme(readme_path, alignment_path, metadata_path, traits_path, mapping_path)

    return {
        "record_count": record_count,
        "alignment": alignment_path,
        "metadata": metadata_path,
        "traits": traits_path,
        "mapping": mapping_path,
        "readme": readme_path,
    }


def parse_args():
    parser = argparse.ArgumentParser(
        description="Prepare aligned FASTA and metadata files for external phylogenetic network tools.",
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
        help="Output directory for network input files.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    outputs = prepare_inputs(args.aligned_fasta, args.excel, args.output_dir)
    print(f"Done: {args.output_dir}")
    print(f"Records processed: {outputs['record_count']}")
    print(f"Alignment: {outputs['alignment']}")
    print(f"Metadata: {outputs['metadata']}")
    print(f"Traits: {outputs['traits']}")


if __name__ == "__main__":
    main()
