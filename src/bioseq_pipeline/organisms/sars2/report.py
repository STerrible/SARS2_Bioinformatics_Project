import argparse
from datetime import datetime
from pathlib import Path

import pandas as pd


PROJECT_DIR = Path(__file__).resolve().parents[4]
DEFAULT_EXCEL = PROJECT_DIR / "results" / "genbank_table.xlsx"
DEFAULT_OUTPUT = PROJECT_DIR / "reports" / "diploma_comparison.md"


def read_sheet(excel_path, sheet_name):
    return pd.read_excel(excel_path, sheet_name=sheet_name, keep_default_na=False)


def text(value):
    if pd.isna(value):
        return ""
    return str(value)


def md(value):
    return text(value).replace("|", "\\|").replace("\n", " ")


def fmt_number(value, digits=2):
    return f"{float(value):.{digits}f}"


def markdown_table(df, columns=None, limit=None):
    if df.empty:
        return ""

    table = df.copy()
    if columns is not None:
        table = table.loc[:, [column for column in columns if column in table.columns]]
    if limit is not None:
        table = table.head(limit)

    headers = list(table.columns)
    lines = [
        "| " + " | ".join(md(header) for header in headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    for _, row in table.iterrows():
        lines.append("| " + " | ".join(md(row[column]) for column in headers) + " |")
    return "\n".join(lines)


def country_distribution(metadata):
    country_counts = metadata["country"].replace("", "unknown").value_counts().reset_index()
    country_counts.columns = ["country", "sample_count"]
    return country_counts


def clade_distribution(nextclade_summary):
    return (
        nextclade_summary[nextclade_summary["section"] == "clade"]
        .loc[:, ["value", "count", "percent"]]
        .rename(columns={"value": "clade"})
        .sort_values(["count", "clade"], ascending=[False, True])
    )


def top_genes(gene_summary):
    return (
        gene_summary
        .sort_values(["mutation_observations", "sample_count", "gene"], ascending=[False, False, True])
        .loc[:, ["gene", "mutation_type", "mutation_observations", "sample_count"]]
    )


def top_amino_acid_substitutions(top_mutations):
    return (
        top_mutations[top_mutations["mutation_type"] == "amino_acid_substitution"]
        .sort_values(["sample_count", "mutation_observations", "gene", "mutation"], ascending=[False, False, True, True])
        .loc[:, ["gene", "mutation", "sample_count", "percent_of_samples"]]
    )


def nucleotide_averages(metadata):
    rows = []
    for base in ["A", "G", "C", "T"]:
        rows.append({
            "base": base,
            "mean_count": round(metadata[f"{base}_count"].mean(), 2),
            "mean_percent": round(metadata[f"{base}_%"].mean(), 2),
        })
    return pd.DataFrame(rows)


def genome_length_summary(metadata):
    shortest = metadata.loc[metadata["length"].idxmin()]
    longest = metadata.loc[metadata["length"].idxmax()]
    return pd.DataFrame([
        {
            "metric": "sample_count",
            "value": len(metadata),
            "accession_version": "",
            "country": "",
            "collection_date": "",
        },
        {
            "metric": "mean_length",
            "value": fmt_number(metadata["length"].mean()),
            "accession_version": "",
            "country": "",
            "collection_date": "",
        },
        {
            "metric": "shortest_genome",
            "value": int(shortest["length"]),
            "accession_version": shortest["accession_version"],
            "country": shortest["country"],
            "collection_date": shortest["collection_date"],
        },
        {
            "metric": "longest_genome",
            "value": int(longest["length"]),
            "accession_version": longest["accession_version"],
            "country": longest["country"],
            "collection_date": longest["collection_date"],
        },
    ])


def d614g_table(mutations):
    rows = mutations[
        (mutations["gene"] == "S")
        & (mutations["mutation"] == "D614G")
        & (mutations["mutation_type"] == "amino_acid_substitution")
    ]
    if rows.empty:
        return pd.DataFrame([{
            "mutation": "S:D614G",
            "sample_count": 0,
            "countries": "",
        }])

    return pd.DataFrame([{
        "mutation": "S:D614G",
        "sample_count": rows["accession_version"].nunique(),
        "countries": ", ".join(sorted(rows["country"].replace("", "unknown").unique())),
    }])


def run_metadata_value(run_metadata, metric):
    if run_metadata.empty or "metric" not in run_metadata.columns or "value" not in run_metadata.columns:
        return ""
    rows = run_metadata[run_metadata["metric"] == metric]
    if rows.empty:
        return ""
    return text(rows.iloc[0]["value"])


def reproducibility_table(run_metadata):
    metrics = [
        ("container_image", "Docker image"),
        ("nextclade_version", "Nextclade CLI"),
        ("nextclade_dataset_name_requested", "Requested dataset"),
        ("nextclade_dataset_tag_requested", "Requested dataset tag"),
        ("nextclade_dataset_version_tag", "Actual dataset tag"),
        ("nextclade_dataset_updated_at", "Dataset updated at"),
        ("nextclade_dataset_reference_name", "Dataset reference"),
        ("nextclade_dataset_reference_accession", "Dataset reference accession"),
        ("nextclade_dataset_cli_compatibility", "Dataset CLI compatibility"),
        ("nextclade_dataset_metadata_status", "Dataset metadata status"),
        ("nextclade_dataset", "Dataset path"),
    ]
    rows = []
    for metric, label in metrics:
        value = run_metadata_value(run_metadata, metric)
        if value:
            rows.append({"metric": label, "value": value})
    return pd.DataFrame(rows, columns=["metric", "value"])


def diploma_comparison(metadata, gene_summary, amino_acid_changes):
    top_gene = top_genes(gene_summary).iloc[0]
    aa_rows = amino_acid_changes[amino_acid_changes["role"] == "reference_amino_acid"].copy()

    def aa_count(name):
        rows = aa_rows[aa_rows["amino_acid_name_ru"].str.lower() == name]
        if rows.empty:
            return 0
        return int(rows.iloc[0]["mutation_observations"])

    return pd.DataFrame([
        {
            "metric": "sample_count",
            "diploma_value": "4000; phylogeny: 8000",
            "pipeline_value": len(metadata),
        },
        {
            "metric": "mean_genome_length",
            "diploma_value": "29870",
            "pipeline_value": fmt_number(metadata["length"].mean()),
        },
        {
            "metric": "top_variable_gene",
            "diploma_value": "S",
            "pipeline_value": f"{top_gene['gene']} ({top_gene['mutation_observations']})",
        },
        {
            "metric": "leucine_observations",
            "diploma_value": "frequent",
            "pipeline_value": aa_count("лейцин"),
        },
        {
            "metric": "threonine_observations",
            "diploma_value": "frequent",
            "pipeline_value": aa_count("треонин"),
        },
        {
            "metric": "histidine_observations",
            "diploma_value": "frequent",
            "pipeline_value": aa_count("гистидин"),
        },
        {
            "metric": "diploma_variable_countries_present",
            "diploma_value": "Egypt; Netherlands; USA",
            "pipeline_value": "; ".join(
                f"{country}:{'yes' if country.lower() in set(metadata['country'].str.lower()) else 'no'}"
                for country in ["Egypt", "Netherlands", "USA"]
            ),
        },
    ])


def build_report(excel_path):
    metadata = read_sheet(excel_path, "Metadata_counts")
    nextclade_summary = read_sheet(excel_path, "Nextclade_Summary")
    gene_summary = read_sheet(excel_path, "Nextclade_Gene_Summary")
    top_mutations = read_sheet(excel_path, "Nextclade_Top_Mutations")
    mutations = read_sheet(excel_path, "Nextclade_Mutations")
    amino_acid_changes = read_sheet(excel_path, "Amino_Acid_Changes")
    run_metadata = read_sheet(excel_path, "Run_Metadata")

    sections = [
        "# Автоматический сравнительный отчет",
        "",
        f"`generated_at`: `{datetime.now().astimezone().isoformat(timespec='seconds')}`",
        f"`source_excel`: `{excel_path}`",
        "",
        "## Reproducibility",
        "",
        markdown_table(reproducibility_table(run_metadata)),
        "",
        "## Diploma_Comparison",
        "",
        markdown_table(diploma_comparison(metadata, gene_summary, amino_acid_changes)),
        "",
        "## Genome_Length_Summary",
        "",
        markdown_table(genome_length_summary(metadata)),
        "",
        "## Mean_Nucleotide_Composition",
        "",
        markdown_table(nucleotide_averages(metadata)),
        "",
        "## Clade_Distribution",
        "",
        markdown_table(clade_distribution(nextclade_summary)),
        "",
        "## Country_Distribution",
        "",
        markdown_table(country_distribution(metadata)),
        "",
        "## Top_Genes_By_Changes",
        "",
        markdown_table(top_genes(gene_summary), limit=10),
        "",
        "## Top_Amino_Acid_Substitutions",
        "",
        markdown_table(top_amino_acid_substitutions(top_mutations), limit=10),
        "",
        "## D614G",
        "",
        markdown_table(d614g_table(mutations)),
    ]

    return "\n".join(sections) + "\n"


def parse_args():
    parser = argparse.ArgumentParser(description="Generate an automatic Markdown comparison report.")
    parser.add_argument("--input", type=Path, default=DEFAULT_EXCEL, help="Input Excel workbook generated by the pipeline.")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT, help="Output Markdown report path.")
    return parser.parse_args()


def main():
    args = parse_args()
    if not args.input.exists():
        raise SystemExit(f"Input Excel workbook was not found: {args.input}")

    report = build_report(args.input)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(report, encoding="utf-8")
    print(f"Report written: {args.output}")


if __name__ == "__main__":
    main()
