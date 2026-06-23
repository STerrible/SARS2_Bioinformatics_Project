import argparse
from datetime import datetime
from pathlib import Path

import pandas as pd


PROJECT_DIR = Path(__file__).resolve().parents[4]
DEFAULT_EXCEL = PROJECT_DIR / "results" / "genbank_table.xlsx"
DEFAULT_OUTPUT = PROJECT_DIR / "reports" / "sars2_analysis_report.md"

COLUMN_LABELS = {
    "metric": "Показатель",
    "value": "Значение",
    "field": "Поле",
    "filled": "Заполнено",
    "missing": "Пусто",
    "filled_percent": "Заполнено, %",
    "count": "Количество",
    "percent": "Доля, %",
    "country": "Страна",
    "sample_count": "Образцов",
    "base": "Нуклеотид",
    "mean_count": "Среднее число",
    "mean_percent": "Средняя доля, %",
    "gene": "Ген",
    "mutation_type": "Тип мутации",
    "mutation_observations": "Наблюдений",
    "mutation_label": "Мутация",
    "percent_of_samples": "Доля образцов, %",
    "role": "Роль",
    "amino_acid": "Аминокислота",
    "amino_acid_name": "Название",
}

PERCENT_COLUMNS = {
    "filled_percent",
    "mean_percent",
    "percent",
    "percent_of_samples",
}


def read_sheet(excel_path, sheet_name):
    return pd.read_excel(excel_path, sheet_name=sheet_name, keep_default_na=False)


def text(value):
    if pd.isna(value):
        return ""
    return str(value)


def md(value):
    return text(value).replace("|", "\\|").replace("\n", " ")


def fmt_number(value, digits=2):
    if pd.isna(value):
        return ""
    return f"{float(value):.{digits}f}"


def fmt_percent(value):
    return f"{float(value):.2f}%"


def fmt_coverage(value):
    if pd.isna(value):
        return ""
    value = float(value)
    if abs(value) <= 1.5:
        value *= 100
    return fmt_percent(value)


def sample_percent(count, total):
    return round(count / total * 100, 2) if total else 0


def format_table_cell(column, value):
    if column in PERCENT_COLUMNS and text(value) != "":
        return fmt_percent(value)
    return value


def markdown_table(df, columns=None, limit=None):
    if df.empty:
        return "_Нет данных._"

    table = df.copy()
    if columns is not None:
        table = table.loc[:, [column for column in columns if column in table.columns]]
    if limit is not None:
        table = table.head(limit)

    headers = list(table.columns)
    display_headers = [COLUMN_LABELS.get(header, header) for header in headers]
    lines = [
        "| " + " | ".join(md(header) for header in display_headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    for _, row in table.iterrows():
        lines.append("| " + " | ".join(md(format_table_cell(column, row[column])) for column in headers) + " |")
    return "\n".join(lines)


def dataframe_from_rows(rows, columns=("metric", "value")):
    return pd.DataFrame(rows, columns=list(columns))


def first_distribution_value(summary_df, section):
    distribution = summary_distribution(summary_df, section)
    if distribution.empty:
        return "", 0, 0
    row = distribution.iloc[0]
    return row["value"], int(row["count"]), float(row["percent"])


def summary_distribution(summary_df, section):
    if summary_df.empty or "section" not in summary_df.columns:
        return pd.DataFrame(columns=["value", "count", "percent"])

    rows = summary_df[summary_df["section"] == section].copy()
    if rows.empty:
        return pd.DataFrame(columns=["value", "count", "percent"])

    rows = rows.loc[:, ["value", "count", "percent"]]
    rows["count"] = pd.to_numeric(rows["count"], errors="coerce").fillna(0).astype(int)
    rows["percent"] = pd.to_numeric(rows["percent"], errors="coerce").fillna(0).round(2)
    return rows.sort_values(["count", "value"], ascending=[False, True]).reset_index(drop=True)


def country_distribution(metadata):
    total = len(metadata)
    country_counts = metadata["country"].replace("", "unknown").value_counts(dropna=False).reset_index()
    country_counts.columns = ["country", "sample_count"]
    country_counts["percent"] = country_counts["sample_count"].map(lambda count: sample_percent(count, total))
    return country_counts


def metadata_completeness(metadata):
    rows = []
    for column, label in [
        ("country", "Страна"),
        ("collection_date", "Дата сбора"),
        ("host", "Хозяин"),
        ("isolate", "Изолят"),
    ]:
        if column not in metadata.columns:
            continue
        filled = metadata[column].astype(str).str.strip().ne("").sum()
        rows.append({
            "field": label,
            "filled": int(filled),
            "missing": int(len(metadata) - filled),
            "filled_percent": sample_percent(int(filled), len(metadata)),
        })
    return pd.DataFrame(rows)


def parse_collection_dates(metadata):
    if "collection_date" not in metadata.columns:
        return pd.Series(dtype="datetime64[ns]")
    return pd.to_datetime(metadata["collection_date"].replace("", pd.NA), errors="coerce")


def genome_length_summary(metadata):
    lengths = pd.to_numeric(metadata["length"], errors="coerce").dropna()
    shortest = metadata.loc[pd.to_numeric(metadata["length"], errors="coerce").idxmin()]
    longest = metadata.loc[pd.to_numeric(metadata["length"], errors="coerce").idxmax()]
    gc_percent = pd.to_numeric(metadata["G_%"], errors="coerce") + pd.to_numeric(metadata["C_%"], errors="coerce")
    n_counts = pd.to_numeric(metadata["N_count"], errors="coerce").fillna(0)

    return dataframe_from_rows([
        ("Количество последовательностей", len(metadata)),
        ("Средняя длина генома", fmt_number(lengths.mean())),
        ("Медианная длина генома", fmt_number(lengths.median())),
        ("Минимальная длина", f"{int(shortest['length'])} ({shortest['accession_version']})"),
        ("Максимальная длина", f"{int(longest['length'])} ({longest['accession_version']})"),
        ("Средняя доля GC", fmt_percent(gc_percent.mean())),
        ("Последовательностей с N", int((n_counts > 0).sum())),
        ("Среднее число N на последовательность", fmt_number(n_counts.mean())),
    ])


def nucleotide_averages(metadata):
    rows = []
    for base in ["A", "G", "C", "T", "N"]:
        count_column = f"{base}_count"
        percent_column = f"{base}_%"
        if count_column not in metadata.columns or percent_column not in metadata.columns:
            continue
        rows.append({
            "base": base,
            "mean_count": round(pd.to_numeric(metadata[count_column], errors="coerce").mean(), 2),
            "mean_percent": round(pd.to_numeric(metadata[percent_column], errors="coerce").mean(), 2),
        })
    return pd.DataFrame(rows)


def qc_overview(nextclade_qc):
    if nextclade_qc.empty:
        return dataframe_from_rows([])

    coverage = pd.to_numeric(nextclade_qc.get("coverage", pd.Series(dtype=float)), errors="coerce")
    warnings = nextclade_qc.get("warnings", pd.Series([""] * len(nextclade_qc))).astype(str).str.strip()
    errors = nextclade_qc.get("errors", pd.Series([""] * len(nextclade_qc))).astype(str).str.strip()
    missing = pd.to_numeric(nextclade_qc.get("totalMissing", pd.Series(dtype=float)), errors="coerce")

    return dataframe_from_rows([
        ("Средний coverage", fmt_coverage(coverage.mean())),
        ("Минимальный coverage", fmt_coverage(coverage.min())),
        ("Максимальный coverage", fmt_coverage(coverage.max())),
        ("Среднее число missing bases", fmt_number(missing.mean())),
        ("Последовательностей с warnings", int(warnings.ne("").sum())),
        ("Последовательностей с errors", int(errors.ne("").sum())),
    ])


def top_genes(gene_summary):
    if gene_summary.empty:
        return pd.DataFrame()
    return (
        gene_summary
        .sort_values(["mutation_observations", "sample_count", "gene"], ascending=[False, False, True])
        .loc[:, ["gene", "mutation_type", "mutation_observations", "sample_count"]]
        .reset_index(drop=True)
    )


def top_mutation_table(top_mutations, mutation_type=None, limit=15):
    if top_mutations.empty:
        return pd.DataFrame()

    rows = top_mutations.copy()
    if mutation_type is not None:
        rows = rows[rows["mutation_type"] == mutation_type]
    if rows.empty:
        return pd.DataFrame()

    rows["mutation_label"] = rows.apply(
        lambda row: f"{row['gene']}:{row['mutation']}" if str(row.get("gene", "")).strip() else row["mutation"],
        axis=1,
    )
    return (
        rows
        .sort_values(["sample_count", "mutation_observations", "mutation_label"], ascending=[False, False, True])
        .loc[:, ["mutation_type", "mutation_label", "sample_count", "percent_of_samples", "mutation_observations"]]
        .head(limit)
        .reset_index(drop=True)
    )


def amino_acid_change_summary(amino_acid_changes):
    if amino_acid_changes.empty:
        return pd.DataFrame()

    rows = amino_acid_changes.copy()
    rows = rows[rows["role"].isin(["reference_amino_acid", "sample_amino_acid"])]
    if rows.empty:
        return pd.DataFrame()

    rows["role"] = rows["role"].replace({
        "reference_amino_acid": "в референсе",
        "sample_amino_acid": "в образцах",
    })
    return (
        rows
        .sort_values(["role", "mutation_observations", "sample_count", "amino_acid"], ascending=[True, False, False, True])
        .loc[:, ["role", "amino_acid", "amino_acid_name", "mutation_observations", "sample_count"]]
        .head(20)
        .reset_index(drop=True)
    )


def d614g_summary(mutations, total_samples):
    rows = mutations[
        (mutations["gene"] == "S")
        & (mutations["mutation"] == "D614G")
        & (mutations["mutation_type"] == "amino_acid_substitution")
    ].copy()
    if rows.empty:
        return dataframe_from_rows([
            ("Мутация", "S:D614G"),
            ("Образцов", 0),
            ("Доля образцов", "0.00%"),
        ])

    sample_count = rows["accession_version"].nunique()
    countries = (
        rows.assign(country=rows["country"].replace("", "unknown"))
        .groupby("country")["accession_version"]
        .nunique()
        .sort_values(ascending=False)
    )
    top_countries = ", ".join(f"{country} ({count})" for country, count in countries.head(12).items())
    return dataframe_from_rows([
        ("Мутация", "S:D614G"),
        ("Образцов", sample_count),
        ("Доля образцов", fmt_percent(sample_percent(sample_count, total_samples))),
        ("Стран с мутацией", len(countries)),
        ("Топ стран", top_countries),
    ])


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
        ("nextclade_dataset_name_requested", "Запрошенный dataset"),
        ("nextclade_dataset_tag_requested", "Запрошенный dataset tag"),
        ("nextclade_dataset_version_tag", "Фактический dataset tag"),
        ("nextclade_dataset_updated_at", "Dataset обновлен"),
        ("nextclade_dataset_reference_name", "Референс"),
        ("nextclade_dataset_reference_accession", "Accession референса"),
        ("nextclade_dataset_cli_compatibility", "Совместимость CLI"),
        ("nextclade_dataset_metadata_status", "Статус metadata"),
        ("nextclade_dataset", "Путь к dataset"),
    ]
    rows = []
    for metric, label in metrics:
        value = run_metadata_value(run_metadata, metric)
        if value:
            rows.append({"metric": label, "value": value})
    return pd.DataFrame(rows, columns=["metric", "value"])


def report_highlights(metadata, nextclade_summary, top_mutations, mutations):
    total_samples = len(metadata)
    countries = metadata["country"].replace("", pd.NA).dropna().nunique()
    dates = parse_collection_dates(metadata).dropna()
    date_range = ""
    if not dates.empty:
        date_range = f"{dates.min().date()} - {dates.max().date()}"

    top_clade, top_clade_count, top_clade_percent = first_distribution_value(nextclade_summary, "clade")
    top_lineage, top_lineage_count, top_lineage_percent = first_distribution_value(nextclade_summary, "Nextclade_pango")
    top_mutations_table = top_mutation_table(top_mutations, limit=1)
    top_mutation = ""
    if not top_mutations_table.empty:
        row = top_mutations_table.iloc[0]
        top_mutation = f"{row['mutation_label']} ({int(row['sample_count'])}; {fmt_percent(row['percent_of_samples'])})"

    d614g_rows = mutations[
        (mutations["gene"] == "S")
        & (mutations["mutation"] == "D614G")
        & (mutations["mutation_type"] == "amino_acid_substitution")
    ]
    d614g_count = d614g_rows["accession_version"].nunique() if not d614g_rows.empty else 0

    return dataframe_from_rows([
        ("Образцов в анализе", total_samples),
        ("Стран/регионов с указанной страной", countries),
        ("Диапазон дат сбора", date_range or "не определен"),
        ("Самый частый clade", f"{top_clade} ({top_clade_count}; {fmt_percent(top_clade_percent)})" if top_clade else ""),
        ("Самый частый lineage", f"{top_lineage} ({top_lineage_count}; {fmt_percent(top_lineage_percent)})" if top_lineage else ""),
        ("Самая частая мутация", top_mutation),
        ("S:D614G", f"{d614g_count} образцов; {fmt_percent(sample_percent(d614g_count, total_samples))}"),
    ])


def build_report(excel_path):
    metadata = read_sheet(excel_path, "Metadata_counts")
    nextclade_summary = read_sheet(excel_path, "Nextclade_Summary")
    nextclade_qc = read_sheet(excel_path, "Nextclade_QC")
    gene_summary = read_sheet(excel_path, "Nextclade_Gene_Summary")
    top_mutations = read_sheet(excel_path, "Nextclade_Top_Mutations")
    mutations = read_sheet(excel_path, "Nextclade_Mutations")
    amino_acid_changes = read_sheet(excel_path, "Amino_Acid_Changes")
    run_metadata = read_sheet(excel_path, "Run_Metadata")

    total_countries = metadata["country"].replace("", pd.NA).dropna().nunique()

    sections = [
        "# Отчет по анализу SARS-CoV-2",
        "",
        f"`generated_at`: `{datetime.now().astimezone().isoformat(timespec='seconds')}`",
        f"`source_excel`: `{excel_path}`",
        "",
        "## Краткая сводка",
        "",
        markdown_table(report_highlights(metadata, nextclade_summary, top_mutations, mutations)),
        "",
        "## Качество и полнота данных",
        "",
        "### Полнота метаданных",
        "",
        markdown_table(metadata_completeness(metadata)),
        "",
        "### Nextclade QC",
        "",
        markdown_table(qc_overview(nextclade_qc)),
        "",
        "### Распределение QC-статусов",
        "",
        markdown_table(summary_distribution(nextclade_summary, "qc.overallStatus")),
        "",
        "## Геномы и нуклеотидный состав",
        "",
        markdown_table(genome_length_summary(metadata)),
        "",
        "### Средний нуклеотидный состав",
        "",
        markdown_table(nucleotide_averages(metadata)),
        "",
        "## География выборки",
        "",
        f"Всего стран/регионов с указанной страной: **{total_countries}**. Ниже показаны первые 30 по числу образцов.",
        "",
        markdown_table(country_distribution(metadata), limit=30),
        "",
        "## Clade и lineage",
        "",
        "### Clade",
        "",
        markdown_table(summary_distribution(nextclade_summary, "clade"), limit=20),
        "",
        "### Nextclade lineage / Pango",
        "",
        markdown_table(summary_distribution(nextclade_summary, "Nextclade_pango"), limit=20),
        "",
        "## Мутации",
        "",
        "### Гены с наибольшим числом изменений",
        "",
        markdown_table(top_genes(gene_summary), limit=15),
        "",
        "### Самые частые мутации",
        "",
        markdown_table(top_mutation_table(top_mutations, limit=15)),
        "",
        "### Самые частые аминокислотные замены",
        "",
        markdown_table(top_mutation_table(top_mutations, mutation_type="amino_acid_substitution", limit=15)),
        "",
        "### S:D614G",
        "",
        markdown_table(d614g_summary(mutations, len(metadata))),
        "",
        "### Аминокислоты, чаще всего участвующие в заменах",
        "",
        markdown_table(amino_acid_change_summary(amino_acid_changes)),
        "",
        "## Воспроизводимость",
        "",
        markdown_table(reproducibility_table(run_metadata)),
        "",
        "## Ограничения интерпретации",
        "",
        "- Отчет использует уже готовые результаты Nextclade и не пересчитывает мутации самостоятельно.",
        "- Наличие страны, даты сбора, host и isolate зависит от полноты исходных GenBank-записей.",
        "- Распределения отражают только текущий локальный набор последовательностей, а не глобальную частоту вариантов.",
        "- Филогенетическая сеть и деревья строятся внешними инструментами, не этим отчетом.",
    ]

    return "\n".join(sections) + "\n"


def parse_args():
    parser = argparse.ArgumentParser(description="Generate a SARS-CoV-2 Markdown analysis report.")
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
