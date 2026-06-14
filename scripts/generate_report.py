import argparse
from datetime import datetime
from pathlib import Path

import pandas as pd


PROJECT_DIR = Path(__file__).resolve().parents[1]
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


def markdown_table(df, columns=None, limit=10):
    if df.empty:
        return "_Нет строк для отображения._"

    table = df.copy()
    if columns is not None:
        table = table.loc[:, [column for column in columns if column in table.columns]]
    table = table.head(limit)

    headers = list(table.columns)
    lines = [
        "| " + " | ".join(md(header) for header in headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    for _, row in table.iterrows():
        lines.append("| " + " | ".join(md(row[column]) for column in headers) + " |")
    return "\n".join(lines)


def shortest_longest(metadata):
    columns = ["accession_version", "country", "collection_date", "length"]
    sorted_metadata = metadata.sort_values("length")
    shortest = sorted_metadata.iloc[[0]].loc[:, columns]
    longest = sorted_metadata.iloc[[-1]].loc[:, columns]
    return shortest, longest


def nucleotide_averages(metadata):
    rows = []
    for base in ["A", "G", "C", "T"]:
        count_column = f"{base}_count"
        pct_column = f"{base}_%"
        rows.append({
            "base": base,
            "mean_count": round(metadata[count_column].mean(), 2),
            "mean_percent": round(metadata[pct_column].mean(), 2),
        })
    return pd.DataFrame(rows)


def clade_distribution(nextclade_summary):
    return (
        nextclade_summary[nextclade_summary["section"] == "clade"]
        .loc[:, ["value", "count", "percent"]]
        .rename(columns={"value": "clade"})
        .sort_values(["count", "clade"], ascending=[False, True])
    )


def country_distribution(metadata):
    country_counts = (
        metadata["country"]
        .replace("", "unknown")
        .value_counts()
        .reset_index()
    )
    country_counts.columns = ["country", "sample_count"]
    return country_counts


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


def d614g_status(mutations):
    rows = mutations[
        (mutations["gene"] == "S")
        & (mutations["mutation"] == "D614G")
        & (mutations["mutation_type"] == "amino_acid_substitution")
    ]
    if rows.empty:
        return "Не обнаружена в текущей выборке."

    sample_count = rows["accession_version"].nunique()
    countries = ", ".join(sorted(country for country in rows["country"].replace("", "unknown").unique()))
    return f"Обнаружена: {sample_count} образцов; страны: {countries}."


def build_report(excel_path):
    metadata = read_sheet(excel_path, "Metadata_counts")
    nextclade_summary = read_sheet(excel_path, "Nextclade_Summary")
    gene_summary = read_sheet(excel_path, "Nextclade_Gene_Summary")
    top_mutations = read_sheet(excel_path, "Nextclade_Top_Mutations")
    mutations = read_sheet(excel_path, "Nextclade_Mutations")

    sample_count = len(metadata)
    mean_length = metadata["length"].mean()
    shortest, longest = shortest_longest(metadata)
    nucleotide_means = nucleotide_averages(metadata)
    clades = clade_distribution(nextclade_summary)
    countries = country_distribution(metadata)
    genes = top_genes(gene_summary)
    aa_substitutions = top_amino_acid_substitutions(top_mutations)

    lines = [
        "# Краткий отчет по анализу SARS-CoV-2",
        "",
        f"Сгенерировано: `{datetime.now().astimezone().isoformat(timespec='seconds')}`.",
        f"Источник: `{excel_path}`.",
        "",
        "## Выборка",
        "",
        f"- Образцов в текущем анализе: **{sample_count}**.",
        "- В дипломной работе для сравнения фигурируют 4000 нуклеотидных последовательностей и 8000 последовательностей в филогенетическом блоке.",
        f"- Средняя длина генома в текущей выборке: **{fmt_number(mean_length)} н.**",
        "- Средняя длина генома в выводах диплома: **29870 н.**",
        "",
        "Самый короткий геном:",
        "",
        markdown_table(shortest),
        "",
        "Самый длинный геном:",
        "",
        markdown_table(longest),
        "",
        "## Средний состав A/G/C/T",
        "",
        markdown_table(nucleotide_means),
        "",
        "## Clade в текущей выборке",
        "",
        markdown_table(clades),
        "",
        "## Страны в текущей выборке",
        "",
        markdown_table(countries, limit=20),
        "",
        "## Топ генов по числу изменений",
        "",
        markdown_table(genes, limit=10),
        "",
        "## Топ аминокислотных замен",
        "",
        markdown_table(aa_substitutions, limit=10),
        "",
        "## Проверка D614G",
        "",
        d614g_status(mutations),
        "",
        "## Короткое сравнение с дипломом",
        "",
        "- В дипломе наиболее вариабельным указан ген **S**.",
        f"- В текущей выборке по листу `Nextclade_Gene_Summary` на первом месте: **{md(genes.iloc[0]['gene'])}**.",
        "- Это не противоречие само по себе: текущая выборка намного меньше и имеет другой состав.",
        "- Сильная сторона текущего проекта: отчет обновляется автоматически после пересборки Excel.",
    ]

    return "\n".join(lines) + "\n"


def parse_args():
    parser = argparse.ArgumentParser(description="Generate a concise Markdown report from the pipeline Excel workbook.")
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
