import argparse
from datetime import datetime
from pathlib import Path

import pandas as pd


PROJECT_DIR = Path(__file__).resolve().parents[1]
DEFAULT_EXCEL = PROJECT_DIR / "results" / "genbank_table.xlsx"
DEFAULT_OUTPUT = PROJECT_DIR / "reports" / "diploma_comparison.md"

DIPLOMA_RESULTS = {
    "data_scope": "4000 нуклеотидных последовательностей; в разделе филогенетики указаны 8000 последовательностей",
    "average_genome_length": "29870 нуклеотидов",
    "similarity": "85,99%",
    "top_variable_gene": "S",
    "top_variable_gene_note": "в гене S наблюдается наибольшее количество единичных замен (23)",
    "frequent_amino_acids": ["лейцин", "треонин", "гистидин"],
    "ancestor": "NC_045512 (China, 2019/12)",
    "variable_countries": ["Египет", "Нидерланды", "США"],
}


def read_sheet(excel_path, sheet_name, required=True):
    try:
        return pd.read_excel(excel_path, sheet_name=sheet_name, keep_default_na=False)
    except ValueError:
        if required:
            raise
        return pd.DataFrame()


def clean_value(value):
    if pd.isna(value):
        return ""
    return str(value)


def md_escape(value):
    return clean_value(value).replace("|", "\\|").replace("\n", " ")


def markdown_table(df, columns=None, limit=10):
    if df.empty:
        return "_Нет данных._"

    table = df.copy()
    if columns is not None:
        table = table.loc[:, [column for column in columns if column in table.columns]]
    table = table.head(limit)

    headers = list(table.columns)
    lines = [
        "| " + " | ".join(md_escape(header) for header in headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    for _, row in table.iterrows():
        lines.append("| " + " | ".join(md_escape(row[column]) for column in headers) + " |")
    return "\n".join(lines)


def metric_value(run_metadata, metric, default=""):
    if run_metadata.empty or "metric" not in run_metadata.columns or "value" not in run_metadata.columns:
        return default
    rows = run_metadata[run_metadata["metric"] == metric]
    if rows.empty:
        return default
    return clean_value(rows.iloc[0]["value"])


def extract_samples_total(nextclade_summary, metadata_counts):
    if not nextclade_summary.empty:
        rows = nextclade_summary[
            (nextclade_summary["section"] == "samples")
            & (nextclade_summary["value"] == "total")
        ]
        if not rows.empty:
            return int(rows.iloc[0]["count"])
    return len(metadata_counts)


def top_gene(nextclade_gene_summary):
    if nextclade_gene_summary.empty:
        return "", ""
    row = nextclade_gene_summary.iloc[0]
    gene = clean_value(row.get("gene", ""))
    observations = clean_value(row.get("mutation_observations", ""))
    mutation_type = clean_value(row.get("mutation_type", ""))
    return gene, f"{observations} наблюдений; тип: {mutation_type}"


def top_reference_amino_acids(amino_acid_changes):
    if amino_acid_changes.empty:
        return pd.DataFrame()
    rows = amino_acid_changes[amino_acid_changes["role"] == "reference_amino_acid"].copy()
    return rows.sort_values(["mutation_observations", "sample_count"], ascending=[False, False])


def diploma_amino_acid_presence(amino_acid_changes):
    if amino_acid_changes.empty:
        return pd.DataFrame(columns=["amino_acid_name_ru", "present_in_current_data", "mutation_observations", "sample_count"])

    rows = top_reference_amino_acids(amino_acid_changes)
    result = []
    for amino_acid in DIPLOMA_RESULTS["frequent_amino_acids"]:
        hit = rows[rows["amino_acid_name_ru"].str.lower() == amino_acid.lower()]
        if hit.empty:
            result.append({
                "amino_acid_name_ru": amino_acid,
                "present_in_current_data": "нет",
                "mutation_observations": 0,
                "sample_count": 0,
            })
        else:
            row = hit.iloc[0]
            result.append({
                "amino_acid_name_ru": amino_acid,
                "present_in_current_data": "да",
                "mutation_observations": int(row["mutation_observations"]),
                "sample_count": int(row["sample_count"]),
            })
    return pd.DataFrame(result)


def country_mutation_summary(country_mutations):
    if country_mutations.empty:
        return pd.DataFrame()
    return (
        country_mutations
        .groupby("country", dropna=False)
        .agg(
            mutation_observations=("mutation_observations", "sum"),
            unique_mutations=("raw_value", "nunique"),
            country_sample_count=("country_sample_count", "max"),
        )
        .reset_index()
        .sort_values(["mutation_observations", "unique_mutations"], ascending=[False, False])
    )


def diploma_country_presence(country_summary):
    if country_summary.empty:
        return pd.DataFrame(columns=["diploma_country", "present_in_current_data", "current_sample_count"])

    sample_rows = country_summary[country_summary["section"] == "samples"].copy()
    sample_rows["country_lower"] = sample_rows["country"].str.lower()
    aliases = {
        "Египет": ["egypt", "египет"],
        "Нидерланды": ["netherlands", "нидерланды"],
        "США": ["usa", "united states", "сша"],
    }

    result = []
    for country in DIPLOMA_RESULTS["variable_countries"]:
        possible = aliases.get(country, [country.lower()])
        hits = sample_rows[sample_rows["country_lower"].isin(possible)]
        if hits.empty:
            result.append({
                "diploma_country": country,
                "present_in_current_data": "нет",
                "current_sample_count": 0,
            })
        else:
            result.append({
                "diploma_country": country,
                "present_in_current_data": "да",
                "current_sample_count": int(hits.iloc[0]["count"]),
            })
    return pd.DataFrame(result)


def build_report(excel_path):
    metadata_counts = read_sheet(excel_path, "Metadata_counts")
    nextclade_results = read_sheet(excel_path, "Nextclade_results")
    nextclade_summary = read_sheet(excel_path, "Nextclade_Summary", required=False)
    nextclade_gene_summary = read_sheet(excel_path, "Nextclade_Gene_Summary", required=False)
    nextclade_top_mutations = read_sheet(excel_path, "Nextclade_Top_Mutations", required=False)
    country_summary = read_sheet(excel_path, "Country_Summary", required=False)
    country_mutations = read_sheet(excel_path, "Country_Mutations", required=False)
    amino_acid_changes = read_sheet(excel_path, "Amino_Acid_Changes", required=False)
    amino_acid_changes_by_gene = read_sheet(excel_path, "Amino_Acid_Changes_By_Gene", required=False)
    run_metadata = read_sheet(excel_path, "Run_Metadata", required=False)

    sample_count = extract_samples_total(nextclade_summary, metadata_counts)
    current_top_gene, current_top_gene_note = top_gene(nextclade_gene_summary)
    top_aa = top_reference_amino_acids(amino_acid_changes)
    aa_presence = diploma_amino_acid_presence(amino_acid_changes)
    country_top = country_mutation_summary(country_mutations)
    country_presence = diploma_country_presence(country_summary)

    generated_at = datetime.now().astimezone().isoformat(timespec="seconds")
    nextclade_version = metric_value(run_metadata, "nextclade_version", "не указано")
    dataset = metric_value(run_metadata, "nextclade_dataset", "не указано")
    run_timestamp = metric_value(run_metadata, "run_timestamp", "не указано")

    lines = [
        "# Сравнение автоматизированного анализа с выводами дипломной работы",
        "",
        f"Отчет сгенерирован: `{generated_at}`.",
        "",
        "## Назначение отчета",
        "",
        "Этот отчет сопоставляет результаты текущего автоматизированного pipeline с контрольными выводами исходной дипломной работы. "
        "Сравнение носит ориентировочный характер: размер и состав текущей выборки могут отличаться от выборки диплома.",
        "",
        "## Источники данных",
        "",
        f"- Итоговый Excel: `{excel_path}`",
        f"- Время запуска pipeline по `Run_Metadata`: `{run_timestamp}`",
        f"- Версия Nextclade: `{nextclade_version}`",
        f"- Датасет Nextclade: `{dataset}`",
        f"- Мутации, clade, lineage и QC берутся только из листа `Nextclade_results`.",
        "",
        "## Масштаб выборки",
        "",
        "| Показатель | Диплом | Текущий pipeline |",
        "| --- | --- | --- |",
        f"| Объем данных | {md_escape(DIPLOMA_RESULTS['data_scope'])} | {sample_count} образцов |",
        f"| Средняя длина генома | {md_escape(DIPLOMA_RESULTS['average_genome_length'])} | см. `Metadata_counts` |",
        f"| Процент сходства | {md_escape(DIPLOMA_RESULTS['similarity'])} | не пересчитывается; проект не делает собственное множественное выравнивание |",
        "",
        "## Сравнение ключевых выводов",
        "",
        "| Тезис диплома | Текущий результат | Комментарий |",
        "| --- | --- | --- |",
        (
            f"| Наиболее вариабельный ген: `{DIPLOMA_RESULTS['top_variable_gene']}` | "
            f"`{md_escape(current_top_gene)}` ({md_escape(current_top_gene_note)}) | "
            "Текущий результат рассчитан по готовым аминокислотным изменениям Nextclade и зависит от состава выборки. |"
        ),
        (
            f"| Частые аминокислоты: {', '.join(DIPLOMA_RESULTS['frequent_amino_acids'])} | "
            "см. таблицу ниже | Сравнение выполнено по листу `Amino_Acid_Changes`. |"
        ),
        (
            f"| Предковый штамм: `{DIPLOMA_RESULTS['ancestor']}` | "
            "не определяется | Pipeline не строит филогенетическое дерево и не выявляет предковый штамм. |"
        ),
        (
            f"| Наиболее вариабельные страны: {', '.join(DIPLOMA_RESULTS['variable_countries'])} | "
            "см. таблицу ниже | Текущая выборка меньше и может не содержать эти страны. |"
        ),
        "",
        "## Топ генов по аминокислотным изменениям",
        "",
        markdown_table(nextclade_gene_summary, limit=12),
        "",
        "## Топ мутаций по текущей выборке",
        "",
        markdown_table(nextclade_top_mutations, limit=15),
        "",
        "## Проверка аминокислот из выводов диплома",
        "",
        markdown_table(aa_presence, limit=10),
        "",
        "## Частоты аминокислот в текущей выборке",
        "",
        markdown_table(top_aa, limit=15),
        "",
        "## Частоты аминокислот по gene",
        "",
        markdown_table(amino_acid_changes_by_gene, limit=20),
        "",
        "## Страны из выводов диплома в текущей выборке",
        "",
        markdown_table(country_presence, limit=10),
        "",
        "## Страны с наибольшим числом наблюдений мутаций в текущей выборке",
        "",
        markdown_table(country_top, limit=15),
        "",
        "## QC и распределение clade/lineage",
        "",
        markdown_table(nextclade_summary, limit=20),
        "",
        "## Методические ограничения",
        "",
        "- Текущая выборка существенно меньше выборки диплома, поэтому совпадения и расхождения нельзя трактовать как окончательный биологический вывод.",
        "- Pipeline не выполняет собственное множественное выравнивание, pairwise alignment или расчет p-distance.",
        "- Pipeline не загружает данные из NCBI и анализирует только локально подготовленные `sequence.gb` и `sequence.fasta`.",
        "- Мутации и QC берутся из Nextclade; листы отчета только агрегируют уже готовые результаты.",
        "- Филогенетика, определение предкового штамма и Pangolin не входят в текущий pipeline.",
        "",
        "## Автоматический вывод",
        "",
        (
            f"Текущий pipeline успешно автоматизирует табличную часть анализа для {sample_count} образцов: "
            "извлекает GenBank-метаданные, запускает Nextclade, сохраняет raw-результаты и строит сводные листы. "
            "По сравнению с дипломной работой наиболее сильная сторона проекта - воспроизводимость и скорость обновления отчета. "
            "Главное ограничение - меньший размер выборки и отсутствие филогенетического блока."
        ),
    ]

    return "\n".join(lines) + "\n"


def parse_args():
    parser = argparse.ArgumentParser(description="Generate Markdown report comparing pipeline results with diploma conclusions.")
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
