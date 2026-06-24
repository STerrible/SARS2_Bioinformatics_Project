from datetime import datetime
from pathlib import Path

import pandas as pd

from bioseq_pipeline.core.excel import write_dataframes_to_excel
from bioseq_pipeline.organisms.tuberculosis.config import (
    DEFAULT_FASTTREE_LOG,
    DEFAULT_FASTTREE_TREE,
    DEFAULT_GENE_SUMMARY_TSV,
    DEFAULT_INPUT_MANIFEST_TSV,
    DEFAULT_ITOL_COUNTRY_STRIP,
    DEFAULT_ITOL_README,
    DEFAULT_ITOL_TARGET_HEATMAP,
    DEFAULT_ITOL_VARIANT_BARS,
    DEFAULT_MUTATION_REPORT,
    DEFAULT_MUTATION_SUMMARY_XLSX,
    DEFAULT_SAMPLE_SUMMARY_TSV,
    DEFAULT_SNIPPY_CORE_ALN,
    DEFAULT_SNIPPY_CORE_TAB,
    DEFAULT_SNIPPY_CORE_TXT,
    DEFAULT_SNIPPY_RUNS_DIR,
    DEFAULT_SNIPPY_SUMMARY_REPORT,
    DEFAULT_TARGET_MATRIX_TSV,
    DEFAULT_TB_REPORT,
    DEFAULT_VARIANTS_TSV,
    TB_GENE_MUTATION_SUMMARY_SHEET,
    TB_SAMPLE_MUTATION_SUMMARY_SHEET,
    TB_TARGET_GENE_MATRIX_SHEET,
    TB_TARGET_GENE_SUMMARY_SHEET,
    TB_VARIANTS_SHEET,
)


TARGET_GENE_GROUPS = {
    "adaptive_pe_ppe": [
        "PE_PGRS47",
        "PPE18",
        "PE_PGRS33",
        "PE_PGRS62",
        "PE_PGRS11",
    ],
    "drug_resistance_candidate": [
        "katG",
        "rpoB",
        "rpoC",
        "gyrA",
        "gyrB",
        "embB",
        "pncA",
        "inhA",
        "fabG1",
        "rrs",
        "eis",
        "ethA",
        "folC",
        "alr",
        "gid",
    ],
}

TARGET_GENES = [gene for genes in TARGET_GENE_GROUPS.values() for gene in genes]
DRUG_GENE_TO_DRUG = {
    "katG": "isoniazid",
    "inhA": "isoniazid",
    "fabG1": "isoniazid",
    "rpoB": "rifampicin",
    "rpoC": "rifampicin_compensatory",
    "gyrA": "fluoroquinolones",
    "gyrB": "fluoroquinolones",
    "embB": "ethambutol",
    "pncA": "pyrazinamide",
    "rrs": "aminoglycosides",
    "eis": "kanamycin",
    "ethA": "ethionamide",
    "folC": "para_aminosalicylic_acid",
    "alr": "cycloserine",
    "gid": "streptomycin",
}

COUNTRY_COLORS = [
    "#2563eb",
    "#dc2626",
    "#16a34a",
    "#ca8a04",
    "#9333ea",
    "#0891b2",
    "#ea580c",
    "#4f46e5",
    "#be123c",
    "#0f766e",
    "#7c2d12",
    "#475569",
]


def write_text_lf(path, text):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(text)


def safe_read_tsv(path):
    return pd.read_csv(path, sep="\t", keep_default_na=False, dtype=str)


def read_manifest(manifest_path=DEFAULT_INPUT_MANIFEST_TSV):
    manifest_path = Path(manifest_path)
    if not manifest_path.exists():
        raise ValueError(f"TB input manifest was not found: {manifest_path}")
    return safe_read_tsv(manifest_path)


def normalize_gene(row):
    gene = str(row.get("GENE", "")).strip()
    locus_tag = str(row.get("LOCUS_TAG", "")).strip()
    if gene:
        return gene
    if locus_tag:
        return locus_tag
    return "intergenic"


def effect_class(effect):
    effect = str(effect).lower()
    if not effect:
        return "intergenic_or_unannotated"
    if "frameshift" in effect:
        return "frameshift"
    if "stop_gained" in effect:
        return "stop_gained"
    if "missense" in effect:
        return "missense"
    if "synonymous" in effect:
        return "synonymous"
    if "upstream" in effect or "downstream" in effect:
        return "regulatory_or_flanking"
    return "other_annotated"


def gene_group(gene_key, product):
    gene_key = str(gene_key)
    product = str(product).lower()
    if gene_key in TARGET_GENE_GROUPS["adaptive_pe_ppe"]:
        return "adaptive_pe_ppe_target"
    if gene_key in TARGET_GENE_GROUPS["drug_resistance_candidate"]:
        return "drug_resistance_candidate"
    if gene_key == "intergenic":
        return "intergenic"
    if gene_key.startswith(("PE_", "PPE")) or "pe-pgrs" in product or "ppe family" in product:
        return "pe_ppe_family"
    return "other_coding_or_annotated"


def read_snippy_variants(runs_dir=DEFAULT_SNIPPY_RUNS_DIR):
    runs_dir = Path(runs_dir)
    if not runs_dir.exists():
        raise ValueError(f"Snippy runs directory was not found: {runs_dir}")

    frames = []
    for snps_tab in sorted(runs_dir.glob("*/snps.tab")):
        if not snps_tab.exists() or snps_tab.stat().st_size == 0:
            continue
        df = safe_read_tsv(snps_tab)
        df.insert(0, "sample_id", snps_tab.parent.name)
        df.insert(1, "snps_tab_path", str(snps_tab))
        frames.append(df)

    if not frames:
        raise ValueError(f"No non-empty snps.tab files were found under: {runs_dir}")

    variants = pd.concat(frames, ignore_index=True)
    for column in ["CHROM", "POS", "TYPE", "REF", "ALT", "EFFECT", "LOCUS_TAG", "GENE", "PRODUCT"]:
        if column not in variants.columns:
            variants[column] = ""
    variants["POS"] = pd.to_numeric(variants["POS"], errors="coerce").fillna(0).astype("int64")
    variants["gene_key"] = variants.apply(normalize_gene, axis=1)
    variants["effect_class"] = variants["EFFECT"].map(effect_class)
    variants["gene_group"] = variants.apply(lambda row: gene_group(row["gene_key"], row["PRODUCT"]), axis=1)
    variants["drug_or_trait"] = variants["gene_key"].map(lambda gene: DRUG_GENE_TO_DRUG.get(gene, ""))
    variants["is_target_gene"] = variants["gene_key"].isin(TARGET_GENES)
    return variants


def manifest_metadata(manifest_df):
    columns = [
        "sample_id",
        "assembly_accession",
        "assembly_name",
        "country",
        "geo_loc_name",
        "collection_date",
        "host",
        "strain",
        "length",
        "GC_percent",
    ]
    available = [column for column in columns if column in manifest_df.columns]
    return manifest_df.loc[:, available].copy()


def build_gene_summary(variants_df, sample_count):
    grouped = (
        variants_df.groupby(["gene_key", "gene_group", "PRODUCT"], dropna=False)
        .agg(
            variant_rows=("sample_id", "size"),
            samples_with_variant=("sample_id", "nunique"),
            positions=("POS", "nunique"),
            missense_rows=("effect_class", lambda values: int((values == "missense").sum())),
            frameshift_rows=("effect_class", lambda values: int((values == "frameshift").sum())),
            stop_gained_rows=("effect_class", lambda values: int((values == "stop_gained").sum())),
            synonymous_rows=("effect_class", lambda values: int((values == "synonymous").sum())),
        )
        .reset_index()
    )
    grouped["sample_frequency_%"] = grouped["samples_with_variant"].map(lambda value: round(value / sample_count * 100, 2))
    grouped = grouped.sort_values(["samples_with_variant", "variant_rows"], ascending=False)
    return grouped


def build_sample_summary(variants_df, manifest_df):
    summary = (
        variants_df.groupby("sample_id")
        .agg(
            variant_rows=("sample_id", "size"),
            variant_positions=("POS", "nunique"),
            coding_or_annotated_rows=("gene_key", lambda values: int((values != "intergenic").sum())),
            intergenic_rows=("gene_key", lambda values: int((values == "intergenic").sum())),
            missense_rows=("effect_class", lambda values: int((values == "missense").sum())),
            frameshift_rows=("effect_class", lambda values: int((values == "frameshift").sum())),
            stop_gained_rows=("effect_class", lambda values: int((values == "stop_gained").sum())),
            synonymous_rows=("effect_class", lambda values: int((values == "synonymous").sum())),
            pe_ppe_rows=("gene_group", lambda values: int(values.astype(str).str.contains("pe_ppe").sum())),
            target_gene_rows=("is_target_gene", lambda values: int(values.sum())),
            target_genes_mutated=("gene_key", lambda values: ";".join(sorted(set(values) & set(TARGET_GENES)))),
        )
        .reset_index()
    )
    metadata = manifest_metadata(manifest_df)
    if not metadata.empty:
        summary = metadata.merge(summary, on="sample_id", how="right")
    return summary.fillna("")


def build_target_gene_matrix(variants_df, manifest_df):
    samples = pd.DataFrame({"sample_id": sorted(variants_df["sample_id"].unique())})
    metadata = manifest_metadata(manifest_df)
    if not metadata.empty:
        samples = metadata.merge(samples, on="sample_id", how="right")

    target = variants_df[variants_df["gene_key"].isin(TARGET_GENES)]
    counts = target.groupby(["sample_id", "gene_key"]).size().unstack(fill_value=0)
    for gene in TARGET_GENES:
        if gene not in counts.columns:
            counts[gene] = 0
    counts = counts.loc[:, TARGET_GENES].reset_index()
    matrix = samples.merge(counts, on="sample_id", how="left")
    for gene in TARGET_GENES:
        matrix[gene] = pd.to_numeric(matrix[gene], errors="coerce").fillna(0).astype("int64")
    return matrix.fillna("")


def build_target_gene_summary(variants_df, sample_count):
    target = variants_df[variants_df["gene_key"].isin(TARGET_GENES)].copy()
    if target.empty:
        return pd.DataFrame(columns=["gene_key", "target_group", "drug_or_trait", "variant_rows", "samples_with_variant", "sample_frequency_%"])

    def target_group(gene):
        for group, genes in TARGET_GENE_GROUPS.items():
            if gene in genes:
                return group
        return "target"

    summary = (
        target.groupby("gene_key")
        .agg(
            product=("PRODUCT", lambda values: next((value for value in values if value), "")),
            variant_rows=("sample_id", "size"),
            samples_with_variant=("sample_id", "nunique"),
            positions=("POS", "nunique"),
            missense_rows=("effect_class", lambda values: int((values == "missense").sum())),
            frameshift_rows=("effect_class", lambda values: int((values == "frameshift").sum())),
            stop_gained_rows=("effect_class", lambda values: int((values == "stop_gained").sum())),
        )
        .reset_index()
    )
    summary["target_group"] = summary["gene_key"].map(target_group)
    summary["drug_or_trait"] = summary["gene_key"].map(lambda gene: DRUG_GENE_TO_DRUG.get(gene, "adaptive_pe_ppe"))
    summary["sample_frequency_%"] = summary["samples_with_variant"].map(lambda value: round(value / sample_count * 100, 2))
    return summary.sort_values(["target_group", "samples_with_variant", "variant_rows"], ascending=[True, False, False])


def write_mutation_report(report_path, sample_summary, gene_summary, target_summary, variants_df):
    top_genes = gene_summary.head(15)
    target_lines = []
    for _, row in target_summary.iterrows():
        target_lines.append(
            f"| {row['gene_key']} | {row['target_group']} | {row['samples_with_variant']} | {row['sample_frequency_%']} |"
        )

    lines = [
        "# TB assembly-based mutation summary",
        "",
        "This report uses Snippy SNP annotations from assembled genomes. TB-Profiler was not used because raw FASTQ reads are not available in this training workflow.",
        "",
        "## Summary",
        "",
        f"- Generated at: {datetime.now().astimezone().isoformat(timespec='seconds')}",
        f"- Samples: {sample_summary['sample_id'].nunique()}",
        f"- Variant annotation rows: {len(variants_df)}",
        f"- Genes/loci with variants: {gene_summary['gene_key'].nunique()}",
        f"- Target genes tracked: {len(TARGET_GENES)}",
        "",
        "## Target Genes",
        "",
        "| gene | group | samples with variant | sample frequency, % |",
        "| --- | --- | ---: | ---: |",
        *target_lines,
        "",
        "## Top Mutated Genes/Loci",
        "",
        "| gene/locus | group | samples with variant | variant rows |",
        "| --- | --- | ---: | ---: |",
    ]
    for _, row in top_genes.iterrows():
        lines.append(f"| {row['gene_key']} | {row['gene_group']} | {row['samples_with_variant']} | {row['variant_rows']} |")

    lines.extend(
        [
            "",
            "## Methodological Limitation",
            "",
            "Without raw reads this workflow cannot estimate allele fractions, read support, heteroresistance, or formal TB-Profiler resistance classes. Drug-resistance genes are therefore reported as candidate gene mutation summaries, not clinical resistance calls.",
        ]
    )
    write_text_lf(report_path, "\n".join(lines) + "\n")


def write_mutation_summary(
    runs_dir=DEFAULT_SNIPPY_RUNS_DIR,
    manifest_path=DEFAULT_INPUT_MANIFEST_TSV,
    variants_tsv=DEFAULT_VARIANTS_TSV,
    gene_summary_tsv=DEFAULT_GENE_SUMMARY_TSV,
    sample_summary_tsv=DEFAULT_SAMPLE_SUMMARY_TSV,
    target_matrix_tsv=DEFAULT_TARGET_MATRIX_TSV,
    output_xlsx=DEFAULT_MUTATION_SUMMARY_XLSX,
    report_path=DEFAULT_MUTATION_REPORT,
):
    manifest_df = read_manifest(manifest_path)
    variants_df = read_snippy_variants(runs_dir)
    sample_count = variants_df["sample_id"].nunique()
    gene_summary = build_gene_summary(variants_df, sample_count)
    sample_summary = build_sample_summary(variants_df, manifest_df)
    target_matrix = build_target_gene_matrix(variants_df, manifest_df)
    target_summary = build_target_gene_summary(variants_df, sample_count)

    for df, path in [
        (variants_df, variants_tsv),
        (gene_summary, gene_summary_tsv),
        (sample_summary, sample_summary_tsv),
        (target_matrix, target_matrix_tsv),
    ]:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(path, sep="\t", index=False)

    write_dataframes_to_excel(
        output_xlsx,
        [
            (TB_SAMPLE_MUTATION_SUMMARY_SHEET, sample_summary),
            (TB_GENE_MUTATION_SUMMARY_SHEET, gene_summary),
            (TB_TARGET_GENE_SUMMARY_SHEET, target_summary),
            (TB_TARGET_GENE_MATRIX_SHEET, target_matrix),
            (TB_VARIANTS_SHEET, variants_df),
        ],
    )
    write_mutation_report(report_path, sample_summary, gene_summary, target_summary, variants_df)

    return {
        "samples": sample_count,
        "variant_rows": len(variants_df),
        "genes": gene_summary["gene_key"].nunique(),
        "target_genes_with_variants": int((target_summary["samples_with_variant"] > 0).sum()) if not target_summary.empty else 0,
        "variants_tsv": Path(variants_tsv),
        "gene_summary_tsv": Path(gene_summary_tsv),
        "sample_summary_tsv": Path(sample_summary_tsv),
        "target_matrix_tsv": Path(target_matrix_tsv),
        "output_xlsx": Path(output_xlsx),
        "report": Path(report_path),
    }


def country_palette(values):
    labels = sorted({str(value).strip() or "Unknown" for value in values})
    return {label: COUNTRY_COLORS[index % len(COUNTRY_COLORS)] for index, label in enumerate(labels)}


def write_itol_colorstrip(output_path, sample_summary):
    palette = country_palette(sample_summary.get("country", pd.Series(dtype=str)))
    legend_labels = list(palette.keys())
    lines = [
        "DATASET_COLORSTRIP",
        "SEPARATOR TAB",
        "DATASET_LABEL\tCountry",
        "COLOR\t#2563eb",
        f"LEGEND_TITLE\tCountry",
        "LEGEND_SHAPES\t" + "\t".join(["1"] * len(legend_labels)),
        "LEGEND_COLORS\t" + "\t".join(palette[label] for label in legend_labels),
        "LEGEND_LABELS\t" + "\t".join(legend_labels),
        "DATA",
    ]
    for _, row in sample_summary.iterrows():
        country = str(row.get("country", "")).strip() or "Unknown"
        lines.append(f"{row['sample_id']}\t{palette[country]}\t{country}")
    write_text_lf(output_path, "\n".join(lines) + "\n")


def write_itol_simplebar(output_path, sample_summary):
    lines = [
        "DATASET_SIMPLEBAR",
        "SEPARATOR TAB",
        "DATASET_LABEL\tSnippy variant rows",
        "COLOR\t#4f46e5",
        "DATA",
    ]
    for _, row in sample_summary.iterrows():
        lines.append(f"{row['sample_id']}\t{row['variant_rows']}")
    write_text_lf(output_path, "\n".join(lines) + "\n")


def write_itol_target_heatmap(output_path, target_matrix):
    heatmap_genes = TARGET_GENE_GROUPS["adaptive_pe_ppe"] + TARGET_GENE_GROUPS["drug_resistance_candidate"][:7]
    colors = ["#dc2626" if gene in TARGET_GENE_GROUPS["drug_resistance_candidate"] else "#16a34a" for gene in heatmap_genes]
    lines = [
        "DATASET_HEATMAP",
        "SEPARATOR TAB",
        "DATASET_LABEL\tTarget gene mutations",
        "COLOR\t#111827",
        "FIELD_LABELS\t" + "\t".join(heatmap_genes),
        "FIELD_COLORS\t" + "\t".join(colors),
        "DATA",
    ]
    for _, row in target_matrix.iterrows():
        values = ["1" if int(row.get(gene, 0)) > 0 else "0" for gene in heatmap_genes]
        lines.append(f"{row['sample_id']}\t" + "\t".join(values))
    write_text_lf(output_path, "\n".join(lines) + "\n")


def write_itol_readme(readme_path, tree_path, country_strip, variant_bars, target_heatmap):
    lines = [
        "# TB iTOL annotations",
        "",
        "Upload the Newick tree first, then add these dataset files in iTOL.",
        "",
        f"- Tree: `{tree_path}`",
        f"- Country strip: `{country_strip}`",
        f"- Variant count bars: `{variant_bars}`",
        f"- Target gene heatmap: `{target_heatmap}`",
        "",
        "The heatmap is based on assembly-derived Snippy variants, not TB-Profiler resistance calls.",
    ]
    write_text_lf(readme_path, "\n".join(lines) + "\n")


def export_itol_annotations(
    sample_summary_tsv=DEFAULT_SAMPLE_SUMMARY_TSV,
    target_matrix_tsv=DEFAULT_TARGET_MATRIX_TSV,
    tree_path=DEFAULT_FASTTREE_TREE,
    country_strip=DEFAULT_ITOL_COUNTRY_STRIP,
    variant_bars=DEFAULT_ITOL_VARIANT_BARS,
    target_heatmap=DEFAULT_ITOL_TARGET_HEATMAP,
    readme_path=DEFAULT_ITOL_README,
):
    sample_summary = safe_read_tsv(sample_summary_tsv)
    target_matrix = safe_read_tsv(target_matrix_tsv)
    for column in TARGET_GENES:
        if column in target_matrix.columns:
            target_matrix[column] = pd.to_numeric(target_matrix[column], errors="coerce").fillna(0).astype("int64")

    write_itol_colorstrip(country_strip, sample_summary)
    write_itol_simplebar(variant_bars, sample_summary)
    write_itol_target_heatmap(target_heatmap, target_matrix)
    write_itol_readme(readme_path, tree_path, country_strip, variant_bars, target_heatmap)
    return {
        "country_strip": Path(country_strip),
        "variant_bars": Path(variant_bars),
        "target_heatmap": Path(target_heatmap),
        "readme": Path(readme_path),
    }


def metric_lookup_from_markdown_table(path):
    path = Path(path)
    metrics = {}
    metric_aliases = {
        "завершенных Snippy-запусков": "complete_runs",
        "строк вариантов в `core.tab`": "core_tab_variant_rows",
    }
    if not path.exists():
        return metrics
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith("|") or "---" in line or "metric" in line or "показатель" in line:
            continue
        parts = [part.strip() for part in line.strip("|").split("|")]
        if len(parts) >= 2:
            metrics[parts[0]] = parts[1]
            if parts[0] in metric_aliases:
                metrics[metric_aliases[parts[0]]] = parts[1]
    return metrics


def write_final_report(
    report_path=DEFAULT_TB_REPORT,
    manifest_path=DEFAULT_INPUT_MANIFEST_TSV,
    snippy_summary_report=DEFAULT_SNIPPY_SUMMARY_REPORT,
    mutation_report=DEFAULT_MUTATION_REPORT,
    sample_summary_tsv=DEFAULT_SAMPLE_SUMMARY_TSV,
    gene_summary_tsv=DEFAULT_GENE_SUMMARY_TSV,
    target_matrix_tsv=DEFAULT_TARGET_MATRIX_TSV,
    tree_path=DEFAULT_FASTTREE_TREE,
    fasttree_log=DEFAULT_FASTTREE_LOG,
    core_aln=DEFAULT_SNIPPY_CORE_ALN,
    core_tab=DEFAULT_SNIPPY_CORE_TAB,
    core_txt=DEFAULT_SNIPPY_CORE_TXT,
):
    manifest = read_manifest(manifest_path)
    sample_summary = safe_read_tsv(sample_summary_tsv)
    gene_summary = safe_read_tsv(gene_summary_tsv)
    target_matrix = safe_read_tsv(target_matrix_tsv)
    snippy_metrics = metric_lookup_from_markdown_table(snippy_summary_report)

    countries = manifest["country"].replace("", "Unknown").value_counts().head(12) if "country" in manifest.columns else pd.Series(dtype=int)
    top_genes = gene_summary.head(12)
    target_counts = []
    for gene in TARGET_GENES:
        if gene in target_matrix.columns:
            count = int(pd.to_numeric(target_matrix[gene], errors="coerce").fillna(0).gt(0).sum())
            target_counts.append((gene, count, round(count / len(target_matrix) * 100, 2) if len(target_matrix) else 0))

    lines = [
        "# Отчет об assembly-based анализе туберкулеза",
        "",
        f"Сгенерировано: {datetime.now().astimezone().isoformat(timespec='seconds')}",
        "",
        "## Область анализа",
        "",
        "Проект воспроизводит основную вычислительную структуру дипломного workflow без TB-Profiler и без raw FASTQ reads. Используются геномы из NCBI Assembly, аннотации вариантов Snippy, филогенетическое дерево FastTree и файлы аннотаций, готовые для загрузки в iTOL.",
        "",
        "## Входные данные",
        "",
        f"- Сборок в manifest: {len(manifest)}",
        f"- Образцов со сводками мутаций: {sample_summary['sample_id'].nunique()}",
        f"- Core SNP alignment: `{core_aln}`",
        f"- Core SNP matrix: `{core_tab}`",
        f"- Сводка Snippy core: `{core_txt}`",
        "",
        "## Snippy и FastTree",
        "",
        f"- Завершенные Snippy-запуски: {snippy_metrics.get('complete_runs', 'нет данных')}",
        f"- Строк вариантов в core SNP matrix: {snippy_metrics.get('core_tab_variant_rows', 'нет данных')}",
        f"- Newick-дерево FastTree: `{tree_path}`",
        f"- Лог FastTree: `{fasttree_log}`",
        "",
        "## Географический охват",
        "",
        "| страна | образцы |",
        "| --- | ---: |",
    ]
    for country, count in countries.items():
        lines.append(f"| {country} | {count} |")

    lines.extend(["", "## Частоты мутаций в целевых генах", "", "| ген | образцов с вариантом | частота, % |", "| --- | ---: | ---: |"])
    for gene, count, frequency in target_counts:
        lines.append(f"| {gene} | {count} | {frequency} |")

    lines.extend(["", "## Наиболее изменчивые гены и локусы", "", "| ген/локус | группа | образцов с вариантом | строк вариантов |", "| --- | --- | ---: | ---: |"])
    for _, row in top_genes.iterrows():
        lines.append(f"| {row['gene_key']} | {row['gene_group']} | {row['samples_with_variant']} | {row['variant_rows']} |")

    lines.extend(
        [
            "",
            "## Файлы для интерпретации",
            "",
            "- `results/tuberculosis_data/variants/tb_mutation_summary.xlsx` содержит таблицы по образцам, генам, целевым генам и исходным вариантам.",
            "- `results/tuberculosis_data/phylogenetics/tb_fasttree.nwk` готов для загрузки в iTOL.",
            "- `results/tuberculosis_data/itol/` содержит iTOL datasets для страны, числа вариантов и целевых генов.",
            "",
            "## Ограничение",
            "",
            "Так как raw reads отсутствуют, workflow не строит формальные TB-Profiler классы устойчивости, allele fractions, read support и оценки heteroresistance. Гены лекарственной устойчивости суммируются только как candidate mutation loci.",
        ]
    )
    write_text_lf(report_path, "\n".join(lines) + "\n")
    return {
        "report": Path(report_path),
        "samples": sample_summary["sample_id"].nunique(),
        "genes": gene_summary["gene_key"].nunique(),
        "target_genes": len(target_counts),
    }
