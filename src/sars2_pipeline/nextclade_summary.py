import os

import pandas as pd


AMINO_ACID_NAMES = {
    "A": "alanine",
    "R": "arginine",
    "N": "asparagine",
    "D": "aspartic acid",
    "C": "cysteine",
    "Q": "glutamine",
    "E": "glutamic acid",
    "G": "glycine",
    "H": "histidine",
    "I": "isoleucine",
    "L": "leucine",
    "K": "lysine",
    "M": "methionine",
    "F": "phenylalanine",
    "P": "proline",
    "S": "serine",
    "T": "threonine",
    "W": "tryptophan",
    "Y": "tyrosine",
    "V": "valine",
    "*": "stop",
    "-": "deletion",
}

AMINO_ACID_NAMES_RU = {
    "A": "аланин",
    "R": "аргинин",
    "N": "аспарагин",
    "D": "аспарагиновая кислота",
    "C": "цистеин",
    "Q": "глутамин",
    "E": "глутаминовая кислота",
    "G": "глицин",
    "H": "гистидин",
    "I": "изолейцин",
    "L": "лейцин",
    "K": "лизин",
    "M": "метионин",
    "F": "фенилаланин",
    "P": "пролин",
    "S": "серин",
    "T": "треонин",
    "W": "триптофан",
    "Y": "тирозин",
    "V": "валин",
    "*": "стоп-кодон",
    "-": "делеция",
}

NEXTCLADE_QC_COLUMNS = [
    "accession",
    "accession_version",
    "country",
    "collection_date",
    "index",
    "seqName",
    "clade",
    "Nextclade_pango",
    "qc.overallScore",
    "qc.overallStatus",
    "coverage",
    "cdsCoverage",
    "totalMissing",
    "totalNonACGTNs",
    "totalSubstitutions",
    "totalDeletions",
    "totalInsertions",
    "totalFrameShifts",
    "totalAminoacidSubstitutions",
    "totalAminoacidDeletions",
    "totalAminoacidInsertions",
    "totalUnknownAa",
    "qc.missingData.status",
    "qc.missingData.score",
    "qc.mixedSites.status",
    "qc.mixedSites.score",
    "qc.privateMutations.status",
    "qc.privateMutations.score",
    "qc.snpClusters.status",
    "qc.snpClusters.score",
    "qc.frameShifts.status",
    "qc.frameShifts.score",
    "qc.stopCodons.status",
    "qc.stopCodons.score",
    "warnings",
    "errors",
]

NEXTCLADE_MUTATION_COLUMNS = [
    ("substitutions", "nucleotide_substitution"),
    ("deletions", "nucleotide_deletion"),
    ("insertions", "nucleotide_insertion"),
    ("frameShifts", "frameshift"),
    ("aaSubstitutions", "amino_acid_substitution"),
    ("aaDeletions", "amino_acid_deletion"),
    ("aaInsertions", "amino_acid_insertion"),
]


def existing_columns(df, columns):
    return [column for column in columns if column in df.columns]


def build_nextclade_qc(nextclade_df):
    columns = existing_columns(nextclade_df, NEXTCLADE_QC_COLUMNS)
    return nextclade_df.loc[:, columns].copy()


def split_nextclade_list(value):
    if pd.isna(value):
        return []

    text = str(value).strip()
    if text in ("", "N/A"):
        return []

    return [item.strip() for item in text.split(",") if item.strip()]


def parse_amino_acid_item(item):
    if ":" not in item:
        return "", item

    gene, mutation = item.split(":", 1)
    return gene, mutation


def extract_accession_version(seq_name):
    return str(seq_name).split()[0] if str(seq_name).strip() else ""


def build_nextclade_mutations(nextclade_df):
    id_columns = existing_columns(nextclade_df, ["index", "seqName", "clade", "Nextclade_pango"])
    rows = []

    for _, row in nextclade_df.iterrows():
        base = {column: row[column] for column in id_columns}

        for source_column, mutation_type in NEXTCLADE_MUTATION_COLUMNS:
            if source_column not in nextclade_df.columns:
                continue

            for item in split_nextclade_list(row[source_column]):
                gene = ""
                mutation = item
                if mutation_type.startswith("amino_acid"):
                    gene, mutation = parse_amino_acid_item(item)

                rows.append({
                    **base,
                    "source_column": source_column,
                    "mutation_type": mutation_type,
                    "gene": gene,
                    "mutation": mutation,
                    "raw_value": item,
                })

    return pd.DataFrame(
        rows,
        columns=[
            *id_columns,
            "source_column",
            "mutation_type",
            "gene",
            "mutation",
            "raw_value",
        ],
    )


def distribution_rows(df, column, section):
    if column not in df.columns:
        return []

    total = len(df)
    counts = df[column].replace("", "empty").value_counts(dropna=False)
    return [
        {
            "section": section,
            "value": value,
            "count": int(count),
            "percent": round(count / total * 100, 2) if total else 0,
        }
        for value, count in counts.items()
    ]


def build_nextclade_summary(nextclade_df):
    rows = [{
        "section": "samples",
        "value": "total",
        "count": len(nextclade_df),
        "percent": 100 if len(nextclade_df) else 0,
    }]

    rows.extend(distribution_rows(nextclade_df, "qc.overallStatus", "qc.overallStatus"))
    rows.extend(distribution_rows(nextclade_df, "clade", "clade"))
    rows.extend(distribution_rows(nextclade_df, "Nextclade_pango", "Nextclade_pango"))

    return pd.DataFrame(rows, columns=["section", "value", "count", "percent"])


def metadata_frame(metadata_rows):
    columns = ["accession", "accession_version", "country", "collection_date"]
    metadata_df = pd.DataFrame(metadata_rows)
    for column in columns:
        if column not in metadata_df.columns:
            metadata_df[column] = ""
    return metadata_df.loc[:, columns].copy()


def enrich_nextclade_with_metadata(nextclade_df, metadata_rows):
    enriched = nextclade_df.copy()
    enriched["accession_version"] = enriched["seqName"].map(extract_accession_version) if "seqName" in enriched else ""
    return enriched.merge(metadata_frame(metadata_rows), on="accession_version", how="left")


def enrich_mutations_with_metadata(mutations_df, metadata_rows):
    mutations = mutations_df.copy()
    if "accession_version" not in mutations.columns:
        mutations["accession_version"] = mutations["seqName"].map(extract_accession_version)
    enriched = mutations.merge(metadata_frame(metadata_rows), on="accession_version", how="left")
    metadata_columns = ["accession", "accession_version", "country", "collection_date"]
    other_columns = [column for column in enriched.columns if column not in metadata_columns]
    return enriched.loc[:, metadata_columns + other_columns]


def sample_percent(count, total):
    return round(count / total * 100, 2) if total else 0


def unique_sample_count(series):
    return series.replace("", pd.NA).dropna().nunique()


def build_nextclade_gene_summary(mutations_df):
    amino_acid_mutations = mutations_df[
        mutations_df["mutation_type"].isin([
            "amino_acid_substitution",
            "amino_acid_deletion",
            "amino_acid_insertion",
        ])
    ].copy()

    if amino_acid_mutations.empty:
        return pd.DataFrame(columns=[
            "gene",
            "mutation_type",
            "mutation_observations",
            "sample_count",
        ])

    summary = (
        amino_acid_mutations
        .groupby(["gene", "mutation_type"], dropna=False)
        .agg(
            mutation_observations=("mutation", "size"),
            sample_count=("accession_version", unique_sample_count),
        )
        .reset_index()
        .sort_values(["mutation_observations", "sample_count", "gene"], ascending=[False, False, True])
    )
    return summary


def build_nextclade_top_mutations(mutations_df, total_samples):
    if mutations_df.empty:
        return pd.DataFrame(columns=[
            "mutation_type",
            "gene",
            "mutation",
            "raw_value",
            "sample_count",
            "percent_of_samples",
            "mutation_observations",
        ])

    summary = (
        mutations_df
        .groupby(["mutation_type", "gene", "mutation", "raw_value"], dropna=False)
        .agg(
            sample_count=("accession_version", unique_sample_count),
            mutation_observations=("mutation", "size"),
        )
        .reset_index()
    )
    summary["percent_of_samples"] = summary["sample_count"].map(lambda count: sample_percent(count, total_samples))
    return summary.sort_values(
        ["sample_count", "mutation_observations", "mutation_type", "gene", "mutation"],
        ascending=[False, False, True, True, True],
    )


def build_country_summary(enriched_nextclade_df):
    rows = []
    if enriched_nextclade_df.empty or "country" not in enriched_nextclade_df.columns:
        return pd.DataFrame(columns=["country", "section", "value", "count", "percent_within_country"])

    country_counts = enriched_nextclade_df["country"].replace("", "unknown").value_counts(dropna=False)
    for country, total in country_counts.items():
        rows.append({
            "country": country,
            "section": "samples",
            "value": "total",
            "count": int(total),
            "percent_within_country": 100,
        })

    for column, section in [
        ("qc.overallStatus", "qc.overallStatus"),
        ("clade", "clade"),
        ("Nextclade_pango", "Nextclade_pango"),
    ]:
        if column not in enriched_nextclade_df.columns:
            continue
        grouped = (
            enriched_nextclade_df
            .assign(country=enriched_nextclade_df["country"].replace("", "unknown"))
            .assign(value=enriched_nextclade_df[column].replace("", "empty"))
            .groupby(["country", "value"], dropna=False)
            .size()
            .reset_index(name="count")
        )
        for _, row in grouped.iterrows():
            country_total = int(country_counts[row["country"]])
            rows.append({
                "country": row["country"],
                "section": section,
                "value": row["value"],
                "count": int(row["count"]),
                "percent_within_country": sample_percent(int(row["count"]), country_total),
            })

    return pd.DataFrame(rows, columns=["country", "section", "value", "count", "percent_within_country"])


def build_country_mutations(enriched_mutations_df, enriched_nextclade_df):
    if enriched_mutations_df.empty:
        return pd.DataFrame(columns=[
            "country",
            "mutation_type",
            "gene",
            "mutation",
            "raw_value",
            "sample_count",
            "country_sample_count",
            "percent_of_country_samples",
            "mutation_observations",
        ])

    country_sample_counts = (
        enriched_nextclade_df
        .assign(country=enriched_nextclade_df["country"].replace("", "unknown"))
        .groupby("country")["accession_version"]
        .nunique()
    )

    grouped = (
        enriched_mutations_df
        .assign(country=enriched_mutations_df["country"].replace("", "unknown"))
        .groupby(["country", "mutation_type", "gene", "mutation", "raw_value"], dropna=False)
        .agg(
            sample_count=("accession_version", unique_sample_count),
            mutation_observations=("mutation", "size"),
        )
        .reset_index()
    )
    grouped["country_sample_count"] = grouped["country"].map(country_sample_counts).fillna(0).astype(int)
    grouped["percent_of_country_samples"] = grouped.apply(
        lambda row: sample_percent(row["sample_count"], row["country_sample_count"]),
        axis=1,
    )
    return grouped.sort_values(
        ["country", "sample_count", "mutation_observations", "mutation_type", "gene", "mutation"],
        ascending=[True, False, False, True, True, True],
    )


def parse_amino_acid_substitution(mutation):
    text = str(mutation)
    if len(text) < 3:
        return "", "", ""

    aa_from = text[0]
    aa_to = text[-1]
    position = text[1:-1]
    if not position.isdigit():
        return "", "", ""

    return aa_from, position, aa_to


def build_amino_acid_changes(mutations_df):
    substitutions = mutations_df[mutations_df["mutation_type"] == "amino_acid_substitution"].copy()
    rows = []

    for _, row in substitutions.iterrows():
        aa_from, position, aa_to = parse_amino_acid_substitution(row["mutation"])
        if not aa_from or not aa_to:
            continue

        for role, amino_acid in [("reference_amino_acid", aa_from), ("sample_amino_acid", aa_to)]:
            rows.append({
                "role": role,
                "amino_acid": amino_acid,
                "amino_acid_name": AMINO_ACID_NAMES.get(amino_acid, ""),
                "amino_acid_name_ru": AMINO_ACID_NAMES_RU.get(amino_acid, ""),
                "gene": row["gene"],
                "position": position,
                "mutation": row["mutation"],
                "accession_version": row["accession_version"],
            })

    changes = pd.DataFrame(rows)
    if changes.empty:
        return pd.DataFrame(columns=[
            "role",
            "amino_acid",
            "amino_acid_name",
            "amino_acid_name_ru",
            "mutation_observations",
            "sample_count",
        ])

    return (
        changes
        .groupby(["role", "amino_acid", "amino_acid_name", "amino_acid_name_ru"], dropna=False)
        .agg(
            mutation_observations=("mutation", "size"),
            sample_count=("accession_version", unique_sample_count),
        )
        .reset_index()
        .sort_values(["role", "mutation_observations", "sample_count", "amino_acid"], ascending=[True, False, False, True])
    )


def build_amino_acid_changes_by_gene(mutations_df):
    substitutions = mutations_df[mutations_df["mutation_type"] == "amino_acid_substitution"].copy()
    rows = []

    for _, row in substitutions.iterrows():
        aa_from, position, aa_to = parse_amino_acid_substitution(row["mutation"])
        if not aa_from or not aa_to:
            continue

        for role, amino_acid in [("reference_amino_acid", aa_from), ("sample_amino_acid", aa_to)]:
            rows.append({
                "gene": row["gene"],
                "role": role,
                "amino_acid": amino_acid,
                "amino_acid_name": AMINO_ACID_NAMES.get(amino_acid, ""),
                "amino_acid_name_ru": AMINO_ACID_NAMES_RU.get(amino_acid, ""),
                "mutation": row["mutation"],
                "accession_version": row["accession_version"],
            })

    changes = pd.DataFrame(rows)
    if changes.empty:
        return pd.DataFrame(columns=[
            "gene",
            "role",
            "amino_acid",
            "amino_acid_name",
            "amino_acid_name_ru",
            "mutation_observations",
            "sample_count",
        ])

    return (
        changes
        .groupby(["gene", "role", "amino_acid", "amino_acid_name", "amino_acid_name_ru"], dropna=False)
        .agg(
            mutation_observations=("mutation", "size"),
            sample_count=("accession_version", unique_sample_count),
        )
        .reset_index()
        .sort_values(
            ["gene", "role", "mutation_observations", "sample_count", "amino_acid"],
            ascending=[True, True, False, False, True],
        )
    )


def build_run_metadata(
    nextclade_df,
    metadata_rows,
    cds_rows,
    input_gb,
    input_fasta,
    output_xlsx,
    nextclade_output_tsv,
    nextclade_aligned_fasta,
    nextclade_exe,
    nextclade_dataset,
    nextclade_version,
    run_timestamp,
    nextclade_dataset_info=None,
):
    nextclade_dataset_info = nextclade_dataset_info or {}
    rows = [
        ("run_timestamp", run_timestamp),
        ("nextclade_version", nextclade_version),
        ("nextclade_exe", str(nextclade_exe)),
        ("nextclade_dataset", str(nextclade_dataset)),
        ("nextclade_dataset_name_requested", os.environ.get("NEXTCLADE_DATASET_NAME", "")),
        ("nextclade_dataset_tag_requested", os.environ.get("NEXTCLADE_DATASET_TAG", "")),
        ("container_image", os.environ.get("SARS2_PIPELINE_CONTAINER_IMAGE", "")),
        ("input_genbank", str(input_gb)),
        ("input_fasta", str(input_fasta)),
        ("output_excel", str(output_xlsx)),
        ("nextclade_tsv", str(nextclade_output_tsv)),
        ("nextclade_aligned_fasta", str(nextclade_aligned_fasta)),
        ("genbank_record_count", len(metadata_rows)),
        ("cds_count", len(cds_rows)),
        ("nextclade_result_rows", len(nextclade_df)),
        ("nextclade_result_columns", len(nextclade_df.columns)),
        ("mutation_source", "Nextclade_results"),
        ("self_mutation_analysis", "disabled"),
        ("pairwise_alignment_in_project_code", "disabled"),
    ]
    rows.extend(
        (key, value)
        for key, value in nextclade_dataset_info.items()
    )
    return pd.DataFrame(rows, columns=["metric", "value"])


def build_nextclade_summary_tables(nextclade_df):
    return (
        build_nextclade_qc(nextclade_df),
        build_nextclade_mutations(nextclade_df),
        build_nextclade_summary(nextclade_df),
    )
