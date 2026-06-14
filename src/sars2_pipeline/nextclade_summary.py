import pandas as pd


NEXTCLADE_QC_COLUMNS = [
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


def build_nextclade_summary_tables(nextclade_df):
    return (
        build_nextclade_qc(nextclade_df),
        build_nextclade_mutations(nextclade_df),
        build_nextclade_summary(nextclade_df),
    )
