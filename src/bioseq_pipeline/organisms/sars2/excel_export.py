from pathlib import Path

import pandas as pd

from bioseq_pipeline.core.excel import write_dataframes_to_excel
from bioseq_pipeline.organisms.sars2.config import (
    ALIGNED_FASTA_SHEET,
    CDS_SHEET,
    AMINO_ACID_CHANGES_SHEET,
    AMINO_ACID_CHANGES_BY_GENE_SHEET,
    COUNTRY_MUTATIONS_SHEET,
    COUNTRY_SUMMARY_SHEET,
    METADATA_SHEET,
    NEXTCLADE_GENE_SUMMARY_SHEET,
    NEXTCLADE_MUTATIONS_SHEET,
    NEXTCLADE_QC_SHEET,
    NEXTCLADE_SHEET,
    NEXTCLADE_SUMMARY_SHEET,
    NEXTCLADE_TOP_MUTATIONS_SHEET,
    QC_SUMMARY_SHEET,
    RUN_METADATA_SHEET,
    SEQUENCES_SHEET,
)


def write_excel(
    output_xlsx,
    metadata_rows,
    cds_rows,
    sequence_rows,
    qc_summary_rows,
    nextclade_df=None,
    nextclade_qc_df=None,
    nextclade_mutations_df=None,
    nextclade_summary_df=None,
    nextclade_gene_summary_df=None,
    nextclade_top_mutations_df=None,
    country_summary_df=None,
    country_mutations_df=None,
    amino_acid_changes_df=None,
    amino_acid_changes_by_gene_df=None,
    aligned_fasta_df=None,
    run_metadata_df=None,
):
    output_xlsx = Path(output_xlsx)

    metadata_df = pd.DataFrame(metadata_rows)
    cds_df = pd.DataFrame(cds_rows)
    sequences_df = pd.DataFrame(sequence_rows)
    qc_summary_df = pd.DataFrame(qc_summary_rows)

    sheets = [
        (METADATA_SHEET, metadata_df),
        (CDS_SHEET, cds_df),
        (SEQUENCES_SHEET, sequences_df),
        (QC_SUMMARY_SHEET, qc_summary_df),
        (NEXTCLADE_SHEET, nextclade_df),
        (NEXTCLADE_QC_SHEET, nextclade_qc_df),
        (NEXTCLADE_MUTATIONS_SHEET, nextclade_mutations_df),
        (NEXTCLADE_SUMMARY_SHEET, nextclade_summary_df),
        (NEXTCLADE_GENE_SUMMARY_SHEET, nextclade_gene_summary_df),
        (NEXTCLADE_TOP_MUTATIONS_SHEET, nextclade_top_mutations_df),
        (COUNTRY_SUMMARY_SHEET, country_summary_df),
        (COUNTRY_MUTATIONS_SHEET, country_mutations_df),
        (AMINO_ACID_CHANGES_SHEET, amino_acid_changes_df),
        (AMINO_ACID_CHANGES_BY_GENE_SHEET, amino_acid_changes_by_gene_df),
        (ALIGNED_FASTA_SHEET, aligned_fasta_df),
        (RUN_METADATA_SHEET, run_metadata_df),
    ]
    write_dataframes_to_excel(output_xlsx, sheets)
