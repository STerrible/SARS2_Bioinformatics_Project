from pathlib import Path

from openpyxl.utils import get_column_letter
import pandas as pd

from sars2_pipeline.config import (
    CDS_SHEET,
    AMINO_ACID_CHANGES_SHEET,
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


def adjust_column_widths(worksheet, max_width=60):
    for column_cells in worksheet.columns:
        max_length = 0
        for cell in column_cells:
            if cell.value is None:
                continue
            max_length = max(max_length, len(str(cell.value)))

        column_letter = get_column_letter(column_cells[0].column)
        worksheet.column_dimensions[column_letter].width = min(max_length + 2, max_width)


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
    run_metadata_df=None,
):
    output_xlsx = Path(output_xlsx)

    metadata_df = pd.DataFrame(metadata_rows)
    cds_df = pd.DataFrame(cds_rows)
    sequences_df = pd.DataFrame(sequence_rows)
    qc_summary_df = pd.DataFrame(qc_summary_rows)

    output_xlsx.parent.mkdir(parents=True, exist_ok=True)

    with pd.ExcelWriter(output_xlsx, engine="openpyxl") as writer:
        metadata_df.to_excel(writer, sheet_name=METADATA_SHEET, index=False)
        cds_df.to_excel(writer, sheet_name=CDS_SHEET, index=False)
        sequences_df.to_excel(writer, sheet_name=SEQUENCES_SHEET, index=False)
        qc_summary_df.to_excel(writer, sheet_name=QC_SUMMARY_SHEET, index=False)
        if nextclade_df is not None:
            nextclade_df.to_excel(writer, sheet_name=NEXTCLADE_SHEET, index=False)
        if nextclade_qc_df is not None:
            nextclade_qc_df.to_excel(writer, sheet_name=NEXTCLADE_QC_SHEET, index=False)
        if nextclade_mutations_df is not None:
            nextclade_mutations_df.to_excel(writer, sheet_name=NEXTCLADE_MUTATIONS_SHEET, index=False)
        if nextclade_summary_df is not None:
            nextclade_summary_df.to_excel(writer, sheet_name=NEXTCLADE_SUMMARY_SHEET, index=False)
        if nextclade_gene_summary_df is not None:
            nextclade_gene_summary_df.to_excel(writer, sheet_name=NEXTCLADE_GENE_SUMMARY_SHEET, index=False)
        if nextclade_top_mutations_df is not None:
            nextclade_top_mutations_df.to_excel(writer, sheet_name=NEXTCLADE_TOP_MUTATIONS_SHEET, index=False)
        if country_summary_df is not None:
            country_summary_df.to_excel(writer, sheet_name=COUNTRY_SUMMARY_SHEET, index=False)
        if country_mutations_df is not None:
            country_mutations_df.to_excel(writer, sheet_name=COUNTRY_MUTATIONS_SHEET, index=False)
        if amino_acid_changes_df is not None:
            amino_acid_changes_df.to_excel(writer, sheet_name=AMINO_ACID_CHANGES_SHEET, index=False)
        if run_metadata_df is not None:
            run_metadata_df.to_excel(writer, sheet_name=RUN_METADATA_SHEET, index=False)

        for worksheet in writer.book.worksheets:
            adjust_column_widths(worksheet)
