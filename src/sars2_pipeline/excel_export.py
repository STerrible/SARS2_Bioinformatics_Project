from pathlib import Path

import pandas as pd

from sars2_pipeline.config import (
    CDS_SHEET,
    METADATA_SHEET,
    MUTATION_COLUMNS,
    MUTATIONS_SHEET,
    QC_SUMMARY_SHEET,
    SEQUENCES_SHEET,
)


def write_excel(output_xlsx, metadata_rows, cds_rows, sequence_rows, mutation_rows, qc_summary_rows):
    output_xlsx = Path(output_xlsx)

    metadata_df = pd.DataFrame(metadata_rows)
    cds_df = pd.DataFrame(cds_rows)
    sequences_df = pd.DataFrame(sequence_rows)
    mutations_df = pd.DataFrame(mutation_rows, columns=MUTATION_COLUMNS)
    qc_summary_df = pd.DataFrame(qc_summary_rows)

    output_xlsx.parent.mkdir(parents=True, exist_ok=True)

    with pd.ExcelWriter(output_xlsx, engine="openpyxl") as writer:
        metadata_df.to_excel(writer, sheet_name=METADATA_SHEET, index=False)
        cds_df.to_excel(writer, sheet_name=CDS_SHEET, index=False)
        sequences_df.to_excel(writer, sheet_name=SEQUENCES_SHEET, index=False)
        mutations_df.to_excel(writer, sheet_name=MUTATIONS_SHEET, index=False)
        qc_summary_df.to_excel(writer, sheet_name=QC_SUMMARY_SHEET, index=False)
