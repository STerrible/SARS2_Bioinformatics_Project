from pathlib import Path

from openpyxl.utils import get_column_letter
import pandas as pd


def adjust_column_widths(worksheet, max_width=60):
    for column_cells in worksheet.columns:
        max_length = 0
        for cell in column_cells:
            if cell.value is None:
                continue
            max_length = max(max_length, len(str(cell.value)))

        column_letter = get_column_letter(column_cells[0].column)
        worksheet.column_dimensions[column_letter].width = min(max_length + 2, max_width)


def write_dataframes_to_excel(output_xlsx, sheets):
    output_xlsx = Path(output_xlsx)
    output_xlsx.parent.mkdir(parents=True, exist_ok=True)

    with pd.ExcelWriter(output_xlsx, engine="openpyxl") as writer:
        for sheet_name, dataframe in sheets:
            if dataframe is None:
                continue
            dataframe.to_excel(writer, sheet_name=sheet_name, index=False)

        for worksheet in writer.book.worksheets:
            adjust_column_widths(worksheet)
