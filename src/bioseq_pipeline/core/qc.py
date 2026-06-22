def build_qc_summary(records, cds_rows):
    records_with_n = sum("N" in str(record.seq).upper() for record in records)

    return [
        {"metric": "record_count", "value": len(records)},
        {"metric": "cds_count", "value": len(cds_rows)},
        {"metric": "records_with_N", "value": records_with_n},
    ]
