def build_qc_summary(records, cds_rows, reference_record, reference_accession,
                     length_mismatch_records, not_comparable_records, mutation_count):
    records_with_n = sum("N" in str(record.seq).upper() for record in records)
    reference_version = reference_record.id if reference_record is not None else ""

    return [
        {"metric": "reference_accession", "value": reference_accession},
        {"metric": "reference_accession_version", "value": reference_version},
        {"metric": "record_count", "value": len(records)},
        {"metric": "cds_count", "value": len(cds_rows)},
        {"metric": "records_with_N", "value": records_with_n},
        {"metric": "records_with_length_different_from_reference", "value": len(length_mismatch_records)},
        {"metric": "mutation_count", "value": mutation_count},
        {"metric": "records_requiring_alignment", "value": len(not_comparable_records)},
    ]
