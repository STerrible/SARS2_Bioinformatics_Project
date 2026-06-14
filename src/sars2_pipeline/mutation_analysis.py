from sars2_pipeline.genbank_parser import first_qualifier, get_accession


def accession_aliases(record):
    accession = get_accession(record)
    aliases = {record.id, record.name, accession}
    aliases.update(record.annotations.get("accessions") or [])
    return {alias for value in aliases for alias in (value, value.split(".")[0]) if alias}


def find_reference_record(records, reference_accession):
    wanted = {reference_accession, reference_accession.split(".")[0]}
    for record in records:
        if accession_aliases(record) & wanted:
            return record
    return None


def require_reference_record(records, reference_accession, input_gb):
    reference_record = find_reference_record(records, reference_accession)
    if reference_record is None:
        available = ", ".join(get_accession(record) for record in records[:10])
        raise ValueError(
            f"Reference accession {reference_accession!r} was not found in {input_gb}. "
            f"First available accessions: {available}"
        )
    return reference_record


def build_reference_cds_index(reference_record):
    cds_index = []
    for feature in reference_record.features:
        if feature.type != "CDS":
            continue
        cds_index.append({
            "location": feature.location,
            "gene": first_qualifier(feature, "gene"),
            "product": first_qualifier(feature, "product"),
            "protein_id": first_qualifier(feature, "protein_id"),
        })
    return cds_index


def annotate_reference_position(position, reference_cds_index):
    zero_based_position = position - 1
    hits = [cds for cds in reference_cds_index if zero_based_position in cds["location"]]

    if not hits:
        return {
            "gene": "",
            "product": "",
            "protein_id": "",
            "is_in_cds": False,
        }

    return {
        "gene": ";".join(dict.fromkeys(hit["gene"] for hit in hits if hit["gene"])),
        "product": ";".join(dict.fromkeys(hit["product"] for hit in hits if hit["product"])),
        "protein_id": ";".join(dict.fromkeys(hit["protein_id"] for hit in hits if hit["protein_id"])),
        "is_in_cds": True,
    }


def build_mutation_rows(records, reference_record, reference_accession):
    reference_sequence = str(reference_record.seq).upper()
    reference_length = len(reference_sequence)
    reference_cds_index = build_reference_cds_index(reference_record)

    mutation_rows = []
    length_mismatch_records = []
    not_comparable_records = []

    for record in records:
        sample_accession = get_accession(record)
        if record is reference_record:
            continue

        sample_sequence = str(record.seq).upper()
        if len(sample_sequence) != reference_length:
            length_mismatch_records.append(sample_accession)
            not_comparable_records.append(sample_accession)
            mutation_rows.append({
                "sample_accession": sample_accession,
                "reference_accession": reference_accession,
                "position": "",
                "ref_base": "",
                "sample_base": "",
                "mutation": "requires_alignment",
                "gene": "",
                "product": "",
                "protein_id": "",
                "is_in_cds": False,
                "requires_alignment": True,
                "comparison_note": (
                    f"sample length {len(sample_sequence)} differs from "
                    f"reference length {reference_length}"
                ),
            })
            continue

        for index, (ref_base, sample_base) in enumerate(zip(reference_sequence, sample_sequence), start=1):
            if ref_base == sample_base:
                continue

            annotation = annotate_reference_position(index, reference_cds_index)
            mutation_rows.append({
                "sample_accession": sample_accession,
                "reference_accession": reference_accession,
                "position": index,
                "ref_base": ref_base,
                "sample_base": sample_base,
                "mutation": f"{ref_base}{index}{sample_base}",
                **annotation,
                "requires_alignment": False,
                "comparison_note": "",
            })

    mutation_count = sum(not row["requires_alignment"] for row in mutation_rows)
    return mutation_rows, length_mismatch_records, not_comparable_records, mutation_count
