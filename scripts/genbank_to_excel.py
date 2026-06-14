from argparse import ArgumentParser
from pathlib import Path

from Bio import SeqIO
import pandas as pd


DEFAULT_REFERENCE_ACCESSION = "NC_045512"


def first_qualifier(feature, key, default=""):
    values = feature.qualifiers.get(key, [])
    return values[0] if values else default


def get_source_feature(record):
    for feature in record.features:
        if feature.type == "source":
            return feature
    return None


def get_source_qualifier(record, key):
    source = get_source_feature(record)
    if source is None:
        return ""

    values = source.qualifiers.get(key)
    return values[0] if values else ""


def extract_country(record):
    country = get_source_qualifier(record, "country")
    geo_loc_name = get_source_qualifier(record, "geo_loc_name")

    if country:
        return country, geo_loc_name

    if geo_loc_name:
        return geo_loc_name.split(":")[0].strip(), geo_loc_name

    return "", ""


def count_nucleotides(seq):
    seq = str(seq).upper()
    length = len(seq)

    a = seq.count("A")
    g = seq.count("G")
    c = seq.count("C")
    t = seq.count("T")
    n = seq.count("N")
    other = length - (a + g + c + t + n)

    def pct(count):
        return round(count / length * 100, 2) if length else 0

    return {
        "length": length,
        "A_count": a,
        "G_count": g,
        "C_count": c,
        "T_count": t,
        "N_count": n,
        "other_count": other,
        "A_%": pct(a),
        "G_%": pct(g),
        "C_%": pct(c),
        "T_%": pct(t),
        "N_%": pct(n),
        "other_%": pct(other),
    }


def get_accession(record):
    accessions = record.annotations.get("accessions") or []
    return accessions[0] if accessions else record.id.split(".")[0]


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


def feature_bounds(feature):
    return int(feature.location.start) + 1, int(feature.location.end)


def feature_strand(feature):
    if feature.location.strand == 1:
        return "+"
    if feature.location.strand == -1:
        return "-"
    return ""


def cds_row_from_feature(record, feature):
    accession = get_accession(record)
    accession_version = record.id
    start, end = feature_bounds(feature)
    gene = first_qualifier(feature, "gene")

    return {
        "accession": accession,
        "accession_version": accession_version,
        "gene": gene,
        "gene_is_missing": gene == "",
        "product": first_qualifier(feature, "product"),
        "protein_id": first_qualifier(feature, "protein_id"),
        "translation_exists": "translation" in feature.qualifiers,
        "start": start,
        "end": end,
        "strand": feature_strand(feature),
        "location": str(feature.location),
    }


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


def parse_genbank_to_excel(input_gb, output_xlsx, reference_accession=DEFAULT_REFERENCE_ACCESSION):
    input_gb = Path(input_gb)
    output_xlsx = Path(output_xlsx)

    records = list(SeqIO.parse(input_gb, "genbank"))
    if not records:
        raise ValueError(f"No GenBank records found in {input_gb}")

    reference_record = find_reference_record(records, reference_accession)
    if reference_record is None:
        available = ", ".join(get_accession(record) for record in records[:10])
        raise ValueError(
            f"Reference accession {reference_accession!r} was not found in {input_gb}. "
            f"First available accessions: {available}"
        )

    metadata_rows = []
    cds_rows = []
    sequence_rows = []

    for record in records:
        accession_version = record.id
        accession = get_accession(record)
        sequence = str(record.seq).upper()

        country, geo_loc_name = extract_country(record)
        counts = count_nucleotides(sequence)

        metadata_rows.append({
            "accession": accession,
            "accession_version": accession_version,
            "description": record.description,
            "organism": record.annotations.get("organism", ""),
            "country": country,
            "geo_loc_name": geo_loc_name,
            "collection_date": get_source_qualifier(record, "collection_date"),
            "host": get_source_qualifier(record, "host"),
            "isolate": get_source_qualifier(record, "isolate"),
            **counts,
        })

        sequence_rows.append({
            "accession": accession,
            "accession_version": accession_version,
            "sequence": sequence,
        })

        for feature in record.features:
            if feature.type == "CDS":
                cds_rows.append(cds_row_from_feature(record, feature))

    mutation_rows, length_mismatch_records, not_comparable_records, mutation_count = build_mutation_rows(
        records,
        reference_record,
        reference_accession,
    )
    qc_summary_rows = build_qc_summary(
        records,
        cds_rows,
        reference_record,
        reference_accession,
        length_mismatch_records,
        not_comparable_records,
        mutation_count,
    )

    metadata_df = pd.DataFrame(metadata_rows)
    cds_df = pd.DataFrame(cds_rows)
    sequences_df = pd.DataFrame(sequence_rows)
    mutations_df = pd.DataFrame(mutation_rows, columns=[
        "sample_accession",
        "reference_accession",
        "position",
        "ref_base",
        "sample_base",
        "mutation",
        "gene",
        "product",
        "protein_id",
        "is_in_cds",
        "requires_alignment",
        "comparison_note",
    ])
    qc_summary_df = pd.DataFrame(qc_summary_rows)

    output_xlsx.parent.mkdir(parents=True, exist_ok=True)

    with pd.ExcelWriter(output_xlsx, engine="openpyxl") as writer:
        metadata_df.to_excel(writer, sheet_name="Metadata_counts", index=False)
        cds_df.to_excel(writer, sheet_name="CDS_features", index=False)
        sequences_df.to_excel(writer, sheet_name="Sequences", index=False)
        mutations_df.to_excel(writer, sheet_name="Mutations", index=False)
        qc_summary_df.to_excel(writer, sheet_name="QC_summary", index=False)

    print(f"Done: {output_xlsx}")
    print(f"Records processed: {len(records)}")
    print(f"CDS processed: {len(cds_rows)}")
    print(f"Mutations found: {mutation_count}")
    print(f"Records requiring alignment: {len(not_comparable_records)}")


def parse_args():
    project_dir = Path(__file__).resolve().parents[1]
    parser = ArgumentParser(description="Convert GenBank records to an Excel workbook.")
    parser.add_argument(
        "--input",
        type=Path,
        default=project_dir / "data" / "raw" / "sequence.gb",
        help="Input GenBank file.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=project_dir / "results" / "genbank_table.xlsx",
        help="Output Excel workbook.",
    )
    parser.add_argument(
        "--reference-accession",
        default=DEFAULT_REFERENCE_ACCESSION,
        help="Reference accession for nucleotide mutation analysis.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    try:
        parse_genbank_to_excel(args.input, args.output, args.reference_accession)
    except ValueError as exc:
        raise SystemExit(f"Error: {exc}") from None


if __name__ == "__main__":
    main()
