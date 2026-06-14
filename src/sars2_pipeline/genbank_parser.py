from pathlib import Path

from Bio import SeqIO


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


def read_genbank_records(input_gb):
    input_gb = Path(input_gb)
    records = list(SeqIO.parse(input_gb, "genbank"))
    if not records:
        raise ValueError(f"No GenBank records found in {input_gb}")
    return records


def extract_genbank_tables(records):
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

    return metadata_rows, cds_rows, sequence_rows
