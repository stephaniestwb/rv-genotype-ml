from pathlib import Path
import pandas as pd
from Bio import SeqIO

from vp7_parser import parse_record, get_vp7_cds
from sequence_qc import check_sequence

raw_dir = Path("data/raw")
processed_dir = Path("data/processed")

processed_dir.mkdir(
    parents = True,
    exist_ok = True
)

def load_records():
    """Load all downloaded GenBank batches."""

    records = []

    batch_files = sorted(
        list(raw_dir.glob("batch_*.gb"))
        + list(raw_dir.glob("accession_batch_*.gb"))
    )

    for batch_file in batch_files:
        print(f"Loading {batch_file.name}...")

        with open(batch_file, "r") as handle:
            records.extend(
                SeqIO.parse(handle, "genbank")
            )

    print(f"\nLoaded {len(records)} records.")

    return records

def process_records(records):
    """Extract metadata and calculate VP7 sequence QC metrics."""

    rows = []

    for record in records:
        metadata = parse_record(record)
        vp7_cds = get_vp7_cds(record)

        if vp7_cds is not None:
            qc_metrics = check_sequence(vp7_cds)
            metadata.update(qc_metrics)

        rows.append(metadata)

    df = pd.DataFrame(rows)

    return df

def generate_qc_report(df):
    """Display a preliminary summary of the dataset."""

    print("\nDataset summary:")
    print(f"Total records: {len(df)}")

    print("\nRecord status:")
    print(df["status"].value_counts(dropna = False))

    print("\nGenotype distribution:")
    print(df["genotype"].value_counts(dropna = False))

    print("\nVP7 CDS length distribution:")
    print(df["cds_length"].value_counts(dropna = False).sort_index())

    print("\nSequence QC summary:")
    print(df[["sequence_length", "ambiguous_count", "ambiguous_fraction"]].describe())

    print("\nRecords requiring QC review:")
    print(df["qc_review"].value_counts(dropna = False))

def save_metadata(df):
    """Save the preliminary metadata and QC table."""

    output_file = (processed_dir / "rva_vp7_all_metadata.csv")

    df.to_csv(output_file, index = False)

    print(f"\nSaved metadata to {output_file}")

def inspect_unsual_lengths(df):
    """Display records with unusually short or long VP7 CDS."""

    unusual = df[(df["cds_length"] < 900) | (df["cds_length"] > 1000)]

    columns = [
        "accession",
        "genotype",
        "status",
        "record_length",
        "cds_length",
        "description",
    ]

    print("\nUnusual CDS annotations:")
    print(f"Records requiring inspection: {len(unusual)}")

    print(unusual[columns].to_string(index = False, max_colwidth = 70))

    return unusual

if __name__ == "__main__":

    records = load_records()
    assert len(records) == 10096

    print("\nProcessing records")

    df = process_records(records)

    generate_qc_report(df)
    unusual = inspect_unsual_lengths(df)
    save_metadata(df)

    print("\nPreliminary processing completed.")
