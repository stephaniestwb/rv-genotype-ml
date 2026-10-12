import pandas as pd

from build_database import load_records
from sequence_qc import evaluate_vp7_cds
from vp7_parser import get_vp7_feature

def run_qc(records):
    """Evaluate the selected VP7 CDS feature for each record."""

    reports = []

    for record in records: 

        feature = get_vp7_feature(record)

        if feature is None:
            continue

        report = evaluate_vp7_cds(record, feature)
        reports.append(report)

    return pd.DataFrame(reports)

if __name__ == "__main__":

    records = load_records()

    qc_df = run_qc(records)

    print(f"CDS features evaluated: {len(qc_df)}")
    print(qc_df["qc_pass"].value_counts())

    print("\nSequence QC flags:")
    print(
        qc_df["sequence_qc_flags"]
        .explode()
        .dropna()
        .value_counts()
    )

    print("\nAnnotation QC flags:")
    print(
        qc_df["annotation_qc_flags"]
        .explode()
        .dropna()
        .value_counts()
    )

    qc_df.to_csv("data/processed/rva_vp7_qc_report.csv", index = False)
    