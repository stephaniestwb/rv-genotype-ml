import pandas as pd

def merge_metadata_qc(metadata, qc):
    """
    Merge VP7 metadata and sequence QC results by accession.

    Preserve records without VP7 CDS features, but reject 
    unexpected missing QC results.
    """

    # Check unique accessions
    if metadata["accession"].duplicated().any():
        raise ValueError("Duplicate accessions in metadata.")

    if qc["accession"].duplicated().any():
        raise ValueError("Duplicate accessions in QC report.")

    # Merge QC results with metadata
    merged = metadata.merge(
        qc,
        on = "accession",
        how = "left", 
        validate = "one_to_one",
        indicator = True
    )

    unmatched = merged["_merge"].eq("left_only")
    unexpected_missing = (unmatched & ~merged["status"].eq("vp7_cds_not_found"))

    if unexpected_missing.any():
        accessions = merged.loc[unexpected_missing, "accession"].tolist()

        raise ValueError(f"Unexpected missing QC results: {accessions}")

    return merged

def assign_preliminary_eligibility(merged):
    """
    Identify sequences meeting preliminary trianing requirements.

    Requirements:
    1. A single accepted G genotype.
    2. Passing VP7 sequence QC.
    """

    merged = merged.copy()

    merged["preliminary_eligible"] = (
        merged["status"].eq("accepted")
        & merged["qc_pass"].eq(True)
    )

    return merged

if __name__ == "__main__":
    
    metadata = pd.read_csv("data/processed/rva_vp7_all_metadata.csv")
    qc = pd.read_csv("data/processed/rva_vp7_qc_report.csv")

    print("Metadata columns")
    print(metadata.columns.tolist())

    print("\nQC columns:")
    print(qc.columns.tolist())

    merged = merge_metadata_qc(metadata, qc)

    print("\nMerge results:")
    print(merged["_merge"].value_counts())

    merged = assign_preliminary_eligibility(merged)

    print("\nPreliminary training eligibility:")
    print(merged["preliminary_eligible"].value_counts())

    print("\nEligible sequences by genotype:")
    print(
        merged.loc[
            merged["preliminary_eligible"],
            "genotype"
        ].value_counts()
    )

    output_path = "data/processed/rva_vp7_merged_qc.csv"

    merged.to_csv("data/processed/rva_vp7_merged_qc.csv", index = False)

    print("\nMerged QC report saved.")


