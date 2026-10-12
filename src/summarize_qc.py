import pandas as pd

metadata = pd.read_csv("data/processed/rva_vp7_all_metadata.csv")

qc = pd.read_csv("data/processed/rva_vp7_qc_report.csv")

print("Metadata columns")
print(metadata.columns.tolist())

print("\nQC columns:")
print(qc.columns.tolist())

# Check whether accession identifiers are unique
print("\nDuplicate accessions in metadata:")
print(metadata["accession"].duplicated().sum())

print("\nDuplicate accessions in QC report:")
print(qc["accession"].duplicated().sum())

# Merge QC results with metadata
merged = metadata.merge(
    qc,
    on = "accession",
    how = "left", 
    validate = "one_to_one",
    indicator = True
    )

print("\nMerge results:")
print(merged["_merge"].value_counts())

# Identify sequences meeting both preliminary requirements
merged["preliminary_eligible"] = (
    merged["status"].eq("accepted")
    & merged["qc_pass"].eq(True)
)

print("\nPreliminary training eligibility:")
print(merged["preliminary_eligible"].value_counts())

print("\nEligible sequences by genotype:")
print(
    merged.loc[
        merged["preliminary_eligible"],
        "genotype"
    ].value_counts()
)

merged.to_csv("data/processed/rva_vp7_merged_qc.csv", index = False)

print("\nMerged QC report saved.")