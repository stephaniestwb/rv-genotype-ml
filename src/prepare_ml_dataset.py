from pathlib import Path
import pandas as pd

def inspect_input_datasets(processed_dir):
    """Inspect the datasets required for ML preparation."""

    processed_dir = Path(processed_dir)

    identity_path = processed_dir / "vp7_sequence_identity.csv"
    qc_path = processed_dir / "rva_vp7_merged_qc.csv"

    identity_df = pd.read_csv(identity_path)
    qc_df = pd.read_csv(qc_path)

    print("\n--- Sequence Identity Dataset ---")
    print(f"Rows: {len(identity_df)}")
    print(f"Unique accessions: {identity_df["accession"].nunique()}")
    print(f"Colums: {identity_df.columns.tolist()}")

    print("\n--- QC Dataset ---")
    print(f"Rows: {len(qc_df)}")
    print(f"Unique accessions: {qc_df["accession"].nunique()}")
    print(f"Columns: {qc_df.columns.tolist()}")

    return identity_df, qc_df

def merge_identity_qc(identity_df, qc_df):
    """Merge VP7 sequence identity information"""
    
    identity_subset = identity_df[["accession", "sequence_hash"]]

    merged_df = qc_df.merge(
        identity_subset,
        on = "accession",
        how = "left",
        validate = "one_to_one",
        indicator = "identity_merge"
    )

    print("\n--- Identity + QC Merge ---")
    print(f"Rows: {len(merged_df)}")
    print(f"Unique accessions: {merged_df["accession"].nunique()}")

    print("\nMerge results:")
    print(merged_df["identity_merge"].value_counts())

    return merged_df

def inspect_unmmatched_records(merged_df):
    """Inspect records without a matching VP7 sequence hash."""

    unmatched_df = merged_df[merged_df["identity_merge"] == "left_only"]

    columns = [
        "accession",
        "genotype",
        "status",
        "cds_length_x",
        "cds_length_y",
        "sequence_qc_flags",
        "annotation_qc_flags",
        "qc_pass",
        "preliminary_eligible",
    ]

    print("\n--- Unmatched Records ---")
    print(unmatched_df[columns].to_string(index = False))

    return unmatched_df

def summarize_preliminary_eligibility(merged_df):
    """Summarize existinf QC and preliminary ML eligibility."""

    print("\n--- Preliminary ML Eligibility ---")
    
    print("\nPreliminary eligibility:")
    print(merged_df["preliminary_eligible"].value_counts(dropna = False))

    print("\nSequence QC:")
    print(merged_df["qc_pass"].value_counts(dropna = False))

    eligible_df = merged_df[merged_df["preliminary_eligible"] == True]

    print("\nEligible records with missing sequence hashes:")
    print(eligible_df["sequence_hash"].isna().sum())

    print("\nGenotype distribution among eligible records:")
    print(eligible_df["genotype"].value_counts().sort_index())

    return eligible_df

def identify_genotype_conflicts(merged_df, processed_dir):
    """Flag records belonging to unresolved genotype-conflict groups."""

    processed_dir = Path(processed_dir)

    conflict_path = processed_dir / "vp7_genotype_conflicts.csv"
    conflict_df = pd.read_csv(conflict_path)
    conflict_hashes = set(conflict_df["sequence_hash"])

    merged_df = merged_df.copy()
    merged_df["genotype_conflict"] = (merged_df["sequence_hash"].isin(conflict_hashes))

    print("\n--- Genotype Conflict Screening ---")
    print(f"Conflicting sequence groups: {len(conflict_hashes)}")
    print("\nConflict status among preliminary eligible records:")

    eligible_df = merged_df[merged_df["preliminary_eligible"] == True]

    print(eligible_df["genotype_conflict"].value_counts()) 
    
    return merged_df

def assign_ml_eligibility(merged_df):
    """Assign ML eligibility based on QC and genotype conflicts."""

    merged_df = merged_df.copy()
    merged_df["ml_eligible"] = (
        (merged_df["preliminary_eligible"] == True)
        & (merged_df["genotype_conflict"] == False)
        & (merged_df["sequence_hash"].notna())
    )

    print("\n--- ML Eligibility ---")
    print(merged_df["ml_eligible"].value_counts())

    eligible_df = merged_df[merged_df["ml_eligible"]]

    print("\nUnique sequences among eligible records:")
    print(eligible_df["sequence_hash"].nunique())

    print("\nGenotypes represented:")
    print(eligible_df["genotype"].nunique())

    return merged_df

def summarize_genotype_diversity(merged_df):
    """Compare GenBank record counts with unique VP7 sequences per genotype."""

    eligible_df = merged_df[merged_df["ml_eligible"]].copy()

    genotype_summary = (
        eligible_df
        .groupby("genotype")
        .agg(
            records = ("accession", "count"),
            unique_sequences = ("sequence_hash", "nunique"),
        )
        .reset_index()
    )

    genotype_summary["duplicate_records"] = (
        genotype_summary["records"] - genotype_summary["unique_sequences"]
    )

    genotype_summary = genotype_summary.sort_values(
        "unique_sequences",
        ascending = False
    )

    print("\n--- Genotype Diversity ---")
    print(genotype_summary.to_string(index = False))

    return genotype_summary

def validate_genotype_consistency(merged_df):
    """Verify that eligible sequence hashes have consistent genotype labels."""

    eligible_df = merged_df[merged_df["ml_eligible"]]

    genotype_counts = (
        eligible_df
        .groupby("sequence_hash")["genotype"]
        .nunique()
    )

    inconsistent_groups = genotype_counts[genotype_counts > 1]

    print("\n--- Genotype Consistency Validation ---")
    print(f"Sequence groups checked: {len(genotype_counts)}")
    print(f"Inconsistent groups: {len(inconsistent_groups)}")

    if len(inconsistent_groups) > 0:
        raise ValueError("Conflicting genotype labels remain in the ML dataset.")

    print("All eligible sequence groups have consistent genotype labels.")

def select_baseline_genotype(merged_df, min_unique_sequences = 100):
    """Select genotype with enough unique VP7 sequences for baseline ML."""

    eligible_df = merged_df[merged_df["ml_eligible"]].copy()

    unique_counts = (
        eligible_df
        .groupby("genotype")["sequence_hash"]
        .nunique()
    )

    selected_genotypes = unique_counts[
        unique_counts >= min_unique_sequences
    ].index.tolist()

    baseline_df = eligible_df[
        eligible_df["genotype"].isin(selected_genotypes)
        ].copy()

    print("\n--- Baseline ML Dataset ---")
    print(f"Minimum unique sequences: {min_unique_sequences}")
    print(f"Selected genotypes: {sorted(selected_genotypes)}")
    print(f"Genotype classes: {len(selected_genotypes)}")
    print(f"Records: {len(baseline_df)}")
    print(f"Unique sequences: {baseline_df["sequence_hash"].nunique()}")

    return baseline_df

def create_unique_sequence_dataset(baseline_df):
    """Retain one representative accession per unique VP7 sequence."""

    unique_df = (
        baseline_df
        .sort_values("accession")
        .drop_duplicates(
            subset = "sequence_hash",
            keep = "first"
        ).copy()
    )

    unique_df = unique_df.reset_index(drop = True)

    print("\n--- Unique-Sequence ML Dataset ---")
    print(f"Original records {len(baseline_df)}")
    print(f"Unique sequences: {len(unique_df)}")
    print(f"Duplicates removed: {len(baseline_df) - len(unique_df)}")

    print("\nGenotype distribution:")
    print(unique_df["genotype"].value_counts().sort_index())

    return unique_df

def save_ml_dataset(unique_df, output_path):
    """Save the deduplicated VP7 ML metadata dataset."""

    output_path = Path(output_path)
    output_path.parent.mkdir(parents = True, exist_ok = True)

    unique_df.to_csv(output_path, index = False)

    print("\n--- Saved ML Dataset ---")
    print(f"Output: {output_path}")
    print(f"Records: {len(unique_df)}")
    print(f"Unique hashes: {unique_df["sequence_hash"].nunique()}")

    return output_path


if __name__ == "__main__":

    identity_df, qc_df = inspect_input_datasets("data/processed")

    merged_df = merge_identity_qc(identity_df, qc_df)

    unmatched_df = inspect_unmmatched_records(merged_df)

    eligible_df = summarize_preliminary_eligibility(merged_df)

    merged_df = identify_genotype_conflicts(merged_df, "data/processed")

    merged_df = assign_ml_eligibility(merged_df)

    genotype_summary = summarize_genotype_diversity(merged_df)

    validate_genotype_consistency(merged_df)

    baseline_df = select_baseline_genotype(merged_df, min_unique_sequences = 100)

    unique_df = create_unique_sequence_dataset(baseline_df)

    save_ml_dataset(unique_df, "data/processed/rva_vp7_ml_unique_metadata.csv")

    