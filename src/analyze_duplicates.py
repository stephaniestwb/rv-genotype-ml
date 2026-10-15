import hashlib
import pandas as pd
from pathlib import Path
from collections import Counter, defaultdict

from Bio import SeqIO
from vp7_parser import get_vp7_cds, parse_record, get_source_qualifiers

def sequence_hash(sequence):
    """
    Generate a reproducible SHA-256 hash for a nucleotide sequence.

    Sequences differing only in letter case receive the same hash.
    """

    normalized = str(sequence).upper()

    return hashlib.sha256(
        normalized.encode("utf-8")
    ).hexdigest()

def load_vp7_sequences(raw_dir):
    """Read all accession batch files and extract VP7 CDS sequences."""

    raw_dir = Path(raw_dir)

    sequences = []
    batch_files = sorted(
        list(raw_dir.glob("accession_batch_*.gb"))
        + list(raw_dir.glob("batch_*.gb"))
    )

    for batch_file in batch_files:
        print(f"Reading {batch_file.name}")

        for record in SeqIO.parse(batch_file, "genbank"):
            vp7_cds = get_vp7_cds(record)

            if vp7_cds is None:
                continue

            metadata = parse_record(record)

            sequences.append({
                "accession": record.id,
                "genotype": metadata["genotype"],
                "status": metadata["status"],
                "sequence": str(vp7_cds).upper(),
                "sequence_hash": sequence_hash(vp7_cds),
            })

    return sequences

def summarize_duplicates(sequences):
    """Summarize exact VP7 nucleotide sequence duplication."""

    hash_counts = Counter(item["sequence_hash"] for item in sequences)
    group_sizes = Counter(hash_counts.values())

    duplicate_groups = sum(count > 1 for count in hash_counts.values())
    records_in_duplicate_groups = sum(count for count in hash_counts.values()
                                      if count > 1)

    print("\n--- VP7 Duplicate Summary ---")
    print("Unique sequence groups:", len(hash_counts))
    print("Groups containing duplicates:", duplicate_groups)
    print("Records in duplicate groups:", records_in_duplicate_groups)
    print("Largest duplicate group:", max(hash_counts.values()))

    print("\nGroup size distribution:")
    for size, number_of_groups in sorted(group_sizes.items()):
        print(f"{size} accession(s): {number_of_groups} groups")

def find_genotype_conflicts(sequences):
    """
    Identify identical VP7 sequences assigned different G genotypes.

    Only accepted genotype annotations are compared.
    """

    groups = defaultdict(list)

    for item in sequences:
        groups[item["sequence_hash"]].append(item)

    conflicts = []

    for seq_hash, records in groups.items():
        genotypes = {
            record["genotype"]
            for record in records
            if record["status"] == "accepted"
            and record["genotype"] is not None 
        }

        if len(genotypes) > 1:
            conflicts.append({
                "sequence_hash": seq_hash,
                "genotypes": sorted(genotypes),
                "records": [
                    {
                        "accession": record["accession"],
                        "genotype": record["genotype"],
                        "status": record["status"]
                    } 
                    for record in records
                ],
                "group_size": len(records),
            })

    return conflicts

def save_sequence_identity(sequences, output_path):
    """
    Save VP7 sequence identity formation to a CSV file.

    Each row represents one GenBank accession.
    """

    records = []

    for item in sequences:
        records.append({
            "accession": item["accession"],
            "genotype": item["genotype"],
            "status": item["status"],
            "cds_length": len(item["sequence"]),
            "sequence_hash": item["sequence_hash"],
        })

    df = pd.DataFrame(records)

    output_path = Path(output_path)
    output_path.parent.mkdir(parents = True, exist_ok = True)

    df.to_csv(output_path, index = False)

    print(f"\nSequence identity table saved: {output_path}")
    print(f"Rows: {len(df)}")
    print(f"Unique accessions: {df["accession"].nunique()}")
    print(f"Unique sequence hashes: {df["sequence_hash"].nunique()}")

def save_genotype_conflicts(conflicts, output_path):
    """
    Save genotype conflicts in a long-format CSV.

    Each row represents one accession belonging yo a 
    sequence group with conflicting accepted G genotypes.
    """

    rows = []

    for conflict in conflicts:
        for record in conflict["records"]:

            rows.append({
                "sequence_hash": conflict["sequence_hash"],
                "group_size": conflict["group_size"],
                "conflicting_genotypes": ",".join(conflict["genotypes"]),
                "accession": record["accession"],
                "genotype": record["genotype"],
                "status": record["status"],
            })

    df = pd.DataFrame(rows)

    output_path = Path(output_path)
    output_path.parent.mkdir(parents = True, exist_ok = True)

    df.to_csv(output_path, index = False)

    print("\n--- Genotype Conflict Report ---")
    print(f"Saved to: {output_path}")
    print(f"Conflicting sequence groups: {df["sequence_hash"].nunique()}")
    print(f"Records in conflicting groups {len(df)}")

    return df    

def inspect_accession_annotations(raw_dir, target_accessions):
    """Inspect GenBank annotations for selected accessions."""

    target_accessions = set(target_accessions)
    
    raw_dir = Path(raw_dir)

    batch_files = sorted(
        list(raw_dir.glob("accession_batch_*.gb"))
        + list(raw_dir.glob("batch_*.gb"))
    )

    for batch_file in batch_files:
        for record in SeqIO.parse(batch_file, "genbank"):

            if record.id not in target_accessions:
                continue

            source = get_source_qualifiers(record)

            vp7_cds = get_vp7_cds(record)

            print(
                "VP7 CDS length:", 
                len(vp7_cds) if vp7_cds is not None else None
            )

            print(f"\nAccession: {record.id}")
            print(f"Description: {record.description}")

            for field in ["serotype", "genotype", "isolate", "strain", "note"]:
                print(f"{field}: {source.get(field)}")

if __name__ == "__main__":
    
    sequences = load_vp7_sequences("data/raw")

    identity_df = save_sequence_identity(sequences, "data/processed/vp7_sequence_identity.csv")

    status_counts = Counter(item["status"] for item in sequences)

    print("\n--- Genotype Annotation Status ---")

    for status, count in sorted(status_counts.items()):
        print(f"{status}: {count}")

    print("\nTotal VP7 CDS sequences extracted:", len(sequences))
    print("Unique accessions:", len(set(item["accession"] for item in sequences)))
    print("Unique sequence hashes:", len(set(item["sequence_hash"] for item in sequences)))

    summarize_duplicates(sequences)

    conflicts = find_genotype_conflicts(sequences)

    conflict_df = save_genotype_conflicts(conflicts, "data/processed/vp7_genotype_conflicts.csv")

    print("\n--- Cross-record Gentoype Conflicts ---")
    print("Conflicting sequence groups:", len(conflicts))
    print("Records in conflicting groups:", sum(group["group_size"] for group in conflicts))

    for conflict in conflicts[:5]:
        
        print(f"Conflicting genotypes: {conflict["genotypes"]}")
        print(f"Group size: {conflict["group_size"]}")

        for record in conflict["records"][:10]:
            print(
                f"{record["accession"]}: "
                f"{record["genotype"]} "
                f"{record["status"]} "
            )

    inspect_accession_annotations("data/raw", ["AB011970.1", "AB039030.1"])
        

