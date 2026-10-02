from Bio import Entrez, SeqIO
from pathlib import Path
import pandas as pd
import re
import time
from sequence_qc import check_sequence

from vp7_parser import parse_record, get_vp7_cds

Entrez.email = "stampwberth.1@osu.edu"

query = (
    '"Rotavirus A"[Organism] '
    'AND VP7[Title] '
    'AND "complete cds"[Title] '
) 

batch_size = 200

raw_dir = Path("data/raw")
processed_dir = Path("data/processed")

raw_dir.mkdir(parents = True, exist_ok = True)
processed_dir.mkdir(parents = True, exist_ok = True)

def search_ncbi(query):
    """Search NCBI and store the results in Web History."""

    print("Searching NCBI...")

    with Entrez.esearch(
        db = "nucleotide",
        term = query,
        usehistory = "y",
        retmax = 0
    ) as handle:
        results = Entrez.read(handle)

    count = int(results["Count"])
    webenv = results["WebEnv"]
    query_key = results["QueryKey"]

    print(f"Found {count} candidate records.")
    
    return count, webenv, query_key

def fetch_accessions(webenv, query_key, count):
    """Retrieve the accession numbers from an NCBI search."""

    print("\nRetrieving accession list...")

    accessions = []

    for start in range(0, count, batch_size):
        current_batch_size = min(batch_size, count - start)
        with Entrez.efetch(
            db = "nucleotide",
            rettype = "acc",
            retmode = "text",
            retstart = 0,
            retmax = count,
            webenv = webenv,
            query_key = query_key
        ) as handle:
            batch_accessions = handle.read().splitlines()

        if len(batch_accessions) != current_batch_size:
            raise RuntimeError(
                f"Expected {current_batch_size} accessions "
                f"at position {start}, "
                f"but received {len(batch_accessions)}."
            ) 

        accessions.extend(batch_accessions)
        print(
            f"Retrieved {len(accessions)} "
            f"of {count} accessions."
        )

    return accessions

def fetch_batch(webenv, query_key, start, batch_size):
    """Retrieve a batch of GenBank records from NCBI Web History."""

    print(f"Fetching records {start + 1} to {start + batch_size}...")

    with Entrez.efetch(
        db = "nucleotide",
        rettype = "gb",
        retmode = "text",
        retstart = start,
        retmax = batch_size,
        webenv = webenv,
        query_key = query_key
    ) as handle: 
        records = list(SeqIO.parse(handle, "genbank"))

    return records

def save_batch(records, start):
    """Save a downloaded batch as a GenBank file."""

    batch_file = raw_dir / f"batch_{start:05d}.gb"

    with open(batch_file, "w") as handle:
        SeqIO.write(records, handle, "genbank")

    return batch_file

def download_all(count, webenv, query_key):
    """Download all matching records in batches."""

    for start in range(0, count, batch_size):

        current_batch_size = min(batch_size, count - start)

        batch_file = raw_dir / f"batch_{start:05d}.gb"

        if batch_file.exists():
            print(f"Skipping existing batch: {batch_file}")
            continue

        records = fetch_batch(
            webenv = webenv,
            query_key = query_key,
            start = start,
            batch_size = current_batch_size
        )

        if len(records) != current_batch_size:
            raise RuntimeError(
                f"Expected {current_batch_size} records, "
                f"but received {len(records)}."
            )
        save_batch(records, start)

        print(f"Saved {len(records)} records to {batch_file}")

        time.sleep(0.4)

    print("\nDownload completed.")

def load_batches():
    """Load all downloaded GenBank batches."""

    records = []

    batch_files = sorted(raw_dir.glob("batch_*.gb"))

    for batch_file in batch_files:
        print(f"Loading {batch_file}...")

        with open(batch_file, "r") as handle:
            batch_records = list(SeqIO.parse(handle, "genbank"))

        records.extend(batch_records)

    print(f"\nLoaded {len(records)} records.")

    return records

if __name__ == "__main__":
    count, webenv, query_key = search_ncbi(query)

    download_all(
        webenv = webenv,
        query_key = query_key,
        count = 400,
    )

    accessions = fetch_accessions(
        webenv = webenv,
        query_key = query_key,
        count = count
    )

    assert len(accessions) == count
    assert len(accessions) == len(set(accessions))

    print("\nAccession list validation passed.")

    print("\nFirst five accessions:")
    for accession in accessions[:5]:
        print(accession)
    records = load_batches()

print("\nParsing downloaded records...")

rows = []

for record in records:

    metadata = parse_record(record)

    vp7_cds = get_vp7_cds(record)

    if vp7_cds is not None:
        qc_metrics = check_sequence(vp7_cds)
        metadata.update(qc_metrics)

    rows.append(metadata)

df = pd.DataFrame(rows)

print("\nQC review summary:")
print(df["qc_review"].value_counts())

print("\nRecords requiring QC review:")
print(
    df[df["qc_review"] == True][
        [
            "accession",
            "genotype",
            "sequence_length",
            "ambiguous_count",
            "frame_valid",
            "status"       
        ]
    ].to_string(index = False)
)

print("\nSequence QC summary:")
print(
    df[
        [
            "sequence_length",
            "ambiguous_count",
            "ambiguous_fraction",
            "frame_valid"
        ]
    ].describe(include = "all")
)

print("\nSequences containing ambiguoous nucleotides:")
print(
    df[df["ambiguous_count"] > 0][
        ["accession", "genotype", "ambiguous_count"]
    ].to_string(index = False)
)

print(
    df[df["frame_valid"] == False][
        ["accession", "genotype", "sequence_length"]
    ].to_string(index = False)
)

print("\nRecord status:")
print(df["status"].value_counts())

print("\nGenotype distribution:")
print(df["genotype"].value_counts())

print("\nVP7 CDS length distribution:")
print(df["cds_length"].value_counts().sort_index())

print("\nRecords with unusual CDS lengths:")

unusual_records = df[
    (df["cds_length"].notna()) &
    (df["cds_length"] != 981)
]

print(
    unusual_records[
        ["accession", "genotype", "cds_length", "status"]
    ].to_string(index = False)
)

raw_file = raw_dir / "rva_vp7_raw.gb"

with open(raw_file, "w") as handle:
    SeqIO.write(records, handle, "genbank")

print(f"\nSaved raw records to {raw_file}")

metadata_file = processed_dir / "rva_vp7_metadata.csv"

accepted_df = df[df["status"] == "accepted"]

accepted_df.to_csv(metadata_file, index = False)

print(f"Saved {len(accepted_df)} accepted records.")

rejected_file = processed_dir / "rva_vp7_rejected.csv"

rejected_df = df[df["status"] != "accepted"]

rejected_df.to_csv(rejected_file, index = False)

print(f"Saved {len(rejected_df)} rejected records.")

fasta_file = processed_dir / "rva_vp7_sequences.fasta"

accepted_accessions = set(accepted_df["accession"])

with open(fasta_file, "w") as handle:
    for record in records:
        if record.id not in accepted_accessions:
            continue

        vp7_cds = get_vp7_cds(record)

        genotype = accepted_df.loc[
            accepted_df["accession"] == record.id,
            "genotype"
        ].iloc[0]

        handle.write(f">{record.id}|{genotype}\n")
        handle.write(f"{str(vp7_cds)}\n")

print(f"Saved {len(accepted_accessions)} VP7 CDS sequences.")

fasta_records = list(SeqIO.parse(fasta_file, "fasta"))

print("\nFASTA validation:")
print(f"Sequences: {len(fasta_records)}")

assert len(fasta_records) == len(accepted_df)

assert all(
    len(record.seq) > 0
    for record in fasta_records
)

print("FASTA validation passed.")

