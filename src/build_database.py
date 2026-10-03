from Bio import Entrez, SeqIO
from pathlib import Path
import pandas as pd
import re
import time
from sequence_qc import check_sequence
from collections import Counter

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
        current_batch_size = min(
            batch_size,
            count - start
        )

        with Entrez.esummary(
            db = "nucleotide",
            retstart = start,
            retmax = current_batch_size,
            retmode = "xml",
            webenv = webenv,
            query_key = query_key
        ) as handle:
            results = Entrez.read(handle)

        batch_accessions = [
            record["AccessionVersion"]
            for record in results
        ]

        if len(batch_accessions) != current_batch_size:
            raise RuntimeError(
                f"Expected {current_batch_size} accessions"
                f"at position {start}, "
                f"but received {len(batch_accessions)}."
            )

        accessions.extend(batch_accessions)

        print(f"Retrieved {len(accessions)} of {count} accessions.")

        time.sleep(0.4)

    return accessions

def save_manifest(accessions):
    """Save and validate the accession manifest."""

    manifest_file = raw_dir / "accession_manifest.txt"

    if len(accessions) != len(set(accessions)):
        raise ValueError(
            "Duplicate accessions found in the manifest."
        )

    with open(manifest_file, "w") as handle:
        
        for accession in accessions:
            handle.write(f"{accession}\n")

    print(f"\nSaved {len(accessions)} accessions.")
    print(f"Manifest: {manifest_file}")
    
    return manifest_file

def load_manifest():
    """Load the saved accession manifest."""

    manifest_file = raw_dir / "accession_manifest.txt"

    with open(manifest_file, "r") as handle:
        accessions = [
            line.strip()
            for line in handle
            if line.strip()
        ]

    return accessions

def compare_existing_records(accessions, records):
    """Compare downloaded records agains the manifest"""

    manifest_set = set(accessions)

    downloaded_set = {
        record.id
        for record in records
    }

    existing = manifest_set & downloaded_set
    missing = manifest_set - downloaded_set
    outdated = downloaded_set - manifest_set

    print("\nAccession comparison:")
    print(f"Manifest accessions: {len(manifest_set)}")
    print(f"Downloaded accessions: {len(downloaded_set)}")
    print(f"Already available: {len(existing)}")
    print(f"Still missing: {len(missing)}")
    print(f"Not in manifest: {len(outdated)}")

    return existing, missing, outdated

def fetch_records_by_accession(accessions):
    """Retrieve GenBank records using accession numbers."""

    with Entrez.efetch(
        db = "nucleotide",
        id = ",".join(accessions),
        rettype = "gb",
        retmode = "text"
    ) as handle:

        records = list(SeqIO.parse(handle, "genbank"))

    return records

def fetch_with_retries(accessions, max_attempts = 3):
    """Retrieve records, retrying after temporary failures."""

    for attempt in range(1, max_attempts + 1):
        try:
            records = fetch_records_by_accession(accessions)

            requested = set(accessions)
            returned = [record.id for record in records]

            if (
                len(returned) != len(requested)
                or set(returned) != requested
            ):
                raise RuntimeError(
                    "Returned accessions do not match the request."
                )

            return records

        except Exception as error:
            print(
                f"Attempt {attempt}/{max_attempts} "
                f"failed: {error}"
            )

            if attempt == max_attempts:
                raise

            wait_time = 2 ** attempt

            print(f"Retrying in {wait_time} seconds...")
            time.sleep(wait_time)

def download_missing_records(missing, max_batches = None):
    """Download missing GenBank records in validated batches."""

    missing_accessions = sorted(missing)

    total = len(missing_accessions)

    existing_batches = list(
        raw_dir.glob("accession_batch_*.gb")
    )

    next_batch_number = max(
        (
            int(path.stem.split("_")[-1])
            for path in existing_batches
        ),
        default = -1
    ) + 1

    print(f"Next batch number: {next_batch_number}")

    print(f"\nDownloading {total} missing records...")

    for start in range(0, total, batch_size):
        batch_index = start // batch_size

        if max_batches is not None and batch_index >= max_batches:
            break

        batch_number = next_batch_number + batch_index

        batch_accessions = missing_accessions[
            start:start + batch_size
        ]

        print(
            f"\nBatch {batch_number + 1}: "
            f"requesting {len(batch_accessions)} records..."
        )

        records = fetch_with_retries(
            batch_accessions
        )

        requested = set(batch_accessions)

        returned = [record.id for record in records]

        if (
            len(returned) != len(requested)
            or set(returned) != requested
        ): 
            raise RuntimeError(
                f"Batch {batch_number + 1} validation failed."
            )

        batch_file = raw_dir / (
            f"accession_batch_{batch_number:05d}.gb"
        )

        if batch_file.exists():
            raise FileExistsError(
                f"Refusing to overwrite {batch_file}"
            )

        temporary_file = batch_file.with_suffix(".gb.tmp")

        with open(temporary_file, "w") as handle:
            SeqIO.write(records, handle, "genbank")

        temporary_file.replace(batch_file)

        print(
            f"Saved {len(records)} records "
            f"to {batch_file}"  
        )

        time.sleep(0.4)

    print("\nDownload run completed.")

def test_accession_boundary(webenv, query_key, count):
    """Test accession retrieval near the 10,000-record boundary"""

    start = 9900

    current_batch_size = min(200, count - start)

    print(
        f"\nTesting retrieval at position {start} "
        f"for {current_batch_size} records..."
    )

    with Entrez.esummary(
        db = "nucleotide",
        webenv = webenv,
        query_key = query_key,
        retstart = start,
        retmax = current_batch_size,
        retmode = "xml"
    ) as handle:
        results = Entrez.read(handle)

    accessions = [
        record["AccessionVersion"]
        for record in results
    ]

    print(f"Retrieved {len(accessions)} accessions.")

    if accessions:
        print(f"First accession: {accessions[0]}")
        print(f"Last accession: {accessions[-1]}")

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

    batch_files = sorted(
        list(raw_dir.glob("batch_*.gb"))
        + list(raw_dir.glob("accession_batch_*.gb"))
        )

    for batch_file in batch_files:
        print(f"Loading {batch_file}...")

        with open(batch_file, "r") as handle:
            batch_records = list(SeqIO.parse(handle, "genbank"))

        records.extend(batch_records)

    print(f"\nLoaded {len(records)} records.")

    return records

def validate_database(accessions, records):
    """Validate downloaded records against the manifest."""

    manifest_set = set(accessions)

    record_ids = [record.id for record in records]
    downloaded_set = set(record_ids)

    counts = Counter(record_ids)

    duplicates = {
        accession: count
        for accession, count in counts.items()
        if count > 1
    }

    missing = manifest_set - downloaded_set
    unexpected = downloaded_set - manifest_set

    print("\nDatabase integrity report:")
    print(f"Manifest entries {len(accessions)}")
    print(f"Loaded records: {len(records)}")
    print(f"Unique records: {len(downloaded_set)}")
    print(f"Duplicate accessions: {len(duplicates)}")
    print(f"Missing accessions: {len(missing)}")
    print(f"Unexpected accessions: {len(unexpected)}")

    assert len(accessions) == len(manifest_set)
    assert len(records) == len(accessions)
    assert not duplicates
    assert not missing
    assert not unexpected

    print("\nDatabase integrity check passed.")

if __name__ == "__main__":
    # count, webenv, query_key = search_ncbi(query)

    accessions = load_manifest()

    records = load_batches()

    validate_database(
        accessions,
        records
    )

    # existing, missing, outdated = compare_existing_records(
    #     accessions,
    #     records
    # )

    # assert len(records) == len(existing)
    # assert len(existing) + len(missing) == len(accessions)
    # assert len(outdated) == 0
    
    # download_missing_records(missing)

#     download_all(
#         webenv = webenv,
#         query_key = query_key,
#         count = 400,
#     )

    # accessions = fetch_accessions(
    #     webenv = webenv,
    #     query_key = query_key,
    #     count = count
    # )

    # assert len(accessions) == count
    # manifest_file = save_manifest(accessions)

    # print("\nManifest creation completed.")

#     print("\nFirst five accessions:")
#     for accession in accessions[:5]:
#         print(accession)
#     records = load_batches()

# print("\nParsing downloaded records...")

# rows = []

# for record in records:

#     metadata = parse_record(record)

#     vp7_cds = get_vp7_cds(record)

#     if vp7_cds is not None:
#         qc_metrics = check_sequence(vp7_cds)
#         metadata.update(qc_metrics)

#     rows.append(metadata)

# df = pd.DataFrame(rows)

# print("\nQC review summary:")
# print(df["qc_review"].value_counts())

# print("\nRecords requiring QC review:")
# print(
#     df[df["qc_review"] == True][
#         [
#             "accession",
#             "genotype",
#             "sequence_length",
#             "ambiguous_count",
#             "frame_valid",
#             "status"       
#         ]
#     ].to_string(index = False)
# )

# print("\nSequence QC summary:")
# print(
#     df[
#         [
#             "sequence_length",
#             "ambiguous_count",
#             "ambiguous_fraction",
#             "frame_valid"
#         ]
#     ].describe(include = "all")
# )

# print("\nSequences containing ambiguoous nucleotides:")
# print(
#     df[df["ambiguous_count"] > 0][
#         ["accession", "genotype", "ambiguous_count"]
#     ].to_string(index = False)
# )

# print(
#     df[df["frame_valid"] == False][
#         ["accession", "genotype", "sequence_length"]
#     ].to_string(index = False)
# )

# print("\nRecord status:")
# print(df["status"].value_counts())

# print("\nGenotype distribution:")
# print(df["genotype"].value_counts())

# print("\nVP7 CDS length distribution:")
# print(df["cds_length"].value_counts().sort_index())

# print("\nRecords with unusual CDS lengths:")

# unusual_records = df[
#     (df["cds_length"].notna()) &
#     (df["cds_length"] != 981)
# ]

# print(
#     unusual_records[
#         ["accession", "genotype", "cds_length", "status"]
#     ].to_string(index = False)
# )

# raw_file = raw_dir / "rva_vp7_raw.gb"

# with open(raw_file, "w") as handle:
#     SeqIO.write(records, handle, "genbank")

# print(f"\nSaved raw records to {raw_file}")

# metadata_file = processed_dir / "rva_vp7_metadata.csv"

# accepted_df = df[df["status"] == "accepted"]

# accepted_df.to_csv(metadata_file, index = False)

# print(f"Saved {len(accepted_df)} accepted records.")

# rejected_file = processed_dir / "rva_vp7_rejected.csv"

# rejected_df = df[df["status"] != "accepted"]

# rejected_df.to_csv(rejected_file, index = False)

# print(f"Saved {len(rejected_df)} rejected records.")

# fasta_file = processed_dir / "rva_vp7_sequences.fasta"

# accepted_accessions = set(accepted_df["accession"])

# with open(fasta_file, "w") as handle:
#     for record in records:
#         if record.id not in accepted_accessions:
#             continue

#         vp7_cds = get_vp7_cds(record)

#         genotype = accepted_df.loc[
#             accepted_df["accession"] == record.id,
#             "genotype"
#         ].iloc[0]

#         handle.write(f">{record.id}|{genotype}\n")
#         handle.write(f"{str(vp7_cds)}\n")

# print(f"Saved {len(accepted_accessions)} VP7 CDS sequences.")

# fasta_records = list(SeqIO.parse(fasta_file, "fasta"))

# print("\nFASTA validation:")
# print(f"Sequences: {len(fasta_records)}")

# assert len(fasta_records) == len(accepted_df)

# assert all(
#     len(record.seq) > 0
#     for record in fasta_records
# )

# print("FASTA validation passed.")

