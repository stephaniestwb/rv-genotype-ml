import pandas as pd

from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from vp7_parser import get_vp7_cds
from analyze_duplicates import sequence_hash

def load_ml_metadata(metadata_path):
    """Load metadata for the deduplicated VP7 ML dataset."""

    metadata_path = Path(metadata_path)

    metadata_df = pd.read_csv(metadata_path)

    required_columns = {
        "accession",
        "genotype",
        "sequence_hash",
    }

    missing_columns = required_columns - set(metadata_df.columns)

    if missing_columns:
        raise ValueError(f"Missing required metadata columns {sorted(missing_columns)}")

    return metadata_df

def find_selected_records(raw_dir, accessions):
    """Find GenBank records matching the selected accessions."""

    raw_dir = Path(raw_dir)
    selected_accessions = set(accessions)
    matched_records = {}

    gb_files = sorted(raw_dir.glob("batch_*.gb"))
    gb_files += sorted(raw_dir.glob("accession_batch_*.gb"))

    for gb_file in gb_files:
        for record in SeqIO.parse(gb_file, "genbank"):
            
            if record.id not in selected_accessions:
                continue

            if record.id in matched_records:
                print(f"\nDuplicate accession: {record.id}") #
                print(f"Found again in: {gb_file}") #

                previous_record = matched_records[record.id] #
                print(
                    "Full record sequences identical:",
                    str(previous_record.seq).upper() == str(record.seq).upper()
                ) #

                raise ValueError(f"Duplicate accession found: {record.id}")

            matched_records[record.id] = record

    return matched_records

def extract_verified_sequences(metadata_df, matched_records):
    """Extract VP7 CDS sequences and verify their stored hashes."""

    verified_sequences = {}

    for row in metadata_df.itertuples(index = False):
        accession = row.accession

        if accession not in matched_records:
            raise ValueError(f"Missing GenBank record: {accession}")

        record = matched_records[accession]
        cds = get_vp7_cds(record)

        if cds is None:
            raise ValueError(f"VP7 CDS not found: {accession}")

        actual_hash = sequence_hash(str(cds))

        if actual_hash != row.sequence_hash:
            raise ValueError(f"Sequence hash mismatch: {accession}")

        verified_sequences[accession] = cds

    return verified_sequences

def save_ml_fasta(metadata_df, verified_sequences, output_path):
    """Save verified VP7 sequences to a FASTA file."""

    output_path = Path(output_path)
    output_path.parent.mkdir(parents = True, exist_ok = True)

    fasta_records = []

    for row in metadata_df.itertuples(index = False):
        accession = row.accession
        genotype = row.genotype

        if accession not in verified_sequences:
            raise ValueError(f"Missing verified sequence: {accession}")

        fasta_record = SeqRecord(
            verified_sequences[accession],
            id = accession,
            description = f"genotype={genotype}"
        )

        fasta_records.append(fasta_record)

    SeqIO.write(fasta_records, output_path, "fasta")

    return output_path

if __name__ == "__main__":

    metadata_path = Path("data/processed/rva_vp7_ml_unique_metadata.csv")
    raw_dir = Path("data/raw")

    metadata_df = load_ml_metadata(metadata_path)
    print(f"Selected accessions: {len(metadata_df)}")

    matched_records = find_selected_records(raw_dir, metadata_df["accession"])
    print(f"Matched GenBank records: {len(matched_records)}")

    missing_accessions = (set(metadata_df["accession"]) - set(matched_records))
    print(f"Missing accessions: {len(missing_accessions)}")

    verified_sequences = extract_verified_sequences(metadata_df, matched_records)
    print(f"Verified VP7 sequences: {len(verified_sequences)}")

    output_path = save_ml_fasta(metadata_df, verified_sequences, "data/processed/rva_vp7_ml_unique.fasta")
    print(f"FASTA saved: {output_path}")