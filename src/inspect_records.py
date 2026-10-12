import pandas as pd
from Bio.Seq import Seq

from build_database import load_records

from vp7_parser import  (
    get_source_qualifiers, 
    extract_g_genotype,
    find_vp7_features
    )
from sequence_qc import (
    check_vp7_annotation,
    check_vp7_cds_integrity, 
    evaluate_vp7_cds
    )

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

def print_failed_record(record):
    """Print useful annotations for records with no confident genotype."""

    source = get_source_qualifiers(record)

    genotype, status = extract_g_genotype(record, source)

    if status != "accepted":
        return
    
    print("\n" + "=" * 60)
    print(f"Accession: {record.id}")
    print(f"Status: {status}")
    print(f"Description: {record.description}")

    print("\nSource qualifiers:")

    for key, values in source.items():
        print(f"{key}: {values}")

def inspect_record_features(record):
    """Display the annotated features of a GenBank record."""

    print("\n" + "=" * 60)
    print(f"Accession: {record.id}")
    print(f"Description: {record.description}")
    print(f"Record length: {len(record.seq)} nt")

    for feature in record.features:
        print(f"\nFeature type: {feature.type}")
        print(f"Location: {feature.location}")
        print(f"Feature length: {len(feature.extract(record.seq))} nt")

        for key in ["gene", "product", "note", "codon_start"]:
            if key in feature.qualifiers:
                print(f"{key}: {feature.qualifiers[key]}")

        if "translation" in feature.qualifiers:
            protein = feature.qualifiers["translation"][0]
            print(f"Protein length: {len(protein)} amino acids")

def inspect_orfs(records):
    """Inspect the three forward reading frames."""

    sequence = record.seq

    print(f"\nORF inspection: {record.id}")

    for frame in range(3):
        coding_sequence = sequence[frame:]

        # Translation requires complete codons
        remainder = len(coding_sequence) % 3

        if remainder:
            coding_sequence = coding_sequence[:-remainder]

        protein = coding_sequence.translate()

        print(f"\nFrame {frame + 1}")
        print(f"Translated length: {len(protein)} aa")
        print(f"Stop codons: {protein.count('*')}")
        print(f"First 50 amino acids: {protein[:50]}")

def find_longest_orf(record):
    """Find the longest stop-free region in each forward reading frame."""

    sequence = record.seq

    print(f"\nLongest stop-free regions: {record.id}")

    for frame in range(3):
        coding_sequence = sequence[frame:]
        remainder = len(coding_sequence) % 3

        if remainder:
            coding_sequence = coding_sequence[:-remainder]

        protein = str(coding_sequence.translate())
        
        regions = protein.split("*")
        longest = max(regions, key = len)

        # Position of the longest region in the translated sequence
        start_aa = protein.find(longest)

        # Convert amino acid coordinates to nucleotide coordinates
        start_nt = frame + start_aa * 3
        end_nt = start_nt + len(longest) * 3

        print(f"\nFrame {frame + 1}")
        print(f"Longest stop-free region: {len(longest)} aa")
        print(f"Start position (0-based): {start_nt}")
        print(f"End position (exclusive): {end_nt}")
        print(f"First 40 amino acids: {longest[:40]}")

def inspect_cds_boundary(record):
    """Inspect the boundary between two candidate coding regions."""

    start = 285
    end = 315

    sequence = record.seq[start:end]

    print(f"\nCDS boundary inspection: {record.id}")
    print(f"Nucleotide positions: {start}-{end}")
    print(f"Sequence: {sequence}")

    for position in range(start, end, 3):
        codon = record.seq[position:position + 3]

        if len(codon) == 3:
            amino_acid = codon.translate()
            print(f"{position}: {codon} -> {amino_acid}")

def inspect_all_qualifiers(record):
    """Display all qualifiers associated with each GenBank feature."""

    print(f"\nFull annotation inspection: {record.id}")

    for feature in record.features:
        print(f"\nFeature: {feature.type}")
        print(f"Location: {feature.location}")

        for key, values in feature.qualifiers.items():

            if key == "translation":
                print(f"translation: {len(values[0])} amino acids")

            else: 
                print(f"{key}: {values}")

if __name__ == "__main__":

    df = pd.read_csv("data/processed/rva_vp7_all_metadata.csv")

    unusual = inspect_unsual_lengths(df)

    records = load_records()

    multiple_candidates = []

    for record in records:
        matches = find_vp7_features(record)

        if len(matches) > 1:
            multiple_candidates.append((record.id, len(matches)))

    print(f"\nRecords with multiple VP7 CDS candidates: {len(multiple_candidates)}")

    for accession, count in multiple_candidates[:20]:
        print(f"{accession}: {count} candidate features")

    test_accessions = {
        "D86277.1",
        "MH829417.1",
        "KP170758.1",
    }

    for record in records:
        if record.id in test_accessions:
            for feature in record.features:
                if feature.type == "CDS":

                    report = evaluate_vp7_cds(record, feature)

                    print(f"\nQC report: {record.id}")

                    for key, value in report.items():
                        print(f"{key}: {value}")