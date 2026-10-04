import pandas as pd

from build_database import load_records
from vp7_parser import  get_source_qualifiers, extract_g_genotype

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

if __name__ == "__main__":

    df = pd.read_csv("data/processed/rva_vp7_all_metadata.csv")

    unusual = inspect_unsual_lengths(df)
