from Bio import Entrez, SeqIO
import re
import pandas as pd

Entrez.email = "stampwberth.1@osu.edu"

query = (
    '"Rotavirus A"[Organism] '
    'AND VP7[Title] '
    'AND "complete cds"[Title] '
) 

handle = Entrez.esearch(
    db = "nucleotide",
    term = query,
    retmax = 10000
)

results = Entrez.read(handle)
handle.close()

all_ids = results["IdList"]

sample_positions = [
    0,
    100,
    500,
    1000,
    2000,
    3000,
    4000,
    5000,
    6000,
    7000,
    8000,
    9000,
    9999,
]

sample_ids = [all_ids[i] for i in sample_positions]

print("Search query:")
print(query)

print("\nNumber of matching records:")
print(results["Count"])

print("\nSampled NCBI IDs:")
for position, ncbi_id in zip(sample_positions, sample_ids):
    print(f"Position {position}: {ncbi_id}")

def get_source_qualifiers(record):
    """Return the qualifiers from a GenBank source feature."""

    for feature in record.features:
        if feature.type == "source":
            return feature.qualifiers
    
    return {}

def get_first(qualifiers, key):
    """Return the first value for a qualifier, or None if absent"""

    values = qualifiers.get(key)

    if values: 
        return values[0]
    
    return None

def extract_g_genotype(record, source):
    """
    Extract the RVA G genotype.

    Multiple annotation fields are searched. All detected G genotypes
    must agree for the record to be accepted.
    """

    candidates = [
        get_first(source, "serotype"),
        get_first(source, "genotype"),
        get_first(source, "isolate"),
        get_first(source, "strain"),
        get_first(source, "note"),
        record.description,
    ]

    # Also inspect relevant CDS annotations
    for feature in record.features:
        if feature.type == "CDS":
            candidates.extend([
                get_first(feature.qualifiers, "note"),
                get_first(feature.qualifiers, "gene"),
                get_first(feature.qualifiers, "product"),
            ])

    genotypes = set()

    for text in candidates:
        if not text:
            continue

        matches = re.findall(
            r"\bG(\d+)(?=P|\b)", 
            text,
            flags = re.IGNORECASE
        )

        for number in matches:
            genotypes.add(f"G{number}")

    if len(genotypes) == 1:
        return genotypes.pop(), "accepted"

    if len(genotypes) == 0:
        return None, "genotype_not_found"
    
    return None, "conflicting_genotypes"

def get_vp7_cds(record):
    """Extract the nucleotide sequence of the annotated VP7 CDS."""

    for feature in record.features:
        if feature.type != "CDS":
            continue

        gene = get_first(feature.qualifiers, "gene")
        product = get_first(feature.qualifiers, "product")

        gene_is_vp7 = gene and gene.upper() == "VP7"
        product_is_vp7 = product and "VP7" in product.upper()

        if gene_is_vp7 or product_is_vp7:
            return feature.extract(record.seq)
    
    return None

def parse_record(record):
    """Convert one GenBank record into structured metadata."""

    source = get_source_qualifiers(record)

    vp7_cds = get_vp7_cds(record)

    genotype, genotype_status = extract_g_genotype(record, source)

    if vp7_cds is None:
        status = "vp7_cds_not_found"
    elif genotype_status != "accepted":
        status = genotype_status
    else:
        status = "accepted"
        
    return {
        "accession": record.id,
        "genotype": genotype,
        "status": status,
        "record_length": len(record.seq),
        "cds_length": len(vp7_cds) if vp7_cds is not None else None,
        "organism": record.annotations.get("organism"),
        "host": get_first(source, "host"),
        "country": get_first(source, "geo_loc_name"),
        "collection_date": get_first(source, "collection_date"),
        "isolate": get_first(source, "isolate"),
        "strain": get_first(source, "strain"),
        "serotype": get_first(source, "serotype"),
        "description": record.description,
    }

def print_failed_record(record):
    """Print useful annotations for records with no confident genotype."""

    source = get_source_qualifiers(record)

    genotype, status = extract_g_genotype(record, source)

    if status != "accepted":
        print("\n" + "=" * 80)
        print("FAILED RECORD")
        print("Accession:", record.id)
        print("Status:", status)
        print("Description:", record.description)

        print("\nRelevant source qualifiers:")
    
    
    return

print("\nFetching sampled GenBank records...\n")

fetch_handle = Entrez.efetch(
    db = "nucleotide",
    id = ",".join(sample_ids),
    rettype = "gb",
    retmode = "text"
)

records = list(SeqIO.parse(fetch_handle, "genbank"))
fetch_handle.close()

rows = []

for record in records:
    rows.append(parse_record(record))

df = pd.DataFrame(rows)

print("\nParsed records:\n")

print(
    df[
        [
            "accession",
            "genotype",
            "status",
            "record_length",
            "cds_length",
            "host",
            "country",
            "collection_date",
            "serotype",
        ]
    ].to_string(index = False)
)

print("\nInspecting failed records:")

for record in records:
    print_failed_record(record)
