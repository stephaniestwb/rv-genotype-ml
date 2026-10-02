import re

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