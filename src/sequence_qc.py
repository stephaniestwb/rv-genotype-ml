from collections import Counter
from Bio.Seq import Seq

def check_sequence(sequence):
    """Calculate basic nucleotide sequence QC metrics."""

    sequence = str(sequence).upper()

    length = len(sequence)

    nucleotide_counts = Counter(sequence)

    ambiguous_count = sum(
        count 
        for nucleotide, count in nucleotide_counts.items()
        if nucleotide not in {"A", "C", "G", "T"}
    )

    ambiguous_fraction = (
        ambiguous_count / length
        if length > 0
        else None
    )

    frame_valid = length % 3 == 0

    qc_review = (
        length < 900
        or ambiguous_count > 0
        or not frame_valid
    )

    return {
        "sequence_length": length,
        "ambiguous_count": ambiguous_count,
        "ambiguous_fraction": ambiguous_fraction,
        "frame_valid": frame_valid,
        "qc_review": qc_review,
    }

def check_vp7_annotation(feature):
    """Identify potentially inconsistent VP7 annotations."""

    gene = " ".join(feature.qualifiers.get("gene", [])).lower()
    product = " ".join(feature.qualifiers.get("product", [])).lower()

    flags = []

    if gene == "vp7" and "spike" in product:
        flags.append("vp7_spike_product_conflict")

    if "vp4" in product:
        flags.append("possible_vp4_annotation")

    return flags

def check_vp7_cds_integrity(sequence):
    """
    Check whether a nucleotide sequence is consistent 
    with a near-full-length VP7 coding sequence.

    Returns a list of QC flags.
    """

    seq = str(sequence).upper()
    flags = []

    # Check approximate VP7 CDS length
    if not 900 <= len(seq) <= 1000:
        flags.append("unexpected_cds_length")

    # Check for ambiguous or invalid nucleotides
    if any(base not in "ACGT" for base in seq):
        flags.append("ambiguous_nucleotides")

    # Check reading frame
    if len(seq) % 3 != 0:
        flags.append("invalid_reading_frame")

    # Check start codon
    if not seq.startswith("ATG"):
        flags.append("noncanonical_start_codon")

    # Check terminal stop codon
    if seq[-3:] not in {"TAA", "TAG", "TGA"}:
        flags.append("missing_terminal_stop")

    # Check internal stop codon
    if len(seq) % 3 == 0 and all(base in "ACGT" for base in seq):
        protein = Seq(seq).translate(to_stop = False)

        if "*" in str(protein)[:-1]:
            flags.append("internal_stop_codon")

    return flags

def evaluate_vp7_cds(record, feature):
    """
    Evaluate a GenBank CDS feature for sequence integrity
    and annotation consistency.

    Returns a structured QC report.
    """

    cds = feature.extract(record.seq)
    
    sequence_flags = check_vp7_cds_integrity(cds)
    annotation_flags = check_vp7_annotation(feature)

    all_flags = sequence_flags + annotation_flags

    return {
        "accession": record.id,
        "cds_length": len(cds),
        "sequence_qc_flags": sequence_flags,
        "annotation_qc_flags": annotation_flags,
        "qc_pass": len(all_flags) == 0,
    }

if __name__ == "__main__":
    test_sequence = "ATGCNNATG"
    result = check_sequence(test_sequence)

    for metric, value in result.items():
        print(f"{metric}: {value}")