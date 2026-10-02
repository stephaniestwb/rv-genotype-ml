from collections import Counter

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

if __name__ == "__main__":
    test_sequence = "ATGCNNATG"
    result = check_sequence(test_sequence)

    for metric, value in result.items():
        print(f"{metric}: {value}")