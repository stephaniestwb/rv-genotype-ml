import sys
import unittest
from pathlib import Path

from Bio.SeqFeature import SeqFeature

# Allow importing modules from src/
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from sequence_qc import (
    check_vp7_cds_integrity,
    check_vp7_annotation,
)

class TestVP7SequenceQC(unittest.TestCase):

    def test_valid_cds(self):
        """A complete, valid CDS should pass QC."""

        sequence = "ATG" + "GCT" * 325 + "TAA"
        flags = check_vp7_cds_integrity(sequence)

        self.assertEqual(flags, [])

    def test_short_cds(self):
        """A truncated CDS should be flagged."""

        sequence = "ATG" + "GCT" * 191 + "TAA"
        flags = check_vp7_cds_integrity(sequence)

        self.assertIn("unexpected_cds_length", flags)

    def test_ambiguous_nucleotides(self):
        """A sequence containing N should be flagged."""

        sequence = "ATG" + "GCT" * 324 + "NNN" + "TAA"
        flags = check_vp7_cds_integrity(sequence)

        self.assertIn("ambiguous_nucleotides", flags)

    def test_internal_stop(self):
        """A premature stop codon should be flagged."""

        sequence = "ATG" + "GCT" * 100 + "TGA" + "GCT" * 224 + "TAA"
        flags = check_vp7_cds_integrity(sequence)

        self.assertIn("internal_stop_codon", flags)

    def test_invalid_reading_frame(self):
        """A CDS with an incomplete codon should be flagged."""

        sequence = "ATG" + "GCT" * 325 + "A" + "TAA"
        flags = check_vp7_cds_integrity(sequence)

        self.assertIn("invalid_reading_frame", flags)

    def test_valid_vp7_annotation(self):
        """A consistent VP7 annotation should not be flagged."""

        feature = SeqFeature(
            type = "CDS",
            qualifiers = {
                "gene": ["VP7"],
                "product": ["outer capsid glycoprotein VP7"],
            },
        )

        flags = check_vp7_annotation(feature)

        self.assertEqual(flags, [])

    def test_vp7_spike_conflict(self):
        """A VP7 gene annotated as a spike protein should be flagged."""

        feature = SeqFeature(
            type = "CDS",
            qualifiers = {
                "gene": ["VP7"],
                "product": ["outer capsid spike protein"],
            },
        )

        flags = check_vp7_annotation(feature)

        self.assertIn("vp7_spike_product_conflict", flags)

    def test_possible_vp4_annotation(self):
        """A product explicitly identifying VP4 should be flagged."""

        feature = SeqFeature(
                    type = "CDS",
                    qualifiers = {
                        "gene": ["VP7"],
                        "product": ["VP4 spike protein"],
                    },
                )
        
        flags = check_vp7_annotation(feature)

        self.assertIn("possible_vp4_annotation", flags)

if __name__ == "__main__":
    unittest.main()