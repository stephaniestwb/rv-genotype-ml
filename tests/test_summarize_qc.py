import sys
import unittest
from pathlib import Path

import pandas as pd

# Allow importing modules from src/
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from summarize_qc import (
    merge_metadata_qc,
    assign_preliminary_eligibility
)

class TestSummarizeQC(unittest.TestCase):

    def setUp(self):
        """Create small synthetic metadata and QC tables."""

        self.metadata = pd.DataFrame({
            "accession": ["TEST001.1", "TEST002.1", "TEST003.1"],
            "genotype": ["G1", "G2", "G3"],
            "status": [
                "accepted",
                "accepted",
                "vp7_cds_not_found"
            ]
        })

        self.qc = pd.DataFrame({
            "accession": ["TEST001.1", "TEST002.1"],
            "qc_pass": [True, False]
        })

    def test_valid_merge(self):
        """Valid records should merge correctly."""

        merged = merge_metadata_qc(self.metadata, self.qc)

        self.assertEqual(len(merged), 3)
        self.assertEqual(merged["_merge"].value_counts()["left_only"], 1)

    def test_duplicate_metadata(self):
        """Duplicate metadata accessions should be rejected."""

        metadata = pd.concat(
            [self.metadata, self.metadata.iloc[[0]]],
            ignore_index = True
        )

        with self.assertRaisesRegex(ValueError, "Duplicate"):
            merge_metadata_qc(metadata, self.qc)

    def test_duplicate_qc(self):
        """Duplicate QC accessions should be rejected."""

        qc = pd.concat(
            [self.qc, self.qc.iloc[[0]]],
            ignore_index = True
        )

        with self.assertRaisesRegex(ValueError, "Duplicate"):
            merge_metadata_qc(self.metadata, qc)

    def test_unexpected_missing_qc(self):
        """Missing QC for an accepted record should be rejected."""

        qc = self.qc.iloc[[0]].copy()

        with self.assertRaisesRegex(ValueError, "Unexpected missing QC"):
            merge_metadata_qc(self.metadata, qc)

    def test_preliminaty_eligibility(self):
         """Only accepted sequences passing QC should be eligible."""

         merged = merge_metadata_qc(self.metadata, self.qc)
         result = assign_preliminary_eligibility(merged)

         self.assertEqual(
              result["preliminary_eligible"].tolist(),
              [True, False, False]
         )

    def test_rejected_genotype_with_passing_qc(self):
        """Passing QC alone is insufficient for eligibility."""

        metadata = self.metadata.copy()
        metadata.loc[0, "status"] = "genotype_not_found"

        merged = merge_metadata_qc(metadata, self.qc)
        result = assign_preliminary_eligibility(merged)

        self.assertFalse(result.loc[0, "preliminary_eligible"])

    def test_eligibility_does_not_modify_input(self):
        """Eligibility calculation should not modify its input."""

        merged = merge_metadata_qc(self.metadata, self.qc)
        result = assign_preliminary_eligibility(merged)

        self.assertNotIn("preliminary_eligible", merged.columns)
        self.assertIn("preliminary_eligible", result.columns)

if __name__ == "__main__":
     unittest.main()