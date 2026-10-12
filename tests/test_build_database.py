import sys
import unittest
from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

# Allow importing modules from src/
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from build_database import process_records

class TestBuildDatabase(unittest.TestCase):

    def setUp(self):
        """Create a synthetic VP7 GenBank-like record."""

        self.record = SeqRecord(
            Seq("ATG" + "GTC" * 325 + "TAA"),
            id = "TEST001.1",
            description = "Rotavirus A G1 VP7 complete cds"
        )

        source = SeqFeature(
            FeatureLocation(0, len(self.record.seq)),
            type = "source",
            qualifiers = {
                "organism": ["Rotavirus A"],
                "strain": ["RVA/Human/Test/G1P[8]"],
            }
        )

        cds = SeqFeature(
            FeatureLocation(0, len(self.record.seq)),
            type = "CDS",
            qualifiers = {
                "gene": ["VP7"],
                "product": ["outer capsid glycoprotein VP7"]
            }
        )

        self.record.features = [source, cds]

    def test_process_valid_record(self):
        """A valid VP7 recorrd should produce metadata and QC metrics."""

        df = process_records([self.record])

        self.assertEqual(len(df), 1)
        self.assertEqual(df.iloc[0]["accession"], "TEST001.1")
        self.assertEqual(df.iloc[0]["genotype"], "G1")

    def test_process_missing_vp7(self):
        """A record without VP7 CDS should retain its metadata."""

        self.record.features = []
        self.assertEqual(len(self.record.features), 0)

        df = process_records([self.record])

        self.assertEqual(len(df), 1)
        self.assertEqual(df.iloc[0]["status"], "vp7_cds_not_found")

    def test_parser_missing_vp7(self):
        """Parser sould reject a record without a VP7 CDS."""

        self.record.features = [
            feature for feature in self.record.features
            if feature.type == "source"
        ]

        from vp7_parser import parse_record

        metadata = parse_record(self.record)

        print("\nMetadata without VP7 CDS:", metadata)

        self.assertEqual(
            metadata["status"],
            "vp7_cds_not_found"
        )

if __name__ == "__main__":
    unittest.main()