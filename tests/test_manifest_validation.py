import sys
import unittest
import tempfile
from pathlib import Path

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# Allow importing modules from src/
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from build_database import validate_records

class TestManifestValidation(unittest.TestCase):

    def setUp(self):
        """Create a temporary accession manifest."""

        self.temp_dir = tempfile.TemporaryDirectory()
        self.manifest_path = Path(self.temp_dir.name) / "accessions.txt"

        self.expected = ["Test001.1", "Test002.1", "Test003.1"]

        self.manifest_path.write_text(
            "\n".join(self.expected) + "\n",
            encoding = "utf-8"
            )

    def tearDown(self):
        """Remove temporary files after each test."""

        self.temp_dir.cleanup()

    def make_records(self, accessions):
        """Create snthetic records with specified accession IDs."""

        return [
            SeqRecord(Seq("ATG"), id = accession)
            for accession in accessions
        ]

    def test_valid_manifest(self):
        """All expected accessions should pass validation."""

        records = self.make_records(self.expected)
        validate_records(records, self.manifest_path)

    def test_duplicate_accession(self):
        """Duplicate records should raise an error."""

        records = self.make_records([
            "TEST001.1",
            "TEST002.1",
            "TEST002.1",
            "TEST003.1"
        ])

        with self.assertRaisesRegex(ValueError, "Duplicate"):
            validate_records(records, self.manifest_path)

    def test_missing_accession(self):
        """Missing records should raise an error."""

        records = self.make_records([
            "TEST001.1", 
            "TEST002.1" 
        ])

        with self.assertRaisesRegex(ValueError, "missing"):
            validate_records(records, self.manifest_path)

    def test_unexpected_accession(self):
        """Unexpected records should raise an error."""

        records = self.make_records([
            "TEST001.1",
            "TEST002.1",
            "TEST003.1",
            "TEST004.1"
        ])

        with self.assertRaisesRegex(ValueError, "unexpected"):
            validate_records(records, self.manifest_path)

if __name__ == "__main__":
    unittest.main()