import sys
import unittest
from pathlib import Path
import tempfile
import pandas as pd

# Allow importing modules from src/
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from analyze_duplicates import (
    sequence_hash,
    find_genotype_conflicts,
    save_sequence_identity,
    save_genotype_conflicts
    )

class TestAnalyzeDuplicates(unittest.TestCase):

    def setUp(self):
        """Create synthetic VP7 nucleotide sequences for testing."""

        self.seq_a = "ATGAAACCT"
        self.seq_b = "atgaaacct"
        self.seq_c = "ATGAAGCCT"

    def test_identical_sequences(self):
        """Identical sequences should produce identical hashes."""

        hash_a = sequence_hash(self.seq_a)
        hash_b = sequence_hash(self.seq_b)

        self.assertEqual(hash_a, hash_b)

    def test_different_sequences(self):
        """A nucleotide substitution should change the hash."""

        hash_a = sequence_hash(self.seq_a)
        hash_c = sequence_hash(self.seq_c)

        self.assertNotEqual(hash_a, hash_c)

    def test_genotype_conflict_detection(self):
        """Identical VP7 sequences with different G labels should be flagged."""

        shared_hash = sequence_hash(self.seq_a)

        sequences = [
            {
                "accession": "TEST001.1",
                "genotype": "G3",
                "status": "accepted",
                "sequence_hash": shared_hash,
            },
            {
                "accession": "TEST002.1",
                "genotype": "G4",
                "status": "accepted",
                "sequence_hash": shared_hash,
            },
            {
                "accession": "TEST003.1",
                "genotype": "G3",
                "status": "accepted",
                "sequence_hash": sequence_hash(self.seq_c),
            }
        ]

        conflicts = find_genotype_conflicts(sequences)

        self.assertEqual(len(conflicts), 1)
        self.assertEqual(conflicts[0]["group_size"], 2)
        self.assertEqual(conflicts[0]["genotypes"], ["G3", "G4"])

    def test_save_Sequence_identity(self):
        """Sequence identity CSV should preserve accession and hash information."""

        sequences = [
            {
                "accession": "TEST001.1",
                "genotype": "G3",
                "status": "accepted",
                "sequence": self.seq_a,
                "sequence_hash": sequence_hash(self.seq_a),  
            },
            {"accession": "TEST002.1",
             "genotype": "G4",
             "status": "accepted",
             "sequence": self.seq_c,
             "sequence_hash": sequence_hash(self.seq_c),
            },
        ]

        with tempfile.TemporaryDirectory() as temp_dir:
            
            output_path = Path(temp_dir) / "sequence_identity.csv"
            save_sequence_identity(sequences, output_path)

            self.assertTrue(output_path.exists())

            df = pd.read_csv(output_path)

            self.assertEqual(len(df), 2)
            self.assertEqual(df["accession"].nunique(), 2)
            self.assertEqual(df["sequence_hash"].nunique(), 2)
            self.assertEqual(df["cds_length"].tolist(), [9, 9])

    def test_save_genotype(self):
        """Conflict CSV should preserve all records in conflicting groups."""

        shared_hash = sequence_hash(self.seq_a)

        sequences = [
            {
                "accession": "TEST001.1",
                "genotype": "G3",
                "status": "accepted",
                "sequence_hash": shared_hash,  
            },
            {
                "accession": "TEST002.1",
                "genotype": "G4",
                "status": "accepted",
                "sequence_hash": shared_hash,
            },
            {
                "accession": "TEST003.1",
                "genotype": None,
                "status": "genotype_not_found",
                "sequence_hash": shared_hash,
            },
        ]

        conflicts = find_genotype_conflicts(sequences)

        with tempfile.TemporaryDirectory() as temp_dir:
                    
                    output_path = Path(temp_dir) / "genotype_conflicts.csv"
                    save_genotype_conflicts(conflicts, output_path)
        
                    self.assertTrue(output_path.exists())
        
                    df = pd.read_csv(output_path)
        
                    self.assertEqual(len(df), 3)
                    self.assertEqual(df["sequence_hash"].nunique(), 1)
                    self.assertEqual(df["group_size"].iloc[0], 3)
                    self.assertEqual(
                        set(df["accession"]),
                        {"TEST001.1", "TEST002.1", "TEST003.1"}
                    )
                    self.assertEqual(df["conflicting_genotypes"].iloc[0], "G3,G4")

if __name__ == "__main__":
    unittest.main()