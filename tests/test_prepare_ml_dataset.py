import sys
import unittest
from pathlib import Path
import pandas as pd
import tempfile

# Allow importing modules from src/
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from prepare_ml_dataset import create_unique_sequence_dataset, save_ml_dataset

class TestPrepareMLData(unittest.TestCase):

    def setUp(self):
        """Create a synthetic dataset with duplicate VP7 sequences."""

        self.baseline_df = pd.DataFrame([
            {
                "accession": "TEST002.1",
                "genotype": "G3",
                "sequence_hash": "hash_A",
            },
            {
                "accession": "TEST001.1",
                "genotype": "G3",
                "sequence_hash": "hash_A",
            },
            {
                "accession": "TEST003.1",
                "genotype": "G4",
                "sequence_hash": "hash_B",
            }
        ])

    def test_unique_Sequence_deduplication(self):
        """Duplicate VP7 sequences should retain one representative."""

        unique_df = create_unique_sequence_dataset(self.baseline_df)
        
        self.assertEqual(len(unique_df), 2)
        self.assertEqual(unique_df["sequence_hash"].nunique(), 2)
        self.assertEqual(
            set(unique_df["accession"]), {"TEST001.1", "TEST003.1"}
        )

    def test_save_ml_dataset(self):
        """Saved ML metadata should retain unique sequences and labels."""

        unique_df = create_unique_sequence_dataset(self.baseline_df)

        with tempfile.TemporaryDirectory() as temp_dir:
            output_path = Path(temp_dir) / "ml_metadata.csv"

            save_ml_dataset(unique_df, output_path)

            self.assertTrue(output_path.exists())

            loaded_df = pd.read_csv(output_path)

            self.assertEqual(len(loaded_df), 2)
            self.assertEqual(loaded_df["sequence_hash"].nunique(), 2)
            self.assertEqual(
                set(loaded_df["accession"]),
                {"TEST001.1", "TEST003.1"}
            )
            self.assertEqual(
                dict(zip(
                    loaded_df["accession"],
                    loaded_df["genotype"]
                )),
                {
                    "TEST001.1": "G3",
                    "TEST003.1": "G4",
                }
            )

if __name__ == "__main__":
    unittest.main()