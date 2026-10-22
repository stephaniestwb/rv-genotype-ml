import sys
import tempfile
import unittest
import pandas as pd

from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

# Allow importing modules from src/
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from export_ml_fasta import (
    load_ml_metadata, 
    find_selected_records, 
    extract_verified_sequences,
    save_ml_fasta,
)

from analyze_duplicates import sequence_hash

class TestExportMLFasta(unittest.TestCase):

    def test_load_ml_metadata(self):
        """Valid metadata should load with the required columns."""

        metadata = pd.DataFrame({
            "accession": ["TEST001.1", "TEST002.1"],
            "genotype": ["G1", "G3"],
            "sequence_hash": ["hash_a", "hash_b"],
        })

        with tempfile.TemporaryDirectory() as temp_dir:
            metadata_path = Path(temp_dir) / "metadata.csv"
            metadata.to_csv(metadata_path, index = False)

            loaded_df = load_ml_metadata(metadata_path)

            self.assertEqual(len(loaded_df), 2)
            self.assertEqual(loaded_df["accession"].tolist(), ["TEST001.1", "TEST002.1"])
            self.assertEqual(loaded_df["genotype"].tolist(), ["G1", "G3"])

    def test_load_ml_metadata_missing_column(self):
        """Missing required column should raise ValueError."""

        metadata = pd.DataFrame({
            "accession": ["TEST001.1"],
            "genotype": ["G1"],
        })
        
        with tempfile.TemporaryDirectory() as temp_dir:
            metadata_path = Path(temp_dir) / "metadata.csv"
            metadata.to_csv(metadata_path, index = False)
        
            with self.assertRaises(ValueError):
                load_ml_metadata(metadata_path)

    def test_find_selected_records(self):
        """Only requested GenBank accessions should be returned."""

        records = [
            SeqRecord(
                Seq("ATGAAATAG"),
                id = "TEST001.1",
                description = "Test record 1",
                annotations = {"molecule_type": "DNA"}
            ),
            SeqRecord(
                Seq("ATGCCCTAG"),
                id = "TEST002.1",
                description = "Test record 2",
                annotations = {"molecule_type": "DNA"}
            ),
            SeqRecord(
                Seq("ATGGGGTAG"),
                id = "TEST003.1",
                description = "Test record 3",
                annotations = {"molecule_type": "DNA"}
            ),
        ]

        with tempfile.TemporaryDirectory() as temp_dir:
            raw_dir = Path(temp_dir)
            gb_path = raw_dir / "batch_00000.gb"

            SeqIO.write(records, gb_path, "genbank")

            selected = {"TEST001.1", "TEST003.1"}

            matched = find_selected_records(raw_dir, selected)

            self.assertEqual(len(matched), 2)
            self.assertEqual(set(matched), selected)
            self.assertNotIn("TEST002.1", matched)

    def test_extract_verified_sequences(self):
        """A VP7 CDS with the correct hash should be accepted."""

        sequence = "ATGAAACCCGGGTAA"

        record = SeqRecord(
            Seq(sequence),
            id = "TEST001.1",
            description = "Synthetic VP7 record"
        )

        record.features = [
            SeqFeature(
                FeatureLocation(0, len(sequence)),
                type = "CDS",
                qualifiers = {"gene": ["VP7"]}
            )
        ]

        metadata = pd.DataFrame({
            "accession": ["TEST001.1"],
            "genotype": ["G1"],
            "sequence_hash": [sequence_hash(sequence)]
        })

        matched_records = {"TEST001.1": record}

        verified = extract_verified_sequences(metadata, matched_records)

        self.assertEqual(len(verified), 1)
        self.assertEqual(str(verified["TEST001.1"]), sequence)

    def test_extract_verified_sequences_hash_mismatch(self):
        """A VP7 CDS with an incorrect hash should be rejected."""

        sequence = "ATGAAACCCGGGTAA"

        record = SeqRecord(
            Seq(sequence),
            id = "TEST0001.1",
            description = "Synthetic VP7 record",
        )
        record.features = [
            SeqFeature(
                FeatureLocation(0, len(sequence)),
                type = "CDS",
                qualifiers = {"gene": ["VP7"]}
            )
        ]

        metadata = pd.DataFrame({
            "accession": ["TEST001.1"],
            "genotype": ["G1"],
            "sequence_hash": [sequence_hash("ATGTTTCCCGGGAAA")]
        })

        matched_records = {"TEST001.1": record}

        with self.assertRaisesRegex(ValueError, "Sequence hash mismatch"): 
            extract_verified_sequences(metadata, matched_records)

    def test_save_ml_fasta(self):
        """FASTA export should preserve accessions, genotypes and sequences."""

        sequences = {
            "TEST001.1": Seq("ATGAAACCCGGGTAA"),
            "TEST002.1": Seq("ATGTTTGGGCCCTAA"),
        }

        metadata = pd.DataFrame({
            "accession": ["TEST001.1", "TEST002.1"],
            "genotype": ["G1", "G3"],
            "sequence_hash": [
                sequence_hash(str(sequences["TEST001.1"])),
                sequence_hash(str(sequences["TEST002.1"]))
            ]
        })

        with tempfile.TemporaryDirectory() as temp_dir:
            output_path = Path(temp_dir) / "test_ml.fasta"

            save_ml_fasta(metadata, sequences, output_path)

            self.assertTrue(output_path.exists())

            loaded_records = list(SeqIO.parse(output_path, "fasta"))

            self.assertEqual(len(loaded_records), 2)
            self.assertEqual(
                [record.id for record in loaded_records],
                ["TEST001.1", "TEST002.1"]
            )
            self.assertEqual(
                [str(record.seq) for record in loaded_records],
                ["ATGAAACCCGGGTAA", "ATGTTTGGGCCCTAA"]
            )
            self.assertEqual(
                [record.description for record in loaded_records],
                ["TEST001.1 genotype=G1", "TEST002.1 genotype=G3"]
            )