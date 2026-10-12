import sys
import unittest
from pathlib import Path

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation

# Allow importing modules from src/
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from vp7_parser import get_vp7_feature, get_vp7_cds

class TestVP7Parser(unittest.TestCase):

    def setUp(self):
        """Create a synthetic GenBank-like record."""

        self.record = SeqRecord(
            Seq("ATG" + "GCT" * 325 + "TAA"),
            id = "TEST001.1",
            description = "Synthetic rotavirus VP7 sequence"
        )

    def test_vp7_gene_selection(self):
        """Select a CDS annotated with gene=VP7."""

        feature = SeqFeature(
            FeatureLocation(0, len(self.record.seq)),
            type = "CDS",
            qualifiers = {"gene": ["VP7"]}
        ) 

        self.record.features = [feature]

        selected = get_vp7_feature(self.record)

        self.assertIs(selected, feature)

    def test_vp7_product_selection(self):
        """Recognize VP7 from the protein product."""

        feature = SeqFeature(
            FeatureLocation(0, len(self.record.seq)),
            type = "CDS",
            qualifiers = {"product": ["outer capsid glycoprotein VP7"]}
        )

        self.record.features = [feature] 

        selected = get_vp7_feature(self.record)

        self.assertIs(selected, feature)

    def test_vp4_not_selected(self):
        """Do not select a CDS explicitly annotated as VP4."""

        feature = SeqFeature(
            FeatureLocation(0, len(self.record.seq)),
            type = "CDS",
            qualifiers = {"gene": ["VP4"],
                          "product": ["outer capsid spike protein"]}
        ) 

        self.record.featires = [feature]

        selected = get_vp7_feature(self.record)

        self.assertIsNone(selected)

    def test_vp7_cds_Extraction(self):
        """Extract the selected CDS nucleotide sequence."""

        feature = SeqFeature(
            FeatureLocation(0, len(self.record.seq)),
            type = "CDS",
            qualifiers = {"gene": ["VP7"]}
        ) 

        self.record.features = [feature]

        cds = get_vp7_cds(self.record)

        self.assertEqual(str(cds), str(self.record.seq))

    def test_missing_vp7_feature(self):
        """Return None when no VP7 CDS is present."""

        self.record.features = []

        self.assertIsNone(get_vp7_feature(self.record))
        self.assertIsNone(get_vp7_cds(self.record))

    def test_multiple_vp7_features(self):
        """Ambiguous VP7 CDS annotations should not be selected."""

        feature1 = SeqFeature(
            FeatureLocation(0, 600),
            type = "CDS",
            qualifiers = {"gene": ["VP7"]}
        )

        feature2 = SeqFeature(
            FeatureLocation(0, 981),
            type = "CDS",
            qualifiers = {"gene": ["VP7"]}
        )

        self.record.features = [feature1, feature2]

        selected = get_vp7_feature(self.record)

        self.assertIsNone(selected)

if __name__ == "__main__":
    unittest.main()