from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import pysam

from src.utils.bio_utils import get_seq_from_fasta_by_coords


class TestGetSeqFromFastaByCoords(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.temp_dir = TemporaryDirectory()
        self.fasta_path = Path(self.temp_dir.name) / "example.fa"
        self.fasta_path.write_text(">chr1\nACGTTGCAACGTNNTA\n", encoding="ascii")
        pysam.faidx(str(self.fasta_path))

    @classmethod
    def tearDownClass(self):
        self.temp_dir.cleanup()

    def test_fetches_a_one_based_inclusive_range(self):
        sequence = get_seq_from_fasta_by_coords(
            self.fasta_path, "chr1", start=2, end=7
        )

        self.assertEqual(sequence, "CGTTGC")

    def test_fetches_reverse_complement_on_minus_strand(self):
        sequence = get_seq_from_fasta_by_coords(
            self.fasta_path, "chr1", start=2, end=7, strand="-"
        )

        self.assertEqual(sequence, "GCAACG")

    def test_rejects_zero_coordinates(self):
        with self.assertRaisesRegex(ValueError, "Coordinates must be 1-based"):
            get_seq_from_fasta_by_coords(self.fasta_path, "chr1", start=0, end=1)

    def test_rejects_reversed_coordinates(self):
        with self.assertRaisesRegex(ValueError, "Coordinates must be 1-based"):
            get_seq_from_fasta_by_coords(self.fasta_path, "chr1", start=5, end=4)

    def test_rejects_an_invalid_strand(self):
        with self.assertRaisesRegex(ValueError, "strand must be"):
            get_seq_from_fasta_by_coords(
                self.fasta_path, "chr1", start=1, end=1, strand="forward"
            )

    def test_rejects_a_missing_fasta(self):
        missing_path = self.fasta_path.parent / "missing.fa"

        with self.assertRaisesRegex(FileNotFoundError, "FASTA file not found"):
            get_seq_from_fasta_by_coords(missing_path, "chr1", start=1, end=1)

    def test_wraps_an_unknown_chromosome_error(self):
        with self.assertRaisesRegex(ValueError, "Could not fetch chrX:1-1"):
            get_seq_from_fasta_by_coords(self.fasta_path, "chrX", start=1, end=1)


if __name__ == "__main__":
    unittest.main()
