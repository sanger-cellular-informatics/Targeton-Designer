from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import MagicMock, patch

import pysam

from src.utils.bio_utils import get_seq_from_fasta_by_coords


class TestGetSeqFromFastaByCoords(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.temp_dir = TemporaryDirectory()
        self.chr_fasta_path = Path(self.temp_dir.name) / "with_chr.fa"
        self.chr_fasta_path.write_text(">chr1\nACGTTGCAACGTNNTA\n", encoding="ascii")
        pysam.faidx(str(self.chr_fasta_path))

        self.no_chr_fasta_path = Path(self.temp_dir.name) / "without_chr.fa"
        self.no_chr_fasta_path.write_text(
            ">1\nACGTTGCAACGTNNTA\n",
            encoding="ascii",
        )
        pysam.faidx(str(self.no_chr_fasta_path))

    @classmethod
    def tearDownClass(self):
        self.temp_dir.cleanup()

    def test_resolves_unprefixed_input_to_chr_prefixed_fasta_contig(self):
        sequence = get_seq_from_fasta_by_coords(
            self.chr_fasta_path, "1", start=2, end=7
        )

        self.assertEqual(sequence, "CGTTGC")

    def test_fetches_from_an_unprefixed_fasta_contig(self):
        sequence = get_seq_from_fasta_by_coords(
            self.no_chr_fasta_path, "1", start=2, end=7
        )

        self.assertEqual(sequence, "CGTTGC")

    def test_resolves_chr_prefixed_input_to_unprefixed_fasta_contig(self):
        sequence = get_seq_from_fasta_by_coords(
            self.no_chr_fasta_path, "chr1", start=2, end=7
        )

        self.assertEqual(sequence, "CGTTGC")

    def test_fetches_reverse_complement_on_minus_strand(self):
        sequence = get_seq_from_fasta_by_coords(
            self.chr_fasta_path, "1", start=2, end=7, strand="-"
        )

        self.assertEqual(sequence, "GCAACG")

    def test_rejects_zero_coordinates(self):
        with self.assertRaisesRegex(ValueError, "Coordinates must be 1-based"):
            get_seq_from_fasta_by_coords(self.chr_fasta_path, "1", start=0, end=1)

    def test_rejects_reversed_coordinates(self):
        with self.assertRaisesRegex(ValueError, "Coordinates must be 1-based"):
            get_seq_from_fasta_by_coords(self.chr_fasta_path, "1", start=5, end=4)

    def test_rejects_an_invalid_strand(self):
        with self.assertRaisesRegex(ValueError, "strand must be"):
            get_seq_from_fasta_by_coords(
                self.chr_fasta_path, "1", start=1, end=1, strand="forward"
            )

    def test_rejects_a_fasta_file_that_cannot_be_opened(self):
        missing_fasta_path = Path(self.temp_dir.name) / "missing.fa"

        with self.assertRaisesRegex(ValueError, "Could not open FASTA file"):
            get_seq_from_fasta_by_coords(
                missing_fasta_path, "1", start=1, end=1
            )

    def test_rejects_chromosome_not_in_fasta_index(self):
        with self.assertRaisesRegex(ValueError, "not found in the FASTA index"):
            get_seq_from_fasta_by_coords(
                self.chr_fasta_path, "X", start=1, end=1
            )

    def test_rejects_a_range_starting_beyond_the_chromosome_length(self):
        # chr1 is 16 bases long in this test FASTA.
        with self.assertRaisesRegex(ValueError, "exceeds the length"):
            get_seq_from_fasta_by_coords(
                self.chr_fasta_path, "1", start=17, end=20
            )

    def test_wraps_an_error_raised_while_fetching_sequence(self):
        mock_fasta = MagicMock()
        mock_fasta.references = ("chr1",)
        mock_fasta.get_reference_length.return_value = 16
        mock_fasta.fetch.side_effect = OSError("Simulated fetch failure")

        with patch(
            "src.utils.bio_utils.pysam.FastaFile",
            return_value=mock_fasta,
        ):
            with self.assertRaisesRegex(ValueError, "Could not fetch 1:2-7"):
                get_seq_from_fasta_by_coords(
                    self.chr_fasta_path, "1", start=2, end=7
                )

        # 1-based inclusive coordinates 2-7 become 0-based, end-exclusive 1-7.
        mock_fasta.fetch.assert_called_once_with("chr1", 1, 7)


if __name__ == "__main__":
    unittest.main()
