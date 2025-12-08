from unittest import TestCase
from unittest.mock import patch

from primer import flanking


class TestGetFlankedCoordinates(TestCase):
    @patch("primer.flanking._get_chromosome_length")
    def test_no_flanking_returns_original_coords(self, mock_get_chr_len):
        """
        If flanking == 0, _get_flanked_coordinates should return the original
        target region and not even look up chromosome length.
        """
        start, end = flanking._get_flanked_coordinates(
            chromosome="19",
            target_region_start=100,
            target_region_end=200,
            flanking=0,
        )

        self.assertEqual(start, 100)
        self.assertEqual(end, 200)
        mock_get_chr_len.assert_not_called()

    @patch("primer.flanking._get_chromosome_length")
    def test_basic_flanking_within_chr_bounds(self, mock_get_chr_len):
        """
        When flanking is > 0 and stays within chromosome bounds, coordinates
        should expand symmetrically without clamping.
        """
        mock_get_chr_len.return_value = 1000000  # plenty large

        start, end = flanking._get_flanked_coordinates(
            chromosome="19",
            target_region_start=10000,
            target_region_end=10100,
            flanking=200,
        )

        self.assertEqual(start, 10000 - 200)
        self.assertEqual(end, 10100 + 200)


class TestClamping(TestCase):
    @patch("primer.flanking.logger")
    def test_left_clamp_logs_warning(self, mock_logger):
        """
        If the flanked start goes below 1, it should be clamped to 1
        and a warning should be emitted.
        """
        flanked_start_unclamped = -50
        flanked_end_unclamped = 100
        chr_len = 1000

        start_clamped, end_clamped = flanking._clamp_flanked_region(
            chromosome="19",
            flanked_start_unclamped=flanked_start_unclamped,
            flanked_end_unclamped=flanked_end_unclamped,
            chr_len=chr_len,
        )

        self.assertEqual(start_clamped, 1)
        self.assertEqual(end_clamped, flanked_end_unclamped)

        warnings = [str(c) for c in mock_logger.warning.call_args_list]
        self.assertTrue(
            any("extends beyond start of chromosome" in w for w in warnings),
            "Expected a warning about flanking extending beyond chromosome start",
        )

    @patch("primer.flanking.logger")
    def test_right_clamp_logs_warning(self, mock_logger):
        """
        If the flanked end exceeds chromosome length, it should
        be clamped to chr_len and a warning should be emitted.
        """
        chr_len = 1000
        flanked_start_unclamped = 900
        flanked_end_unclamped = 1200  # beyond chromosome end

        start_clamped, end_clamped = flanking._clamp_flanked_region(
            chromosome="19",
            flanked_start_unclamped=flanked_start_unclamped,
            flanked_end_unclamped=flanked_end_unclamped,
            chr_len=chr_len,
        )

        self.assertEqual(start_clamped, flanked_start_unclamped)
        self.assertEqual(end_clamped, chr_len)

        warnings = [str(c) for c in mock_logger.warning.call_args_list]
        self.assertTrue(
            any("extends beyond end of chromosome" in w for w in warnings),
            "Expected a warning about flanking expanding beyond chromosome end",
        )


class TestBuildFlankedSlice(TestCase):
    @patch("primer.flanking.get_seq_from_ensembl_by_coords")
    @patch("primer.flanking._get_flanked_coordinates")
    @patch("primer.flanking.logger")
    def test_build_flanked_slice_with_mode_logs_summary(
        self,
        mock_logger,
        mock_get_flanked_coords,
        mock_get_seq,
    ):
        """
        build_flanked_slice should:
          - compute flanked coordinates
          - fetch the sequence
          - log summary if mode is provided
        """
        mock_get_flanked_coords.return_value = (100, 200)
        mock_get_seq.return_value = "ACGT" * 25  # 100 bp

        flanked_start, flanked_end, flanked_seq = flanking.build_flanked_slice(
            chromosome="19",
            target_region_start=120,
            target_region_end=180,
            strand="+",
            flanking=20,
            mode="REGION",
            note="test-note",
        )

        self.assertEqual(flanked_start, 100)
        self.assertEqual(flanked_end, 200)
        self.assertEqual(flanked_seq, "ACGT" * 25)

        # helper called correctly
        mock_get_flanked_coords.assert_called_once_with(
            chromosome="19",
            target_region_start=120,
            target_region_end=180,
            flanking=20,
        )

        # Ensembl sequence requested with flanked coords
        mock_get_seq.assert_called_once_with(
            chromosome="19",
            start=100,
            end=200,
            strand="+",
        )

        self.assertTrue(mock_logger.info.called)
        logged_args, _ = mock_logger.info.call_args
        first_arg = logged_args[0]
        self.assertIn("Region mode", first_arg)
        self.assertIn("test-note", first_arg)


class TestChromosomeLength(TestCase):
    @patch("primer.flanking.CHR_LENGTHS_GRCh38", {"19": 58617616})
    @patch("primer.flanking.logger")
    def test_get_chromosome_length_known_chr(self, mock_logger):
        """
        _get_chromosome_length should return the length for known chromosomes
        and not log any warnings.
        """
        length = flanking._get_chromosome_length("19")
        self.assertEqual(length, 58617616)
        mock_logger.warning.assert_not_called()

    @patch("primer.flanking.CHR_LENGTHS_GRCh38", {})
    @patch("primer.flanking.logger")
    def test_get_chromosome_length_unknown_chr_logs_warning(self, mock_logger):
        """
        For unknown chromosomes, _get_chromosome_length should return None
        and log a warning.
        """
        length = flanking._get_chromosome_length("99")
        self.assertIsNone(length)
        self.assertTrue(mock_logger.warning.called)
        logged_args, _ = mock_logger.warning.call_args
        self.assertIn("Chromosome '99' not found", logged_args[0])


if __name__ == "__main__":
    unittest.main()
