import unittest
import sys
import pandas as pd

from unittest.mock import patch
from unittest import TestCase
from pathlib import Path
from tempfile import TemporaryDirectory

from cli import (
    slicer_command, primer_command,
    scoring_command,
    collate_primer_designer_data_command
)
from utils.arguments_parser import ParsedInputArguments
from utils.write_output_files import write_targeton_csv
from designer.output_data_classes import DesignOutputData
from primer.slice_data import SliceData
from config.config import DesignerConfig
from primer.filter.filter_manager import FilterManager


class TestSlicerIntegration(TestCase):
    def setUp(self):
        self.bed_file_path = r"./tests/integration/fixtures/bed_example.bed"
        self.fasta_file_path = r"./tests/integration/fixtures/fasta_example.fa"

    def test_slicer_output(self):
        with TemporaryDirectory() as tmpdir:
            # Arrange
            # Use unittest patch to mock sys.argv as if given the commands listed via CLI.
            with patch.object(sys, 'argv', ["./designer.sh", "slicer", "--bed", self.bed_file_path, "--fasta", self.fasta_file_path, "--dir", tmpdir]):
                parsed_input = ParsedInputArguments()
                args = parsed_input.get_args()

                # Act
                slicer_result = slicer_command(args)
                path_bed = Path(slicer_result.bed)
                path_fasta = Path(slicer_result.fasta)

                # Assert
                self.assertTrue(path_bed.is_file())
                self.assertTrue(path_fasta.is_file())
                self.assertGreater(path_bed.stat().st_size, 0)
                self.assertGreater(path_fasta.stat().st_size, 0)


class BasePrimerIntegrationTest(TestCase):
    """Base class for primer integration tests with shared setup and assertions."""

    def setUp(self):
        # Common config
        self.config_file_path = r"./tests/config_files/test_user_primer3.config.json"
        self.designer_config = r"./tests/config_files/test_user_designer.config.json"

        self.mock_config = {
            "stringency_vector": [1, 2, 3],
            "flanking_region": 150,
            "exclusion_region": 5,
            "csv_column_order": ["col1", "col2", "col3"],
            "filters": {"duplicates": True, "HAP1_variant": True},
        }

    def run_primer_test(self, cli_args_builder):
        """Helper to run primer_command using a temporary directory"""
        with TemporaryDirectory() as tmpdir:
            # Use the tmpdir in the CLI arguments
            cli_args = cli_args_builder(tmpdir)

            # Use unittest patch to mock sys.argv as given using CLI
            with patch.object(sys, "argv", cli_args):
                parsed_input = ParsedInputArguments()
                args = parsed_input.get_args()
                # Act
                primer_result = primer_command(args=args)

                # Shared assertions
                self._assert_primer_outputs(primer_result, args)

    def _assert_primer_outputs(self, primer_result, args):
        """Reusable assertions for primer command output files"""
        # Check BED and CSV exist
        path_primer_bed = Path(primer_result.bed)
        path_primer_csv = Path(primer_result.csv)

        # Assert
        self.assertTrue(path_primer_bed.is_file())
        self.assertTrue(path_primer_csv.is_file())
        self.assertGreater(path_primer_bed.stat().st_size, 0)
        self.assertGreater(path_primer_csv.stat().st_size, 0)

        # Check CSV headers
        expected_csv_headers = DesignerConfig(args={"conf": self.designer_config}).csv_column_order
        df_primers = pd.read_csv(path_primer_csv)
        csv_headers = list(df_primers.columns)

        self.assertEqual(set(csv_headers), set(expected_csv_headers))

        # Tests that the filters have been applied and written to a file
        path_discarded_csv = Path(primer_result.discarded_csv)

        self.assertTrue(path_discarded_csv.is_file())
        self.assertGreater(path_discarded_csv.stat().st_size, 0)

        expected_discarded_headers = expected_csv_headers.copy()
        expected_discarded_headers.append("discard_reason")
        df_discarded = pd.read_csv(path_discarded_csv)
        discarded_headers = list(df_discarded.columns)

        self.assertEqual(set(discarded_headers), set(expected_discarded_headers))

        expected_num_pre_filter = 120
        num_primers = df_primers.shape[0]
        num_discarded = df_discarded.shape[0]

        self.assertGreater(num_primers, 0)
        self.assertGreater(num_discarded, 0)
        self.assertEqual(num_primers + num_discarded, expected_num_pre_filter)

        filters = FilterManager(self.mock_config["filters"])._filters_to_apply
        expected_discard_reasons = [f.reason_discarded for f in filters]
        discard_reasons = df_discarded["discard_reason"].unique().tolist()

        self.assertTrue(set(discard_reasons).issubset(set(expected_discard_reasons)))

        # Tests that ranker was applied and the three optimal primers are output
        path_optimal_primer_pairs_csv = Path(primer_result.optimal_primer_pairs_csv)

        self.assertTrue(path_optimal_primer_pairs_csv.is_file())
        self.assertGreater(path_optimal_primer_pairs_csv.stat().st_size, 0)

        df_optimal_primers = pd.read_csv(path_optimal_primer_pairs_csv)
        optimal_csv_headers = list(df_optimal_primers.columns)
        self.assertEqual(set(optimal_csv_headers), set(expected_csv_headers))
        num_optimal_primers = df_optimal_primers.shape[0]

        # Maximum optimal primer pairs is 3 (6 primers)
        expected_num_optimal_primers = min(6, num_primers)
        self.assertEqual(num_optimal_primers, expected_num_optimal_primers)


class TestPrimerIntegrationFasta(BasePrimerIntegrationTest):
    """Tests primer output using FASTA input."""

    def setUp(self):
        super().setUp()
        self.fasta_file_path = r"./tests/integration/fixtures/test_mask.fa"

    def test_primer_output_from_fasta(self):
        self.run_primer_test(lambda tmpdir: [
            "./designer.sh",
            "primer",
            "--fasta", self.fasta_file_path,
            "--dir", tmpdir,
            "--conf", self.designer_config,
            "--primer3_params", self.config_file_path,
        ])


class TestPrimerIntegrationRegion(BasePrimerIntegrationTest):
    """Tests primer output using region input."""

    def setUp(self):
        super().setUp()
        self.targeton_id = "STEQ"
        self.strand = "+"
        self.region = "chr7:44490254-44490755"

    def test_primer_output_from_region(self):
        self.run_primer_test(lambda tmpdir: [
            "./designer.sh",
            "primer",
            "--targeton_id", self.targeton_id,
            "--region", self.region,
            "--strand", self.strand,
            "--dir", tmpdir,
            "--conf", self.designer_config,
            "--primer3_params", self.config_file_path,
        ])


class TestTargetonCSVIntegration(TestCase):
    def setUp(self):
        self.ipcress_input_path = r"./tests/integration/fixtures/ipcress_primer_input.txt"
        self.slices = [SliceData('exon1', 100, 250, '+', 'chr1', 'bases', 0, 0)]

    def test_write_targeton_csv_output(self):
        with TemporaryDirectory() as tmpdir:
            # Arrange
            cli_input = ["./designer.sh", "generate_targeton_csv", "--primers",
                self.ipcress_input_path, "--dir", tmpdir, ]
            # Use unittest patch to mock sys.argv as if given the commands listed via CLI.
            with patch.object(sys, 'argv', cli_input):
                parsed_input = ParsedInputArguments()
                args = parsed_input.get_args()

                # Act
                result = write_targeton_csv(args['primers'], self.slices,
                                            args['dir'])
                path_csv = Path(result.csv)

                # Assert
                self.assertTrue(path_csv.is_file())
                self.assertGreater(path_csv.stat().st_size, 0)


class TestPrimerDesignerIntegration(TestCase):
    def setUp(self):
        self.scoring_output_tsv_path = r"./tests/integration/fixtures/scoring_output.tsv"
        self.p3_output_csv_path = r"./tests/integration/fixtures/p3_output.csv"

    def test_primer_designer_output(self):
        with TemporaryDirectory() as tmpdir:
            # Arrange
            # Use unittest patch to mock sys.argv as if given the commands listed via CLI.
            with patch.object(sys, 'argv', ["./designer.sh", "design", "--score_tsv", self.scoring_output_tsv_path, "--dir", tmpdir, "--p3_csv", self.p3_output_csv_path]):
                parsed_input = ParsedInputArguments()
                args = parsed_input.get_args()

                # Act
                design_output_data = DesignOutputData(tmpdir)
                design_output_data.p3_csv = args['p3_csv']
                design_output_data.scoring_tsv = args['score_tsv']
                result = collate_primer_designer_data_command(
                    design_output_data,
                    prefix=args['dir']
                    )
                path_json = Path(result.json)
                path_csv = Path(result.csv)

                # Assert
                self.assertTrue(path_json.is_file())
                self.assertTrue(path_csv.is_file())
                self.assertGreater(path_json.stat().st_size, 0)
                self.assertGreater(path_csv.stat().st_size, 0)


class TestScoringIntegration(TestCase):
    def setUp(self):
        self.ipcress_output_path = r"./tests/integration/fixtures/ipcress_output.txt"

    def test_scoring_output(self):
        with TemporaryDirectory() as tmpdir:
            # Arrange
            scoring_output_path = str(Path(tmpdir) / 'scoring_output.tsv')
            cli_input = [
                "./designer.sh", "scoring",
                "--ipcress_file", self.ipcress_output_path,
                "--scoring_mismatch", "5",
                "--output_tsv", scoring_output_path,
            ]
            # Use unittest patch to mock sys.argv as if given the commands listed via CLI.
            with patch.object(sys, 'argv', cli_input):
                parsed_input = ParsedInputArguments()
                args = parsed_input.get_args()

                # Act
                result = scoring_command(
                    args['ipcress_file'], args['scoring_mismatch'], args['output_tsv']
                )
                path_tsv = Path(result.tsv)

                # Assert
                self.assertTrue(path_tsv.is_file())
                self.assertGreater(path_tsv.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
