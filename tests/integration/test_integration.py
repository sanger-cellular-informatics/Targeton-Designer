import unittest
import sys

from unittest.mock import patch
from unittest import TestCase
from pathlib import Path
from tempfile import TemporaryDirectory

from cli import (
    slicer_command, primer_command, ipcress_command,
    scoring_command, design_command,
    collate_primer_designer_data_command
)
from utils.arguments_parser import ParsedInputArguments
from utils.write_output_files import write_targeton_csv
from designer.output_data_classes import DesignOutputData, PrimerDesignerOutputData, ScoringOutputData
from primer.slice_data import SliceData


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


class TestPrimerIntegration(TestCase):
    def setUp(self):
        self.fasta_file_path = r"./tests/integration/fixtures/test_mask.fa"
        self.config_file_path = r"./tests/primer3_test_config.json"

    def test_primer_output(self):
        with TemporaryDirectory() as tmpdir:
            # Arrange
            # Use unittest patch to mock sys.argv as if given the commands listed via CLI.
            with patch.object(
                sys, 'argv',
                ["./designer.sh", "primer", "--fasta", self.fasta_file_path, "--dir", tmpdir, "--primer3_params", self.config_file_path]
            ):
                parsed_input = ParsedInputArguments()
                args = parsed_input.get_args()

                # Act
                primer_result = primer_command(fasta=args["fasta"], prefix=args["dir"], config=args["primer3_params"])

                path_primer_bed = Path(primer_result.bed)
                path_primer_csv = Path(primer_result.csv)

                # Assert
                self.assertTrue(path_primer_bed.is_file())
                self.assertTrue(path_primer_csv.is_file())
                self.assertGreater(path_primer_bed.stat().st_size, 0)
                self.assertGreater(path_primer_csv.stat().st_size, 0)


class TestIPcressIntegration(TestCase):
    def setUp(self):
        self.fasta_file_path = r"./tests/integration/fixtures/fasta_example.fa"
        self.p3_output_csv_path = r"./tests/integration/fixtures/p3_output.csv"

    def test_ipcress_output(self):
        with TemporaryDirectory() as tmpdir:
            # Arrange
            # Use unittest patch to mock sys.argv as if given the commands listed via CLI.
            with patch.object(sys, 'argv', ["./designer.sh", "ipcress", "--fasta", self.fasta_file_path, "--dir", tmpdir, "--p3_csv", self.p3_output_csv_path]):
                parsed_input = ParsedInputArguments()
                args = parsed_input.get_args()

                # Act
                result = ipcress_command(args)
                path_input = Path(result.input_file)
                path_stnd = Path(result.stnd)
                path_err = Path(result.err)

                # Assert
                self.assertTrue(path_input.is_file())
                self.assertTrue(path_stnd.is_file())
                self.assertTrue(path_err.is_file())
                self.assertGreater(path_input.stat().st_size, 0)
                self.assertGreater(path_stnd.stat().st_size, 0)
                self.assertGreater(path_err.stat().st_size, 0)


class TestTargetonCSVIntegration(TestCase):
    def setUp(self):
        self.ipcress_input_path = r"./tests/integration/fixtures/ipcress_primer_input.txt"
        self.slices = [SliceData('exon1', '100', '250', '+', 'chr1', 'bases')]

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
            with patch.object(sys, 'argv', ["./designer.sh", "primer_designer", "--score_tsv", self.scoring_output_tsv_path, "--dir", tmpdir, "--p3_csv", self.p3_output_csv_path]):
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
