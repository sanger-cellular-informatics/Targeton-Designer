import unittest
import sys

from unittest.mock import patch
from unittest import TestCase
from pathlib import Path
from tempfile import TemporaryDirectory

from cli import slicer_command, primer_command, ipcress_command, scoring_command, design_command
from utils.arguments_parser import ParsedInputArguments


class TestSlicerIntegration(TestCase):
    def setUp(self):
        self.bed_file_path = r"./tests/integration/fixtures/bed_example.bed"
        self.fasta_file_path = r"./tests/integration/fixtures/fasta_example.fa"

    def test_SlicerOutput(self):
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
        self.fasta_file_path = r"./tests/integration/fixtures/slicer_output.fasta"

    def test_PrimerOutput(self):
        with TemporaryDirectory() as tmpdir:
            # Arrange
            # Use unittest patch to mock sys.argv as if given the commands listed via CLI.
            with patch.object(sys, 'argv', ["./designer.sh", "primer", "--fasta", self.fasta_file_path, "--dir", tmpdir]):
                parsed_input = ParsedInputArguments()
                args = parsed_input.get_args()

                # Act
                primer_result = primer_command(fasta = args["fasta"], prefix = args["dir"])
                path_bed = Path(primer_result.bed)
                path_csv = Path(primer_result.csv)

                # Assert
                self.assertTrue(path_bed.is_file())
                self.assertTrue(path_csv.is_file())
                self.assertGreater(path_bed.stat().st_size, 0)
                self.assertGreater(path_csv.stat().st_size, 0)


class TestIPcressIntegration(TestCase):
    def setUp(self):
        self.fasta_file_path = r"./tests/integration/fixtures/fasta_example.fa"
        self.p3_output_csv_path = r"./tests/integration/fixtures/p3_output.csv"

    def test_iPCRessOutput(self):
        with TemporaryDirectory() as tmpdir:
            # Arrange
            # Use unittest patch to mock sys.argv as if given the commands listed via CLI.
            with patch.object(sys, 'argv', ["./designer.sh", "ipcress", "--fasta", self.fasta_file_path, "--dir", tmpdir, "--p3_csv", self.p3_output_csv_path]):
                parsed_input = ParsedInputArguments()
                args = parsed_input.get_args()

                # Act
                result = ipcress_command(args)
                path_stnd = Path(result.stnd)
                path_err = Path(result.err)

                # Assert
                self.assertTrue(path_stnd.is_file())
                self.assertTrue(path_err.is_file())
                self.assertGreater(path_stnd.stat().st_size, 0)
                self.assertGreater(path_err.stat().st_size, 0)

# class TestScoringIntegration(TestCase):
#     def setUp(self):
#         self.fasta_file_path = r"./tests/integration/fixtures/slicer_output.fasta"
#         self.p3_output_csv_path = r"./tests/integration/fixtures/p3_output.csv"
#         self.ipcress_stnd_path = r"./tests/integration/fixtures/p3_output.csv"

#     def test_ScoringOutput(self):
        # with TemporaryDirectory() as tmpdir:
        #     # Use unittest patch to mock sys.argv as if given the commands listed via CLI.
        #     with patch.object(sys, 'argv', ["./designer.sh", "primer", "--fasta", self.fasta_file_path, "--dir", tmpdir]):
        #         parsed_input = ParsedInputArguments()
        #         args = parsed_input.get_args()
        #         path_tsv = scoring_command(
        #             args['ipcress_file'],
        #             args['scoring_mismatch'],
        #             args['output_tsv'],
        #             args['targeton_csv'],
        #         )
        #         path_tsv = Path(path_tsv)
        #         # # Check if the file exist.
        #         self.assertTrue(path_tsv.is_file())
        #         # # Check if the file are empty
        #         self.assertGreater(path_tsv.stat().st_size, 0)


class TestTargetonDesignerIntegration(TestCase):
    def setUp(self):
        self.fasta_file_path = r"./tests/integration/fixtures/fasta_example.fa"
        self.bed_file_path = r"./tests/integration/fixtures/bed_example.bed"

    def test_TDOutput(self):
        with TemporaryDirectory() as tmpdir:
            # Arrange
            # Use unittest patch to mock sys.argv as if given the commands listed via CLI.
            with patch.object(sys, 'argv', ["./designer.sh", "design", "--bed", self.bed_file_path, "--fasta", self.fasta_file_path, "--dir", tmpdir]):
                parsed_input = ParsedInputArguments()
                args = parsed_input.get_args()

                # Act
                result = design_command(args)

                # Assert
                fields = result.fields()
                fields.remove('dir')
                for field in fields:
                    field_value = getattr(result, field)
                    path = Path(field_value)
                    print(f"Checking file {field} -> {path.name}")
                    self.assertTrue(path.is_file())
                    self.assertGreater(path.stat().st_size, 0)


if __name__ == '__main__':
    unittest.main()
