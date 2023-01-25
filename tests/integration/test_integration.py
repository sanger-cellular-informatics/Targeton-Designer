import unittest
import sys 
import os

from unittest.mock import patch
from pyfakefs.fake_filesystem_unittest import TestCase
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
            # Use unittest patch to mock sys.argv as if given the commands listed via CLI.
            with patch.object(sys, 'argv', ["./designer.sh", "slicer", "--bed", self.bed_file_path, "--fasta", self.fasta_file_path, "--dir", tmpdir]):
                parsed_input = ParsedInputArguments()
                args = parsed_input.get_args()
                slicer_result = slicer_command(args)
                path_bed = Path(slicer_result.bed)
                path_fasta = Path(slicer_result.fasta)
                # # Check if the files exist.
                self.assertTrue(path_bed.is_file())
                self.assertTrue(path_fasta.is_file())
                # # Check if the files are empty
                self.assertGreater(path_bed.stat().st_size, 0)
                self.assertGreater(path_fasta.stat().st_size, 0)


class TestPrimerIntegration(TestCase):
    def setUp(self):
        self.fasta_file_path = r"./tests/integration/fixtures/slicer_output.fasta"

    def test_PrimerOutput(self):
        with TemporaryDirectory() as tmpdir:
            # Use unittest patch to mock sys.argv as if given the commands listed via CLI.
            with patch.object(sys, 'argv', ["./designer.sh", "primer", "--fasta", self.fasta_file_path, "--dir", tmpdir]):
                parsed_input = ParsedInputArguments()
                args = parsed_input.get_args()
                primer_result = primer_command(args["fasta"], prefix=args["dir"])
                path_bed = Path(primer_result.bed)
                path_csv = Path(primer_result.csv)
                # # Check if the files exist.
                self.assertTrue(path_bed.is_file())
                self.assertTrue(path_csv.is_file())
                # # Check if the files are empty
                self.assertGreater(path_bed.stat().st_size, 0)
                self.assertGreater(path_csv.stat().st_size, 0)
                
class TestIPcressIntegration(TestCase):
    def setUp(self):
        self.use_homo_sapiens = False
        if self.use_homo_sapiens:
            self.fasta_file_path = r"GRCh38.fa"    
        else:
            # self.fasta_file_path = r"./tests/integration/fixtures/slicer_output.fasta"
            self.fasta_file_path = r"./tests/integration/fixtures/fasta_example.fa"
        
        self.p3_output_csv_path = r"./tests/integration/fixtures/p3_output.csv"
       

    def test_iPCRessOutput(self):
        # tmpdir="./tests/integration/fixtures"
        with TemporaryDirectory() as tmpdir:
            if self.use_homo_sapiens:
                self.fasta_file_path = str((Path(tmpdir)/self.fasta_file_path).absolute())
                os.system("wget -cO - http://ftp.ensembl.org/pub/release-106/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz | gunzip > " + self.fasta_file_path)
            # Use unittest patch to mock sys.argv as if given the commands listed via CLI.
            with patch.object(sys, 'argv', ["./designer.sh", "ipcress", "--fasta", self.fasta_file_path, "--dir", tmpdir,"--p3_csv",self.p3_output_csv_path]):
                parsed_input = ParsedInputArguments()
                args = parsed_input.get_args()
                result = ipcress_command(args)
                path_stnd = Path(result.stnd)
                path_err = Path(result.err)
                print("{0} size: {1}".format(path_stnd,path_stnd.stat().st_size))
                print("{0} size: {1}".format(path_err,path_err.stat().st_size))
                # # Check if the files exist.
                self.assertTrue(path_stnd.is_file())
                self.assertTrue(path_err.is_file())
                # # Check if the files are empty
                self.assertGreater(path_stnd.stat().st_size, 0)
                self.assertGreater(path_err.stat().st_size, 0)
            # path_stnd.unlink(missing_ok=True)
            # path_err.unlink(missing_ok=True)

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
        self.bed_file_path = r"./tests/integration/fixtures/bed_example.bed"
        self.use_homo_sapiens = False
        if self.use_homo_sapiens:
            self.fasta_file_path = r"GRCh38.fa"    
        else:
            # self.fasta_file_path = r"./tests/integration/fixtures/slicer_output.fasta"
            self.fasta_file_path = r"./tests/integration/fixtures/fasta_example.fa"

    def test_TDOutput(self):
        with TemporaryDirectory() as tmpdir:
            if self.use_homo_sapiens:
                self.fasta_file_path = str((Path(tmpdir)/self.fasta_file_path).absolute())
                os.system("wget -cO - http://ftp.ensembl.org/pub/release-106/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz | gunzip > " + self.fasta_file_path)
            # Use unittest patch to mock sys.argv as if given the commands listed via CLI.
            with patch.object(sys, 'argv', ["./designer.sh", "design", "--bed", self.bed_file_path, "--fasta", self.fasta_file_path, "--dir", tmpdir]):
                parsed_input = ParsedInputArguments()
                args = parsed_input.get_args()
                result = design_command(args)
                path_bed = Path(result.bed)
                path_fasta = Path(result.fasta)
                path_csv = Path(result.csv)
                path_stnd = Path(result.stnd)
                path_err = Path(result.err)
                
                # # Check if the files exist.
                self.assertTrue(path_bed.is_file())
                self.assertTrue(path_fasta.is_file())
                self.assertTrue(path_csv.is_file())
                self.assertTrue(path_stnd.is_file())
                self.assertTrue(path_err.is_file())
                
                # # Check if the files are empty
                self.assertGreater(path_bed.stat().st_size, 0)
                self.assertGreater(path_fasta.stat().st_size, 0)
                self.assertGreater(path_csv.stat().st_size, 0)
                self.assertGreater(path_stnd.stat().st_size, 0)
                self.assertGreater(path_err.stat().st_size, 0)

if __name__ == '__main__':
    unittest.main()
