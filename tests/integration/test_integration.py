import unittest
import sys

from unittest.mock import patch
from pyfakefs.fake_filesystem_unittest import TestCase
from pathlib import Path

from slicer.slicer import Slicer
from cli import slicer_command
from utils.arguments_parser import ParsedInputArguments

import os
print(os.getcwd())


class TestSlicerIntegration(TestCase):
    def setUp(self):
        self.slicer = Slicer()
        # self.setUpPyfakefs()
        self.bed_file_path = r"./tests/integration/fixtures/bed_example.bed"
        self.fasta_file_path = r"./tests/integration/fixtures/fasta_example.fa"
        self.output_dir=r"./tests/integration/td_output"
# def slicer_command(args):
#     validate_files(bed = args['bed'], fasta = args['fasta'])
#     slicer = Slicer()
#     slices = slicer.get_slices(args)

#     return write_slicer_output(args['dir'], slices)
    def test_SlicerOutput(self):
        with patch.object(sys,'argv',["./designer.sh","slicer", "--bed", self.bed_file_path, "--fasta", self.fasta_file_path,"--dir",self.output_dir]):
        # sys.argv= ["./designer.sh","slicer", "--bed", r"./tests/integration/fixtures/bed_example.bed", "--fasta", r"./tests/integration/fixtures/fasta_example.fa"]
            parsed_input = ParsedInputArguments()
            args = parsed_input.get_args()
            # args={"bed":self.bed_file_path,"fasta":self.fasta_file_path,"dir":self.output_dir}
            slicer_result=slicer_command(args)
            path_bed=Path(slicer_result.bed)
            path_fasta=Path(slicer_result.fasta)
            # # Check if the files exist.
            self.assertTrue(path_bed.is_file())
            self.assertTrue(path_fasta.is_file())
            # # Check if the files are empty
            self.assertGreater(path_bed.stat().st_size,0)
            self.assertGreater(path_fasta.stat().st_size,0)
        # print("running")
        


if __name__ == '__main__':
    unittest.main()