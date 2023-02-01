import unittest
import sys

from unittest import TestCase
from unittest.mock import patch

from cli import primer_for_ipcress


PRIMER_INPUT_FASTA_PATH = 'tests/test_data/primer_input.fasta'
PARAMS_MIN = '200'
PARAMS_MAX = '300'


class TestCliCommands(TestCase):
    @patch('utils.write_output_files.write_slicer_output')
    def test_primer_for_ipcress(self, write_mock):
        expected_result_primers = [
            'ENSE00003571441_HG6_1_LibAmp_0 CCATCTGGGACTCCCTGGG AAAAGGAAGACTGGGTCCTGG 200 300',
            'ENSE00003571441_HG6_1_LibAmp_1 CCATCTGGGACTCCCTGGG AAAAGGAAGACTGGGTCCTG 200 300',
            'ENSE00003571441_HG6_1_LibAmp_2 CCATCTGGGACTCCCTGG AAAAGGAAGACTGGGTCCTGG 200 300',
            'ENSE00003571441_HG6_1_LibAmp_3 CCATCTGGGACTCCCTGGG AAAAGGAAGACTGGGTCCTGGC 200 300',
            'ENSE00003571441_HG6_1_LibAmp_4 CCATCTGGGACTCCCTGG AAAAGGAAGACTGGGTCCTG 200 300'
        ]

        result = primer_for_ipcress(fasta=PRIMER_INPUT_FASTA_PATH,
                                    prefix='!!test_', 
                                    min=PARAMS_MIN, 
                                    max=PARAMS_MAX
                                    )

        self.assertEqual(result, expected_result_primers)
