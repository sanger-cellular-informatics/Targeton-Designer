import unittest
import sys

from pyfakefs.fake_filesystem_unittest import TestCase
from unittest.mock import patch, Mock

#sys.modules["utils.validate_files"] = Mock()
#sys.modules["slicer.slicer"] = Mock()

from cli import primer_for_ipcress


PRIMER_INPUT_FASTA_PATH = 'tests/test_data/primer_input.fasta'
PARAMS_MIN = '200'
PARAMS_MAX = '300'

class TestCliCommands(TestCase):
    @patch('utils.write_output_files.write_slicer_output')
    def test_primer_for_ipcress(self, write_mock):

        result = primer_for_ipcress(fasta = PRIMER_INPUT_FASTA_PATH, prefix = '!!test_', min = PARAMS_MIN, max = PARAMS_MAX)

        print('TEST result:::::: ')
        print(result)
