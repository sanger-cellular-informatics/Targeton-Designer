import unittest
from unittest.mock import patch, MagicMock
from pyfakefs.fake_filesystem_unittest import TestCase

from tests.test_data.primer3_output_data import primer3_output_data
from primer.primer3 import Primer3
from primer.slice_data import SliceData
from config.config import DesignerConfig


class TestPrimer3(TestCase):
    primer3_output_json_data = primer3_output_data

    def setUp(self):
        self.setUpPyfakefs()

        config = MagicMock(spec=DesignerConfig)
        config.stringency_vector = [""]
        self.config = config
        
        self.p3_params = {
            "PRIMER_TASK": "pick_cloning_primers",
            "PRIMER_PICK_LEFT_PRIMER": 1,
            "PRIMER_PICK_RIGHT_PRIMER": 1,
            "PRIMER_OPT_SIZE": 20,
            "PRIMER_MIN_SIZE": 18,
            "PRIMER_MAX_SIZE": 23,
            "P3_FILE_FLAG": 1,
            "SEQUENCE_INCLUDED_REGION" : [0, 200],
            "PRIMER_EXPLAIN_FLAG": 1
        }

    def create_files(self):
        self.fs.create_file('/fwd_primer3_output.json', contents=self.primer3_output_json_data)
        self.fs.create_file('/rev_primer3_output.json', contents=self.primer3_output_json_data)
        self.fs.create_file('/fasta.fa', contents='>region1_1::chr1:5-10(+)\nCTTTTTTCTCTTTCCTTCTGCTTTTGTTTAAAGCGACAAGATGTTGCTCTTTTCCCAGGCTGGAATACAGTGGCATGATCATAGCTCAAGCTCCTGGGCTCAAGTGATCCTCCCGCCTCAGCCTCTCAAGTAGCTAGGACTACAGGCATATCACCACACCAGCGTTTTCTTTGTAGAGGCAGAGTCTCACTCTGTTGCTCAGGCAGGTGTTGAACTCCTGCCTCAAGCAATCCTCCCACCTCAGCCTCCCAGAGCCCTCAAATTATAAGCCACTGTGCTCGGGGCATCCTTTTTGGGGGGTAATCAGCAAACTGAAAAACCTCTTCTTACAACTCCCTATACATTCTCATTCCCAGTATAGAGGAGACTTTTTGTTTTTAAACACTTCCAAAGAATGCAAATTTATAATCCAGAGTATATACATTCTCACTGAATTATTGTACTGTTTCAG')

    def test_primer3_run_valid_success(self):
        # arrange
        self.create_files()

        input = [SliceData('region1_1', '5', '10', '+', 'chr1',
            'GTGATCGAGGAGTTCTACAACCAGACATGGGTCCACCGCTATGGGGAGAGCATCCTGCCCACCACGCT'
            'CACCACGCTCTGGTCCCTCTCAGTGGCCATCTTTTCTGTTGGGGGCATGATTGGCTCCTTCTCTGTGG'
            'GCCTTTTCGTTAACCGCTTTGGCCGGTAAGTAGGAGAGGTCCTGGCACTGCCCTTGGAGGGCCCATGC'
            'CCTCCT'
        )]


        expected = [{
            'PRIMER_INTERNAL': [],
            'PRIMER_LEFT': [],
            'PRIMER_RIGHT': [],
            'PRIMER_PAIR': [],
            'PRIMER_LEFT_EXPLAIN': 'considered 6, low tm 4, ok 2',
            'PRIMER_RIGHT_EXPLAIN': 'considered 6, high tm 6, ok 0',
            'PRIMER_PAIR_EXPLAIN': 'considered 0, ok 0',
            'PRIMER_LEFT_NUM_RETURNED': 0,
            'PRIMER_RIGHT_NUM_RETURNED': 0,
            'PRIMER_INTERNAL_NUM_RETURNED': 0,
            'PRIMER_PAIR_NUM_RETURNED': 0,
            'stringency': ''
        }]
        # act
        actual = Primer3(self.config, self.p3_params)._primer3_run(input, "")[0].designs

        # assert
        self.assertEqual(actual, expected)


if __name__ == '__main__':
    unittest.main()
