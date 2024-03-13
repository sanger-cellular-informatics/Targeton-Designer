import unittest
from pyfakefs.fake_filesystem_unittest import TestCase

from tests.test_data.primer3_output_data import primer3_output_data
from primer.primer3 import Primer3
from primer.slice_data import SliceData


class TestPrimer3(TestCase):
    primer3_output_json_data = primer3_output_data

    def setUp(self):
        self.setUpPyfakefs()

    def create_files(self):
        self.fs.create_file('/fwd_primer3_output.json', contents=self.primer3_output_json_data)
        self.fs.create_file('/rev_primer3_output.json', contents=self.primer3_output_json_data)
        self.fs.create_file('/fasta.fa', contents='>region1_1::chr1:5-10(+)\nGTGATCGAGGAGTTCTACAACCAGACATGGGTCCACCGCTA'
                                                  'TGGGGAGAGCATCCTGCCCACCACGCTCACCACGCTCTGGTCCCTCTCAGTGGCCATCTTTTCTGTT'
                                                  'GGGGGCATGATTGGCTCCTTCTCTGTGGGCCTTTTCGTTAACCGCTTTGGCCGGTAAGTAGGAGAGG'
                                                  'TCCTGGCACTGCCCTTGGAGGGCCCATGCCCTCCT')

        self.fs.create_file('/primer3_test_config.json', contents='{\
            "PRIMER_TASK": "pick_cloning_primers",\
            "PRIMER_PICK_LEFT_PRIMER": 1,\
            "PRIMER_PICK_RIGHT_PRIMER": 1,\
            "PRIMER_OPT_SIZE": 20,\
            "PRIMER_MIN_SIZE": 18,\
            "PRIMER_MAX_SIZE": 23,\
            "P3_FILE_FLAG": 1,\
            "SEQUENCE_INCLUDED_REGION": [0,200],\
            "PRIMER_EXPLAIN_FLAG": 1\
        }')


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
        actual = Primer3('/primer3_test_config.json')._primer3_run(input, "")[0].designs

        # assert
        self.assertEqual(actual, expected)

    def test_primer3_initialisation_with_no_user_config_file_as_parameter(self):
        primer = Primer3()

        self.assertEqual(primer._p3_config, './src/primer/primer3.config.json')

    def test_primer3_initialisation_with_user_config_file_as_parameter(self):
        user_config_file = "config.json"

        primer = Primer3(user_config_file)

        self.assertEqual(primer._p3_config, user_config_file)


if __name__ == '__main__':
    unittest.main()
