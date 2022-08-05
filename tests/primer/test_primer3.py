from cgi import test
import unittest
import os
import json
import collections

from unittest.mock import patch
from pyfakefs.fake_filesystem_unittest import TestCase

from tests.test_data.primer3_output_data import primer3_output_data
from primer.primer3 import Primer


class TestPrimer3(TestCase):
    primer3_output_json_data = primer3_output_data

    def setUp(self):
        self.primer = Primer()
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
            "SEQUENCE_INCLUDED_REGION": [0,212],\
            "PRIMER_EXPLAIN_FLAG": 1\
        }')

    def test_name_primers_left_fwd_success(self):
        # arrange
        test_input = {'side': 'left'}
        expected = 'LibAmpF'

        # act
        actual = self.primer.name_primers(test_input, '+')

        # assert
        self.assertEqual(actual, expected)

    def test_name_primers_left_rev_success(self):
        # arrange
        test_input = {'side': 'left'}
        expected = 'LibAmpR'

        # act
        actual = self.primer.name_primers(test_input, '-')

        # assert
        self.assertEqual(actual, expected)

    def test_name_primers_right_fwd_success(self):
        # arrange
        test_input = {'side': 'right'}
        expected = 'LibAmpR'

        # act
        actual = self.primer.name_primers(test_input, '+')

        # assert
        self.assertEqual(actual, expected)

    def test_name_primers_right_rev_success(self):
        # arrange
        test_input = {'side': 'right'}
        expected = 'LibAmpF'

        # act
        actual = self.primer.name_primers(test_input, '-')

        # assert
        self.assertEqual(actual, expected)

    def test_capture_primer_details_left_side_success(self):
        # arrange
        test_input = 'primer_left_1_assembly'
        expected = 'left'

        # act
        actual = self.primer.capture_primer_details(test_input)

        # assert
        self.assertEqual(actual['side'], expected)

    def test_capture_primer_details_right_success(self):
        # arrange
        test_input = 'primer_right_3_start'
        expected = 'right'

        # act
        actual = self.primer.capture_primer_details(test_input)

        # assert
        self.assertEqual(actual['side'], expected)

    def test_capture_primer_details_left_pair_success(self):
        # arrange
        test_input = 'primer_left_1_assembly'
        expected = '1'

        # act
        actual = self.primer.capture_primer_details(test_input)

        # assert
        self.assertEqual(actual['pair'], expected)

    def test_capture_primer_details_left_id_success(self):
        # arrange
        test_input = 'primer_left_1_assembly'
        expected = 'primer_left_1'

        # act
        actual = self.primer.capture_primer_details(test_input)

        # assert
        self.assertEqual(actual['id'], expected)

    def test_read_input_fasta_valid_success(self):
        # arrange
        self.create_files()
        fasta = '/fasta.fa'
        expected = [{
            'name'      : 'region1_1',
            'start'     : '5',
            'end'       : '10',
            'strand'    : '+',
            'chrom'     : 'chr1',
            'p3_input'  : {
                            'SEQUENCE_ID'       : 'region1_1',
                            'SEQUENCE_TEMPLATE' : 'GTGATCGAGGAGTTCTACAACCAGACATGGGTCCACCGCTATGGGGAGAGCATCCTGCCCACCACGCT'
                                                  'CACCACGCTCTGGTCCCTCTCAGTGGCCATCTTTTCTGTTGGGGGCATGATTGGCTCCTTCTCTGTGG'
                                                  'GCCTTTTCGTTAACCGCTTTGGCCGGTAAGTAGGAGAGGTCCTGGCACTGCCCTTGGAGGGCCCATGC'
                                                  'CCTCCT'
            }
        }]

        # act
        actual = self.primer.read_input_fasta(fasta)

        # assert
        self.assertEqual(actual, expected)

    def test_primer3_design_valid_success(self):
        # arrange
        self.create_files()

        with open('/primer3_test_config.json', "r") as file:
            config = json.load(file)

        input = [{
            'name'      : 'region1_1',
            'start'     : '5',
            'end'       : '10',
            'strand'    : '+',
            'chrom'     : 'chr1',
            'p3_input'  : {
                            'SEQUENCE_ID'       : 'region1_1',
                            'SEQUENCE_TEMPLATE' : 'GTGATCGAGGAGTTCTACAACCAGACATGGGTCCACCGCTATGGGGAGAGCATCCTGCCCACCACGCT'
                                                  'CACCACGCTCTGGTCCCTCTCAGTGGCCATCTTTTCTGTTGGGGGCATGATTGGCTCCTTCTCTGTGG'
                                                  'GCCTTTTCGTTAACCGCTTTGGCCGGTAAGTAGGAGAGGTCCTGGCACTGCCCTTGGAGGGCCCATGC'
                                                  'CCTCCT'
            }
        }]

        expected = {
                'PRIMER_LEFT_EXPLAIN': 'considered 6, low tm 4, ok 2',
                'PRIMER_RIGHT_EXPLAIN': 'considered 6, high tm 6, ok 0',
                'PRIMER_PAIR_EXPLAIN': 'considered 0, ok 0',
                'PRIMER_LEFT_NUM_RETURNED': 0,
                'PRIMER_RIGHT_NUM_RETURNED': 0,
                'PRIMER_INTERNAL_NUM_RETURNED': 0,
                'PRIMER_PAIR_NUM_RETURNED': 0
            }
        # act
        actual = self.primer.primer3_design(input, config)[0]['design']

        # assert
        self.assertEqual(actual, expected)

    @patch('primer.primer3.Primer.build_primers_dict')
    def test_locate_primers_design_success(self, dict_mock):
        # arrange
        dict_mock.return_value = {'region1_1_libamp_name_2': 'build_primer_dict'}

        input = [{'name': 'slice1', 'design': {'design_key': 'design_val'}},
                 {'name': 'slice2', 'design': {'design_key': 'design_val'}}]

        expected = [{'name': 'slice1', 'primers': {'region1_1_libamp_name_2': 'build_primer_dict'}},
                    {'name': 'slice2', 'primers': {'region1_1_libamp_name_2': 'build_primer_dict'}}]

        # act
        actual = self.primer.locate_primers(input)

        # assert
        self.assertEqual(expected, actual)
        self.assertEqual(dict_mock.call_count, 2)

    @patch('primer.primer3.Primer.build_primer_loci')
    @patch('primer.primer3.Primer.name_primers')
    @patch('primer.primer3.Primer.capture_primer_details')
    def test_build_primers_dict_valid_success(
            self, details_mock, name_mock, loci_mock):
        # arrange
        details_mock.return_value = {'pair': '2'}
        name_mock.return_value = 'libamp_name'
        loci_mock.return_value = 'build_primer_dict'
        expected = {'region1_1_libamp_name_2': 'build_primer_dict'}
        input_design = 'design'
        input_primer_keys = {'key_1': 'value'}
        input_slice_data = {'strand': '+', 'name': 'region1_1'}

        # act
        actual = self.primer.build_primers_dict(
            input_design, input_primer_keys, input_slice_data)

        # assert
        self.assertEqual(expected, actual)
        self.assertEqual(f"{details_mock.call_args}", f"call('key_1')")
        self.assertEqual(f"{name_mock.call_args}", "call({'pair': '2'}, '+')")
        self.assertEqual(
            f"{loci_mock.call_args}", ("call({}, 'key_1', 'design', "
            "{'pair': '2'}, {'strand': '+', 'name': 'region1_1'})"))

    @patch('primer.primer3.Primer.capture_primer_details')
    def test_build_primers_dict_no_details_empty(self, details_mock):
        # arrange
        details_mock.return_value = {}
        expected = {}
        input_design = 'design'
        input_primer_keys = {'key_1': 'value'}
        input_slice_data = {'strand': '+', 'name': 'region1_1'}

        # act
        actual = self.primer.build_primers_dict(input_design, input_primer_keys, input_slice_data)

        # assert
        self.assertEqual(expected, actual)
        self.assertEqual(f"{details_mock.call_args}", "call('key_1')")

    @patch('primer.primer3.Primer.revcom_reverse_primer')
    @patch('primer.primer3.Primer.determine_primer_strands')
    @patch('primer.primer3.Primer.calculate_primer_coords')
    def test_build_primer_loci_with_coords_success(
            self, coords_mock, strands_mock, revcom_mock):
        # arrange
        input_primer = {
            'penalty': 1, 'side': 'primer_side', 'sequence': 'primer_seq'}
        input_key = 'design_key'
        input_design = {'design_key': 'design_value'}
        input_primer_details = {'field': 'coords', 'side': 'primer_side'}
        input_slice_data = {'start': 'slice_start', 'strand': '+'}

        expected = {}
        expected['coords'] = 'design_value'
        expected['side'] = 'primer_side'
        expected['primer_start'] = '100'
        expected['primer_end'] = '250'
        expected['strand'] = 'primer_side_+'
        expected['sequence'] = 'primer_seq'
        expected['penalty'] = 1

        coords_mock.return_value = ['100', '250']
        strands_mock.return_value = 'primer_side_+'
        revcom_mock.return_value = 'primer_seq'

        # act
        actual = self.primer.build_primer_loci(
            input_primer, input_key, input_design,
            input_primer_details, input_slice_data)

        # assert
        self.assertEqual(expected, actual)
        self.assertEqual(
            f"{coords_mock.call_args}",
            "call('primer_side', 'design_value', 'slice_start')")
        self.assertEqual(
                f"{strands_mock.call_args}", "call('primer_side', '+')")
        self.assertEqual(
            f"{revcom_mock.call_args}", "call('primer_seq', 'primer_side_+')")

    def test_build_primer_loci_no_coords_success(self):
        # arrange
        input_primer = {'penalty': 1, 'side': 'primer_side'}
        input_key = 'design_key'
        input_design = {'design_key': 'primer_sequence'}
        input_primer_details = {'field': 'sequence', 'side': 'primer_side'}
        input_slice_data = {'start': 'slice_start', 'strand': '+'}

        expected = {}
        expected['sequence'] = 'primer_sequence'
        expected['side'] = 'primer_side'
        expected['penalty'] = 1

        # act
        actual = self.primer.build_primer_loci(
            input_primer, input_key, input_design,
            input_primer_details, input_slice_data)

        # assert
        self.assertEqual(expected, actual)


if __name__ == '__main__':
    unittest.main()
