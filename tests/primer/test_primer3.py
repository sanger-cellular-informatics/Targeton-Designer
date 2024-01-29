import unittest
import json

from unittest.mock import patch
from pyfakefs.fake_filesystem_unittest import TestCase

from tests.test_data.primer3_output_data import primer3_output_data
from primer.primer3 import Primer3
from utils.parsers import SliceData


class TestPrimer3(TestCase):
    primer3_output_json_data = primer3_output_data

    def setUp(self):
        self.primer = Primer3()
        self.setUpPyfakefs()
        self.input_slice_data = SliceData('slice_name', 'slice_start', 'slice_end', '+', 'slice_chrom', 'bases')

    def create_files(self):
        self.fs.create_file('/fwd_primer3_output.json', contents=self.primer3_output_json_data)
        self.fs.create_file('/rev_primer3_output.json', contents=self.primer3_output_json_data)
        self.fs.create_file(
            '/fasta.fa',
            contents='>region1_1::chr1:5-10(+)\nGTGATCGAGGAGTTCTACAACCAGACATGGGTCCACCGCTA'
            'TGGGGAGAGCATCCTGCCCACCACGCTCACCACGCTCTGGTCCCTCTCAGTGGCCATCTTTTCTGTT'
            'GGGGGCATGATTGGCTCCTTCTCTGTGGGCCTTTTCGTTAACCGCTTTGGCCGGTAAGTAGGAGAGG'
            'TCCTGGCACTGCCCTTGGAGGGCCCATGCCCTCCT',
        )

        self.fs.create_file(
            '/primer3_test_config.json',
            contents='{\
            "PRIMER_TASK": "pick_cloning_primers",\
            "PRIMER_PICK_LEFT_PRIMER": 1,\
            "PRIMER_PICK_RIGHT_PRIMER": 1,\
            "PRIMER_OPT_SIZE": 20,\
            "PRIMER_MIN_SIZE": 18,\
            "PRIMER_MAX_SIZE": 23,\
            "P3_FILE_FLAG": 1,\
            "SEQUENCE_INCLUDED_REGION": [0,212],\
            "PRIMER_EXPLAIN_FLAG": 1\
        }',
        )

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
        expected = [
            SliceData(
                'region1_1',
                '5',
                '10',
                '+',
                'chr1',
                'GTGATCGAGGAGTTCTACAACCAGACATGGGTCCACCGCTATGGGGAGAGCATCCTGCCCACCACGCT'
                'CACCACGCTCTGGTCCCTCTCAGTGGCCATCTTTTCTGTTGGGGGCATGATTGGCTCCTTCTCTGTGG'
                'GCCTTTTCGTTAACCGCTTTGGCCGGTAAGTAGGAGAGGTCCTGGCACTGCCCTTGGAGGGCCCATGC'
                'CCTCCT',
            )
        ]

        # act
        actual = self.primer.read_input_fasta(fasta)

        # assert
        self.assertEqual(actual[0].name, expected[0].name)
        self.assertEqual(actual[0].start, expected[0].start)
        self.assertEqual(actual[0].end, expected[0].end)
        self.assertEqual(actual[0].strand, expected[0].strand)
        self.assertEqual(actual[0].bases, expected[0].bases)

    def test_primer3_design_valid_success(self):
        # arrange
        self.create_files()

        input = [
            SliceData(
                'region1_1',
                '5',
                '10',
                '+',
                'chr1',
                'GTGATCGAGGAGTTCTACAACCAGACATGGGTCCACCGCTATGGGGAGAGCATCCTGCCCACCACGCT'
                'CACCACGCTCTGGTCCCTCTCAGTGGCCATCTTTTCTGTTGGGGGCATGATTGGCTCCTTCTCTGTGG'
                'GCCTTTTCGTTAACCGCTTTGGCCGGTAAGTAGGAGAGGTCCTGGCACTGCCCTTGGAGGGCCCATGC'
                'CCTCCT',
            )
        ]

        expected = {
            'PRIMER_LEFT_EXPLAIN': 'considered 6, low tm 4, ok 2',
            'PRIMER_RIGHT_EXPLAIN': 'considered 6, high tm 6, ok 0',
            'PRIMER_PAIR_EXPLAIN': 'considered 0, ok 0',
            'PRIMER_LEFT_NUM_RETURNED': 0,
            'PRIMER_RIGHT_NUM_RETURNED': 0,
            'PRIMER_INTERNAL_NUM_RETURNED': 0,
            'PRIMER_PAIR_NUM_RETURNED': 0,
        }
        # act

        actual = Primer3('/primer3_test_config.json').primer3_design(input)[0].design

        # assert
        self.assertEqual(actual, expected)

    @patch('primer.primer3.Primer3._build_primers_dict')
    def test_locate_primers_design_success(self, dict_mock):
        # arrange
        dict_mock.return_value = {'region1_1_libamp_name_2': 'build_primer_dict'}

        slice1 = SliceData('slice1', 'start', 'end', 'strand', 'chrom', 'bases')
        slice1.design = {'design_key': 'design_val'}
        slice2 = SliceData('slice2', 'start', 'end', 'strand', 'chrom', 'bases')
        slice2.design = {'design_key': 'design_val'}

        input = [slice1, slice2]

        expected_primers = {'region1_1_libamp_name_2': 'build_primer_dict'}

        # act
        result = self.primer.locate_primers(input)

        # assert
        self.assertEqual(result[0].name, 'slice1')
        self.assertEqual(result[1].name, 'slice2')
        self.assertEqual(result[0].primers, expected_primers)
        self.assertEqual(result[1].primers, expected_primers)
        self.assertEqual(dict_mock.call_count, 2)

    @patch('primer.primer3.Primer3._build_primer_loci')
    @patch('primer.primer3.Primer3.name_primers')
    @patch('primer.primer3.Primer3.capture_primer_details')
    def test_build_primers_dict_valid_success(self, details_mock, name_mock, loci_mock):
        # arrange
        details_mock.return_value = {'pair': '2'}
        name_mock.return_value = 'libamp_name'
        loci_mock.return_value = 'build_primer_dict'
        expected = {'slice_name_libamp_name_2': 'build_primer_dict'}
        input_design = 'design'
        input_primer_keys = {'key_1': 'value'}

        # act
        actual = self.primer._build_primers_dict(input_design, input_primer_keys, self.input_slice_data)

        # assert
        self.assertEqual(expected, actual)

    @patch('primer.primer3.Primer3.capture_primer_details')
    def test_build_primers_dict_no_details_empty(self, details_mock):
        # arrange
        details_mock.return_value = {}
        expected = {}
        input_design = 'design'
        input_primer_keys = {'key_1': 'value'}

        # act
        actual = self.primer._build_primers_dict(input_design, input_primer_keys, self.input_slice_data)

        # assert
        self.assertEqual(expected, actual)
        self.assertEqual(f"{details_mock.call_args}", "call('key_1')")

    # @patch('primer.primer3.Primer3.revcom_reverse_primer')
    @patch('primer.primer3.Primer3.determine_primer_strands')
    @patch('primer.primer3.Primer3.calculate_primer_coords')
    def test_build_primer_loci_with_coords_success(self, coords_mock, strands_mock):
        # arrange
        input_primer = {'penalty': 1, 'side': 'primer_side', 'sequence': 'primer_seq'}
        input_key = 'design_key'
        input_design = {'design_key': 'design_value'}
        input_primer_details = {'field': 'coords', 'side': 'primer_side'}

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
        # revcom_mock.return_value = 'primer_seq'

        # act
        actual = self.primer._build_primer_loci(
            input_primer, input_key, input_design, input_primer_details, self.input_slice_data
        )

        # assert
        self.assertEqual(expected, actual)
        self.assertEqual(f"{coords_mock.call_args}", "call('primer_side', 'design_value', 'slice_start')")
        self.assertEqual(f"{strands_mock.call_args}", "call('primer_side', '+')")

    # self.assertEqual(
    #    f"{revcom_mock.call_args}", "call('primer_seq', 'primer_side_+')")

    def test_build_primer_loci_no_coords_success(self):
        # arrange
        input_primer = {'penalty': 1, 'side': 'primer_side'}
        input_key = 'design_key'
        input_design = {'design_key': 'primer_sequence'}
        input_primer_details = {'field': 'sequence', 'side': 'primer_side'}

        expected = {}
        expected['sequence'] = 'primer_sequence'
        expected['side'] = 'primer_side'
        expected['penalty'] = 1

        # act
        actual = self.primer._build_primer_loci(
            input_primer, input_key, input_design, input_primer_details, self.input_slice_data
        )

        # assert
        self.assertEqual(expected, actual)

    def test_primer3_initialisation_with_no_user_config_file_as_parameter(self):
        primer = Primer3()

        self.assertEqual(primer._config, './src/primer/primer3.config.json')

    def test_primer3_initialisation_with_user_config_file_as_parameter(self):
        user_config_file = "config.json"

        primer = Primer3(user_config_file)

        self.assertEqual(primer._config, user_config_file)


if __name__ == '__main__':
    unittest.main()
