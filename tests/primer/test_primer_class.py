from pyfakefs.fake_filesystem_unittest import TestCase
from parameterized import parameterized
from unittest.mock import patch

from collections import defaultdict

from tests.test_data.primer3_output_data import primer3_output_data
from primer.primer3 import Primer3
from primer.slice_data import SliceData
from primer.primer_class import \
    parse_designs_to_primers, \
    name_primers, \
    capture_primer_details, \
    parse_designs_to_primers, \
    build_primers_dict, \
    build_primer_loci


class TestPrimerClassNamePrimers(TestCase):
    @parameterized.expand([
        ({'side': 'left'}, '+', 'LibAmpF'),
        ({'side': 'left'}, '-', 'LibAmpR'),
        ({'side': 'right'}, '+', 'LibAmpR'),
        ({'side': 'right'}, '-', 'LibAmpF'),
    ])
    def test_name_primers(self, test_input, strand, expected):
        # act
        actual = name_primers(test_input, strand)

        # assert
        self.assertEqual(actual, expected)


class TestPrimerClassCapturePrimerDetails(TestCase):
    @parameterized.expand(
        [('primer_left_1_assembly', 'side', 'left'),
        ('primer_right_3_start', 'side', 'right'),
        ('primer_left_1_assembly', 'pair', '1'),
        ('primer_left_1_assembly', 'id', 'primer_left_1'), ])
    def test_capture_primer_details(
            self,
            test_input,
            field,
            expected_value
    ):
        # act
        actual = capture_primer_details(test_input)

        # assert
        self.assertEqual(actual[field], expected_value)


class TestPrimerClass(TestCase):
    primer3_output_json_data = primer3_output_data

    def setUp(self):
        self.primer = Primer3()
        self.setUpPyfakefs()
        self.input_slice_data = SliceData(
            'slice_name', 'slice_start', 'slice_end', '+', 'slice_chrom', 'bases'
        )

    @patch('primer.primer_class.build_primers_dict')
    def test_parse_designs_to_primers_success(self, dict_mock):
        # arrange
        dict_mock.return_value = {'region1_1_libamp_name_2': 'build_primer_dict'}

        slice1 = SliceData('slice1', 'start', 'end', 'strand', 'chrom', 'bases')
        slice1.designs =  [{'design_key': 'design_val', 'stringency': '0.1'}]
        slice2 = SliceData('slice2', 'start', 'end', 'strand', 'chrom', 'bases')
        slice2.designs = [{'design_key': 'design_val', 'stringency': '0.1'}]

        input = [slice1, slice2]

        expected_primers = {'region1_1_libamp_name_2': 'build_primer_dict'}

        # act
        result = parse_designs_to_primers(input)

        # assert
        self.assertEqual(result[0].name, 'slice1')
        self.assertEqual(result[1].name, 'slice2')
        self.assertEqual(result[0].primers, expected_primers)
        self.assertEqual(result[1].primers, expected_primers)
        self.assertEqual(dict_mock.call_count, 2)

    @patch('primer.primer_class.build_primer_loci')
    @patch('primer.primer_class.name_primers')
    @patch('primer.primer_class.capture_primer_details')
    def test_build_primers_dict_valid_success(
            self, details_mock, name_mock, loci_mock
    ):
        # arrange
        details_mock.return_value = {'pair': '2'}
        name_mock.return_value = 'libamp_name'
        loci_mock.return_value = 'build_primer_dict'

        expected = defaultdict(dict)
        expected['slice_name_libamp_name_2_str'] = 'build_primer_dict'

        input_design = 'design'
        input_primer_keys = {'key_1': 'value'}

        # act
        actual = build_primers_dict(
            input_design, input_primer_keys, self.input_slice_data
        )

        # assert
        self.assertEqual(expected, actual)

    @patch('primer.primer_class.capture_primer_details')
    def test_build_primers_dict_no_details_empty(self, details_mock):
        # arrange
        details_mock.return_value = {}
        expected = {}
        input_design = 'design'
        input_primer_keys = {'key_1': 'value'}

        # act
        actual = build_primers_dict(input_design, input_primer_keys, self.input_slice_data)

        # assert
        self.assertEqual(expected, actual)
        self.assertEqual(f"{details_mock.call_args}", "call('key_1')")


    @patch('primer.primer_class.determine_primer_strands')
    @patch('primer.primer_class.calculate_primer_coords')
    def test_build_primer_loci_with_coords_success(
            self, coords_mock, strands_mock
    ):
        # arrange
        input_primer = {
            'penalty': 1, 'side': 'primer_side', 'sequence': 'primer_seq'}
        input_key = 'design_key'
        input_design = {'design_key': 'design_value'}
        input_primer_details = {'field': 'coords', 'side': 'primer_side'}
        input_name = 'name'
        input_id = 'id'

        expected = {'coords': 'design_value',
            'primer': 'name',
            'pair_id': 'id',
            'penalty': 1,
            'primer_end': '250',
            'primer_start': '100',
            'sequence': 'primer_seq',
            'side': 'primer_side',
            'strand': 'primer_side_+',
        }

        coords_mock.return_value = ['100', '250']
        strands_mock.return_value = 'primer_side_+'

        # act
        actual = build_primer_loci(
            input_primer,
            input_key,
            input_design,
            input_primer_details,
            self.input_slice_data,
            input_name,
            input_id,
        )

        # assert
        self.assertEqual(expected, actual)
        self.assertEqual(
            f"{coords_mock.call_args}",
            "call('primer_side', 'design_value', 'slice_start')"
        )

        self.assertEqual(
            f"{strands_mock.call_args}", "call('primer_side', '+')"
        )

    def test_build_primer_loci_no_coords_success(self):
        # arrange
        input_primer = {'penalty': 1, 'side': 'primer_side'}
        input_key = 'design_key'
        input_design = {'design_key': 'primer_sequence'}
        input_primer_details = {'field': 'sequence', 'side': 'primer_side'}
        primer_name = 'name'
        primer_pair_id = 'id'

        expected = {
            'penalty': 1,
            'side': 'primer_side',
            'primer': 'name',
            'sequence': 'primer_sequence',
            'pair_id': 'id',
        }

        # act
        actual = build_primer_loci(
            input_primer,
            input_key,
            input_design,
            input_primer_details,
            self.input_slice_data,
            primer_name,
            primer_pair_id,
        )

        # assert
        self.assertEqual(expected, actual)


if __name__ == '__main__':
    unittest.main()
