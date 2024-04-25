from pyfakefs.fake_filesystem_unittest import TestCase
from parameterized import parameterized
from unittest.mock import patch

from tests.test_data.primer3_output_data import primer3_output_data
from primer.slice_data import SliceData
from primer.primer_pair import \
    PrimerPair, \
    build_primer_pairs, \
    name_primers, \
    capture_primer_details, \
    build_primer_loci


class TestPrimerPairNamePrimers(TestCase):
    @parameterized.expand([
        ('left', '+', 'LibAmpF'),
        ('left', '-', 'LibAmpR'),
        ('right', '+', 'LibAmpR'),
        ('right', '-', 'LibAmpF'),
    ])
    def test_name_primers(self, test_input, strand, expected):
        # act
        actual = name_primers(test_input, strand)

        # assert
        self.assertEqual(actual, expected)


class TestPrimerPairCapturePrimerDetails(TestCase):
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


class TestPrimerPair(TestCase):
    primer3_output_json_data = primer3_output_data

    def setUp(self):
        self.setUpPyfakefs()
        self.input_slice_data = SliceData(
            'slice_name', 'slice_start', 'slice_end', '+', 'slice_chrom', 'bases'
        )

    @patch('primer.primer_pair.build_primer_loci')
    @patch('primer.primer_pair.name_primers')
    @patch('primer.primer_pair.capture_primer_details')
    def test_build_primers_pairs_valid_success(
            self, details_mock, name_mock, loci_mock
    ):
        # arrange
        details_mock.side_effect = [{'side': 'left', 'pair': '0'}, {'side': 'right', 'pair': '0'}, {}]
        name_mock.return_value = 'libamp_name'
        loci_mock.side_effect = [{'id': 'primer_left_0', 'side': 'left'}, {'id': 'primer_right_0', 'side': 'right'}, {}]

        expected = PrimerPair(pair_id="slice_name_0_str", chromosome="", pre_targeton_start="", pre_targeton_end="", product_size="",
                              targeton_id="")
        expected.forward = {}
        expected.reverse = {}

        input_design = {
            "PRIMER_LEFT_0_SEQUENCE":"CAGTGCCAGGACCTCTCCTA",
            "PRIMER_RIGHT_0_SEQUENCE":"TCCCTCTCAGTGGCCATCTT",
            "PRIMER_PAIR_0_PRODUCT_SIZE":""
        }

        # act
        actual = build_primer_pairs(input_design, self.input_slice_data)

        # assert
        self.assertEqual(expected.id, actual[0].id)
        self.assertEqual(expected.forward, actual[0].forward)
        self.assertEqual(expected.reverse, actual[0].reverse)

    @patch('primer.primer_pair.capture_primer_details')
    def test_build_primer_pairs_no_details_empty_(self, details_mock):
        # arrange
        details_mock.return_value = {}
        expected = []
        input_design = {'key_1': 'value'}

        # act
        actual = build_primer_pairs(input_design, self.input_slice_data)

        # assert
        self.assertEqual(expected, actual)
        self.assertEqual(f"{details_mock.call_args}", "call('key_1')")

    @patch('primer.primer_pair.determine_primer_strands')
    @patch('primer.primer_pair.calculate_primer_coords')
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
