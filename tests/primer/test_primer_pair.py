import unittest

from pyfakefs.fake_filesystem_unittest import TestCase
from parameterized import parameterized
from unittest.mock import patch, Mock

from tests.test_data.primer3_output_data import primer3_output_data
from primer.slice_data import SliceData
from primer.build_primer_pairs import _name_primers, _calculate_primer_coords, build_primer_pairs


class TestPrimerPairNamePrimers(TestCase):
    @parameterized.expand([
        ('left', '+', 'LibAmpF'),
        ('left', '-', 'LibAmpR'),
        ('right', '+', 'LibAmpR'),
        ('right', '-', 'LibAmpF'),
    ])
    def test_name_primers(self, test_input, strand, expected):
        # act
        actual = _name_primers(test_input, strand)

        # assert
        self.assertEqual(actual, expected)


class TestPrimerPair(TestCase):
    primer3_output_json_data = primer3_output_data

    def setUp(self):
        self.setUpPyfakefs()
        self.input_slice_data = SliceData(
            name='slice_name',
            start=100, end=200,
            strand='+',
            chromosome='1',
            bases='bases',
            flanking_region=[0, 5, 0, 5],
            exclusion_region=0
        )

    def test_build_primer_pairs_valid_success(self):
        self.assertEqual(False, True)

    def test_build_primer_pairs_no_primer_pairs(self):
        designs = {
            'PRIMER_INTERNAL': [],
            'PRIMER_INTERNAL_NUM_RETURNED': 0,
            'PRIMER_LEFT': [],
            'PRIMER_LEFT_EXPLAIN': 'considered 8528, not in any ok left region 8528, ok 0',
            'PRIMER_LEFT_NUM_RETURNED': 0,
            'PRIMER_PAIR': [],
            'PRIMER_PAIR_EXPLAIN': 'considered 0, ok 0',
            'PRIMER_PAIR_NUM_RETURNED': 0,
            'PRIMER_RIGHT': [],
            'PRIMER_RIGHT_EXPLAIN': 'considered 8528, not in any ok right region 8528, ok 0',
            'PRIMER_RIGHT_NUM_RETURNED': 0
        }

        result = build_primer_pairs(design=designs, slice_data=Mock(), stringency=0.5)

        self.assertEqual(result, [])


class TestCalculatePrimerCoords(unittest.TestCase):

    @parameterized.expand(
        [
            ('left', [55, 20], 44490254, 44490755, '+', (44490309, 44490328)),
            ('right', [168, 19], 44490254, 44490755, '+', (44490404, 44490422)),
            ('left', [39, 20], 44490254, 44490755, '-', (44490697, 44490716)),
            ('right', [167, 19], 44490254, 44490755, '-', (44490588, 44490606)),
        ])
    def test_calculate_primer_coords(
            self,
            side,
            coords,
            slice_start,
            slice_end,
            strand,
            expected_result
    ):
        result = _calculate_primer_coords(side, coords, slice_start, slice_end, strand)

        self.assertEqual(result, expected_result)


if __name__ == '__main__':
    unittest.main()
