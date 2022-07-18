from cgi import test
import unittest
import os
import json

from unittest.mock import patch
from pyfakefs.fake_filesystem_unittest import TestCase

from tests.test_data.primer3_output_data import primer3_output_data
from src.primer.primer3 import (
    name_primers,
    capture_primer_details
)


class TestPrimer3(TestCase):
    primer3_output_json_data = primer3_output_data

    def setUp(self):
        self.setUpPyfakefs()

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
        os.environ["PRIMER3_CONFIG"] = '/primer3_test_config.json'

        self.fs.create_file('/fwd_primer3_output.json', contents=self.primer3_output_json_data)
        self.fs.create_file('/rev_primer3_output.json', contents=self.primer3_output_json_data)

    def test_name_primers_left_fwd_success(self):
        # arrange
        test_input = {'side': 'left'}
        expected = 'LibAmpF'

        # act
        actual = name_primers(test_input, '+')

        # assert
        self.assertEqual(actual, expected)

    def test_name_primers_left_rev_success(self):
        # arrange
        test_input = {'side': 'left'}
        expected = 'LibAmpR'

        # act
        actual = name_primers(test_input, '-')

        # assert
        self.assertEqual(actual, expected)

    def test_name_primers_right_fwd_success(self):
        # arrange
        test_input = {'side': 'right'}
        expected = 'LibAmpR'

        # act
        actual = name_primers(test_input, '+')

        # assert
        self.assertEqual(actual, expected)

    def test_name_primers_right_rev_success(self):
        # arrange
        test_input = {'side': 'right'}
        expected = 'LibAmpF'

        # act
        actual = name_primers(test_input, '-')

        # assert
        self.assertEqual(actual, expected)

    def test_capture_primer_details_left_side_success(self):
        # arrange
        test_input = 'primer_left_1_assembly'
        expected = 'left'

        # act
        actual = capture_primer_details(test_input)

        # assert
        self.assertEqual(actual['side'], expected)

    def test_capture_primer_details_right_success(self):
        # arrange
        test_input = 'primer_right_3_start'
        expected = 'right'

        # act
        actual = capture_primer_details(test_input)

        # assert
        self.assertEqual(actual['side'], expected)

    def test_capture_primer_details_left_pair_success(self):
        # arrange
        test_input = 'primer_left_1_assembly'
        expected = '1'

        # act
        actual = capture_primer_details(test_input)

        # assert
        self.assertEqual(actual['pair'], expected)

    def test_capture_primer_details_left_id_success(self):
        # arrange
        test_input = 'primer_left_1_assembly'
        expected = 'primer_left_1'

        # act
        actual = capture_primer_details(test_input)

        # assert
        self.assertEqual(actual['id'], expected)

###
# TODO: Refactor to new inputs
    # def test_locate_primers_fwd_success(self):
    #     # arrange
    #     test_input = '/fwd_primer3_output.json'
    #     expected_f = 'slice1_LibAmpF_0'
    #     expected_r = 'slice1_LibAmpR_0'
    #     test_dict = {}
    #
    #     with open(test_input, "r") as p3:
    #         test_dict = json.load(p3)
    #
    #     # act
    #     actual = locate_primers(test_dict, 'slice1', '1', '1')
    #
    #     # assert
    #     self.assertIn(expected_f, actual)
    #     self.assertEqual(actual[expected_f]['sequence'], 'TCCACACAGGATGCCAGG', 'Primer seq left ok')
    #     self.assertEqual(actual[expected_r]['sequence'], 'GGACACTCACCTCAGTTCCTG', 'Primer seq right ok')
    #     self.assertEqual(actual[expected_f]['primer_start'], 1, 'Left start ok')
    #     self.assertEqual(actual[expected_f]['primer_end'], 19, 'Left end ok')
    #     self.assertEqual(actual[expected_r]['primer_start'], 191, 'Right start ok')
    #     self.assertEqual(actual[expected_r]['primer_end'], 212, 'Right end ok')
    #
    # def test_locate_primers_rev_success(self):
    #     # arrange
    #     test_input = '/rev_primer3_output.json'
    #     expected_f = 'slice2_LibAmpF_0'
    #     expected_r = 'slice2_LibAmpR_0'
    #
    #     with open(test_input, "r") as p3:
    #         test_dict = json.load(p3)
    #
    #     # act
    #     actual = locate_primers(test_dict, 'slice2', '-1', '400')
    #
    #     # assert
    #     self.assertIn(expected_f, actual, 'Rename ok')
    #     self.assertEqual(actual[expected_f]['sequence'], 'GGACACTCACCTCAGTTCCTG', 'Primer seq left ok')
    #     self.assertEqual(actual[expected_r]['sequence'], 'TCCACACAGGATGCCAGG', 'Primer seq right ok')
    #     self.assertEqual(actual[expected_f]['primer_start'], 590, 'Left start ok')
    #     self.assertEqual(actual[expected_f]['primer_end'], 611, 'Left end ok')
    #     self.assertEqual(actual[expected_r]['primer_start'], 400, 'Right start ok')


if __name__ == '__main__':
    unittest.main()
