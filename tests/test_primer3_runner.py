from cgi import test
import unittest
import os
import json
from unittest.mock import patch
from pyfakefs.fake_filesystem_unittest import TestCase

from runner.primer3_runner import (
    primer3_runner,
    name_primers,
    capture_primer_details,
    locate_primers
)

class TestRunner(TestCase):

    def setUp(self):
        self.setUpPyfakefs()

        self.fs.create_file('/primer3_test_config.json', contents = '{\
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

        primer3_output_json_data = '{"PRIMER_LEFT_EXPLAIN": "considered 6, high tm 4, ok 2", "PRIMER_RIGHT_EXPLAIN": "considered 6, low tm 1, high tm 1, ok 4", "PRIMER_PAIR_EXPLAIN": "considered 5, ok 5", "PRIMER_LEFT_NUM_RETURNED": 5, "PRIMER_RIGHT_NUM_RETURNED": 5, "PRIMER_INTERNAL_NUM_RETURNED": 0, "PRIMER_PAIR_NUM_RETURNED": 5, "PRIMER_PAIR_0_PENALTY": 4.0820471720029445, "PRIMER_LEFT_0_PENALTY": 3.080640677279632, "PRIMER_RIGHT_0_PENALTY": 1.0014064947233123, "PRIMER_LEFT_0_SEQUENCE": "TCCACACAGGATGCCAGG", "PRIMER_RIGHT_0_SEQUENCE": "GGACACTCACCTCAGTTCCTG", "PRIMER_LEFT_0": [0, 18], "PRIMER_RIGHT_0": [211, 21], "PRIMER_LEFT_0_TM": 58.91935932272037, "PRIMER_RIGHT_0_TM": 59.99859350527669, "PRIMER_LEFT_0_GC_PERCENT": 61.111111111111114, "PRIMER_RIGHT_0_GC_PERCENT": 57.142857142857146, "PRIMER_LEFT_0_SELF_ANY_TH": 0.0, "PRIMER_RIGHT_0_SELF_ANY_TH": 0.0, "PRIMER_LEFT_0_SELF_END_TH": 0.0, "PRIMER_RIGHT_0_SELF_END_TH": 0.0, "PRIMER_LEFT_0_HAIRPIN_TH": 43.068834862821745, "PRIMER_RIGHT_0_HAIRPIN_TH": 0.0, "PRIMER_LEFT_0_END_STABILITY": 4.45, "PRIMER_RIGHT_0_END_STABILITY": 3.86, "PRIMER_PAIR_0_COMPL_ANY_TH": 0.0, "PRIMER_PAIR_0_COMPL_END_TH": 3.4095207963081293, "PRIMER_PAIR_0_PRODUCT_SIZE": 212, "PRIMER_PAIR_1_PENALTY": 4.273770098103569, "PRIMER_LEFT_1_PENALTY": 3.2723636033802563, "PRIMER_RIGHT_1_PENALTY": 1.0014064947233123, "PRIMER_LEFT_1_SEQUENCE": "TCCACACAGGATGCCAGGC", "PRIMER_RIGHT_1_SEQUENCE": "GGACACTCACCTCAGTTCCTG", "PRIMER_LEFT_1": [0, 19], "PRIMER_RIGHT_1": [211, 21], "PRIMER_LEFT_1_TM": 62.272363603380256, "PRIMER_RIGHT_1_TM": 59.99859350527669, "PRIMER_LEFT_1_GC_PERCENT": 63.1578947368421, "PRIMER_RIGHT_1_GC_PERCENT": 57.142857142857146, "PRIMER_LEFT_1_SELF_ANY_TH": 3.6618239726775528, "PRIMER_RIGHT_1_SELF_ANY_TH": 0.0, "PRIMER_LEFT_1_SELF_END_TH": 3.6618239726775528, "PRIMER_RIGHT_1_SELF_END_TH": 0.0, "PRIMER_LEFT_1_HAIRPIN_TH": 43.068834862821745, "PRIMER_RIGHT_1_HAIRPIN_TH": 0.0, "PRIMER_LEFT_1_END_STABILITY": 4.85, "PRIMER_RIGHT_1_END_STABILITY": 3.86, "PRIMER_PAIR_1_COMPL_ANY_TH": 0.0, "PRIMER_PAIR_1_COMPL_END_TH": 0.0, "PRIMER_PAIR_1_PRODUCT_SIZE": 212, "PRIMER_PAIR_2_PENALTY": 4.423709731626673, "PRIMER_LEFT_2_PENALTY": 3.080640677279632, "PRIMER_RIGHT_2_PENALTY": 1.3430690543470405, "PRIMER_LEFT_2_SEQUENCE": "TCCACACAGGATGCCAGG", "PRIMER_RIGHT_2_SEQUENCE": "GGACACTCACCTCAGTTCCT", "PRIMER_LEFT_2": [0, 18], "PRIMER_RIGHT_2": [211, 20], "PRIMER_LEFT_2_TM": 58.91935932272037, "PRIMER_RIGHT_2_TM": 58.65693094565296, "PRIMER_LEFT_2_GC_PERCENT": 61.111111111111114, "PRIMER_RIGHT_2_GC_PERCENT": 55.0, "PRIMER_LEFT_2_SELF_ANY_TH": 0.0, "PRIMER_RIGHT_2_SELF_ANY_TH": 0.0, "PRIMER_LEFT_2_SELF_END_TH": 0.0, "PRIMER_RIGHT_2_SELF_END_TH": 0.0, "PRIMER_LEFT_2_HAIRPIN_TH": 43.068834862821745, "PRIMER_RIGHT_2_HAIRPIN_TH": 0.0, "PRIMER_LEFT_2_END_STABILITY": 4.45, "PRIMER_RIGHT_2_END_STABILITY": 3.36, "PRIMER_PAIR_2_COMPL_ANY_TH": 0.0, "PRIMER_PAIR_2_COMPL_END_TH": 0.0, "PRIMER_PAIR_2_PRODUCT_SIZE": 212, "PRIMER_PAIR_3_PENALTY": 4.615432657727297, "PRIMER_LEFT_3_PENALTY": 3.2723636033802563, "PRIMER_RIGHT_3_PENALTY": 1.3430690543470405, "PRIMER_LEFT_3_SEQUENCE": "TCCACACAGGATGCCAGGC", "PRIMER_RIGHT_3_SEQUENCE": "GGACACTCACCTCAGTTCCT", "PRIMER_LEFT_3": [0, 19], "PRIMER_RIGHT_3": [211, 20], "PRIMER_LEFT_3_TM": 62.272363603380256, "PRIMER_RIGHT_3_TM": 58.65693094565296, "PRIMER_LEFT_3_GC_PERCENT": 63.1578947368421, "PRIMER_RIGHT_3_GC_PERCENT": 55.0, "PRIMER_LEFT_3_SELF_ANY_TH": 3.6618239726775528, "PRIMER_RIGHT_3_SELF_ANY_TH": 0.0, "PRIMER_LEFT_3_SELF_END_TH": 3.6618239726775528, "PRIMER_RIGHT_3_SELF_END_TH": 0.0, "PRIMER_LEFT_3_HAIRPIN_TH": 43.068834862821745, "PRIMER_RIGHT_3_HAIRPIN_TH": 0.0, "PRIMER_LEFT_3_END_STABILITY": 4.85, "PRIMER_RIGHT_3_END_STABILITY": 3.36, "PRIMER_PAIR_3_COMPL_ANY_TH": 0.0, "PRIMER_PAIR_3_COMPL_END_TH": 0.0, "PRIMER_PAIR_3_PRODUCT_SIZE": 212, "PRIMER_PAIR_4_PENALTY": 6.997713081139523, "PRIMER_LEFT_4_PENALTY": 3.080640677279632, "PRIMER_RIGHT_4_PENALTY": 3.9170724038598905, "PRIMER_LEFT_4_SEQUENCE": "TCCACACAGGATGCCAGG", "PRIMER_RIGHT_4_SEQUENCE": "GGACACTCACCTCAGTTCC", "PRIMER_LEFT_4": [0, 18], "PRIMER_RIGHT_4": [211, 19], "PRIMER_LEFT_4_TM": 58.91935932272037, "PRIMER_RIGHT_4_TM": 57.08292759614011, "PRIMER_LEFT_4_GC_PERCENT": 61.111111111111114, "PRIMER_RIGHT_4_GC_PERCENT": 57.89473684210526, "PRIMER_LEFT_4_SELF_ANY_TH": 0.0, "PRIMER_RIGHT_4_SELF_ANY_TH": 0.0, "PRIMER_LEFT_4_SELF_END_TH": 0.0, "PRIMER_RIGHT_4_SELF_END_TH": 0.0, "PRIMER_LEFT_4_HAIRPIN_TH": 43.068834862821745, "PRIMER_RIGHT_4_HAIRPIN_TH": 0.0, "PRIMER_LEFT_4_END_STABILITY": 4.45, "PRIMER_RIGHT_4_END_STABILITY": 3.62, "PRIMER_PAIR_4_COMPL_ANY_TH": 0.0, "PRIMER_PAIR_4_COMPL_END_TH": 0.0, "PRIMER_PAIR_4_PRODUCT_SIZE": 212}'
        self.fs.create_file('/fwd_primer3_output.json', contents = primer3_output_json_data)
        self.fs.create_file('/rev_primer3_output.json', contents = primer3_output_json_data)


    def test_runner(self):
        # arrange
        design_input = {
            "SEQUENCE_ID": "ENSE00000893952",
            "SEQUENCE_TEMPLATE": "TCCACACAGGATGCCAGGCCAAGGTGGAGCAAGCGGTGGAGACAGAGCCGGAGCCCGAGCTGCGCCAGCAGACCGAGTGGCAGAGCGGCCAGCGCTGGGAACTGGCACTGGGTCGCTTTTGGGATTACCTGCGCTGGGTGCAGACACTGTCTGAGCAGGTGCAGGAGGAGCTGCTCAGCTCCCAGGTCACCCAGGAACTGAGGTGAGTGTCC"
        }

        # act
        with open(os.devnull, 'w') as devnull:
            with patch('sys.stdout', devnull):
                result = primer3_runner(design_input, '1')

        # assert
        self.assertEqual(result['ENSE00000893952_LibAmpF_0']['sequence'], 'TCCACACAGGATGCCAGG', 'Left primer ok')
        self.assertEqual(result['ENSE00000893952_LibAmpR_0']['sequence'], 'GGACACTCACCTCAGTTCCTG', 'Right primer ok')
    
    def test_name_primers_left_fwd_success(self):
        # arrange
        test_input = { 'side': 'left' }
        expected = 'LibAmpF'

        # act
        actual = name_primers(test_input, '1')

        # assert
        self.assertEqual(actual, expected)

    def test_name_primers_left_rev_success(self):
        # arrange
        test_input = { 'side': 'left' }
        expected = 'LibAmpR'

        # act
        actual = name_primers(test_input, '-1')

        # assert
        self.assertEqual(actual, expected)
    
    def test_name_primers_right_fwd_success(self):
        # arrange
        test_input = { 'side': 'right' }
        expected = 'LibAmpR'

        # act
        actual = name_primers(test_input, '1')

        # assert
        self.assertEqual(actual, expected)

    def test_name_primers_right_rev_success(self):
        # arrange
        test_input = { 'side': 'right' }
        expected = 'LibAmpF'

        # act
        actual = name_primers(test_input, '-1')

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

    def test_locate_primers_fwd_success(self):
        # arrange
        test_input = '/fwd_primer3_output.json'
        expected_f = 'slice1_LibAmpF_0'
        expected_r = 'slice1_LibAmpR_0'
        
        test_dict = {}

        with open(test_input, "r") as p3:
            test_dict = json.load(p3)

        # act        
        actual = locate_primers(test_dict, 'slice1', '1', '1')

        # assert
        self.assertIn(expected_f, actual)
        self.assertEqual(actual[expected_f]['sequence'], 'TCCACACAGGATGCCAGG', 'Primer seq left ok')
        self.assertEqual(actual[expected_r]['sequence'], 'GGACACTCACCTCAGTTCCTG', 'Primer seq right ok')
        self.assertEqual(actual[expected_f]['primer_start'], 1, 'Left start ok')
        self.assertEqual(actual[expected_f]['primer_end'], 19, 'Left end ok')
        self.assertEqual(actual[expected_r]['primer_start'], 191, 'Right start ok')
        self.assertEqual(actual[expected_r]['primer_end'], 212, 'Right end ok')

    def test_locate_primers_rev_success(self):
        # arrange
        test_input = '/rev_primer3_output.json'
        expected_f = 'slice2_LibAmpF_0'
        expected_r = 'slice2_LibAmpR_0'

        with open(test_input, "r") as p3:
            test_dict = json.load(p3)

        # act
        actual = locate_primers(test_dict, 'slice2', '-1', '400')

        # assert 
        self.assertIn(expected_f, actual, 'Rename ok')
        self.assertEqual(actual[expected_f]['sequence'], 'GGACACTCACCTCAGTTCCTG', 'Primer seq left ok')
        self.assertEqual(actual[expected_r]['sequence'], 'TCCACACAGGATGCCAGG', 'Primer seq right ok')
        self.assertEqual(actual[expected_f]['primer_start'], 590, 'Left start ok')
        self.assertEqual(actual[expected_f]['primer_end'], 611, 'Left end ok')
        self.assertEqual(actual[expected_r]['primer_start'], 400, 'Right start ok')
        self.assertEqual(actual[expected_r]['primer_end'], 418, 'Right end ok')

if __name__ == '__main__':
    unittest.main()
