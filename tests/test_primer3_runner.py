import unittest
import os
import json

from runner.primer3_runner import (
    primer3_runner,
    name_primers,
    capture_primer_details,
    locate_primers
)

class TestRunner(unittest.TestCase):
    def test_runner(self):
        os.environ["PRIMER3_CONFIG"] = "../tests/primer3_test_config.json"

        design_input = {
            "SEQUENCE_ID": "ENSE00000893952",
            "SEQUENCE_TEMPLATE": "TCCACACAGGATGCCAGGCCAAGGTGGAGCAAGCGGTGGAGACAGAGCCGGAGCCCGAGCTGCGCCAGCAGACCGAGTGGCAGAGCGGCCAGCGCTGGGAACTGGCACTGGGTCGCTTTTGGGATTACCTGCGCTGGGTGCAGACACTGTCTGAGCAGGTGCAGGAGGAGCTGCTCAGCTCCCAGGTCACCCAGGAACTGAGGTGAGTGTCC"
        }

        result = primer3_runner(design_input, '1')
        self.assertEqual(result['ENSE00000893952_LibAmpF_0']['sequence'], 'TCCACACAGGATGCCAGG', 'Left primer ok')
        self.assertEqual(result['ENSE00000893952_LibAmpR_0']['sequence'], 'GGACACTCACCTCAGTTCCTG', 'Right primer ok')
    
    def test_primer_orientation(self):
        left = { 'side' : 'left' }
        right = { 'side' : 'right' }

        left_fwd = name_primers(left, '1')
        left_rev = name_primers(left, '-1')
        right_fwd = name_primers(right, '1')
        right_rev = name_primers(right, '-1')

        self.assertEqual(left_fwd, 'LibAmpF', 'Left fwd orientation ok')
        self.assertEqual(left_rev, 'LibAmpR', 'Left rev orientation ok')
        self.assertEqual(right_fwd, 'LibAmpR', 'Right fwd orientation ok')
        self.assertEqual(right_rev, 'LibAmpF', 'Right rev orientation ok')

    def test_capturing_primer_details(self):
        left_value = 'primer_left_1_assembly'
        right_value = 'primer_right_3_start'

        left_output = capture_primer_details(left_value)
        right_output = capture_primer_details(right_value)

        self.assertEqual(left_output['side'], 'left', 'Captured left ok')
        self.assertEqual(right_output['side'], 'right', 'Captured right ok')
        self.assertEqual(left_output['pair'], '1', 'Captured pair ok')
        self.assertEqual(left_output['id'], 'primer_left_1', 'Captured id ok')

    def test_locate_primers(self):
        fwd = './test_data/fwd_primer3_output.json'
        path = os.path.dirname(__file__)
        
        test_dict = {}

        with open(os.path.join(path, fwd), "r") as p3:
            test_dict = json.load(p3)
        
        result = locate_primers(test_dict, 'slice1', '1', '1')    
        self.assertIn('slice1_LibAmpF_0', result, 'Rename ok')
        self.assertEqual(result['slice1_LibAmpF_0']['sequence'], 'TCCACACAGGATGCCAGG', 'Primer seq left ok')
        self.assertEqual(result['slice1_LibAmpR_0']['sequence'], 'GGACACTCACCTCAGTTCCTG', 'Primer seq right ok')
        self.assertEqual(result['slice1_LibAmpF_0']['primer_start'], 1, 'Left start ok')
        self.assertEqual(result['slice1_LibAmpF_0']['primer_end'], 19, 'Left end ok')
        self.assertEqual(result['slice1_LibAmpR_0']['primer_start'], 191, 'Right start ok')
        self.assertEqual(result['slice1_LibAmpR_0']['primer_end'], 212, 'Right end ok')
        
        rev = './test_data/rev_primer3_output.json'
        with open(os.path.join(path, rev), "r") as p3:
            test_dict = json.load(p3)

        result = locate_primers(test_dict, 'slice2', '-1', '400')    
        self.assertIn('slice2_LibAmpF_0', result, 'Rename ok')
        self.assertEqual(result['slice2_LibAmpF_0']['sequence'], 'GGACACTCACCTCAGTTCCTG', 'Primer seq left ok')
        self.assertEqual(result['slice2_LibAmpR_0']['sequence'], 'TCCACACAGGATGCCAGG', 'Primer seq right ok')
        self.assertEqual(result['slice2_LibAmpF_0']['primer_start'], 590, 'Left start ok')
        self.assertEqual(result['slice2_LibAmpF_0']['primer_end'], 611, 'Left end ok')
        self.assertEqual(result['slice2_LibAmpR_0']['primer_start'], 400, 'Right start ok')
        self.assertEqual(result['slice2_LibAmpR_0']['primer_end'], 418, 'Right end ok')

if __name__ == '__main__':
    unittest.main()
