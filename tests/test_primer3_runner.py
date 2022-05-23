import unittest
import os

from runner.primer3 import primer3_runner, name_primers, capture_primer_details

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

if __name__ == '__main__':
    unittest.main()
