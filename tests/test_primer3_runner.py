import unittest
import os

from runner.primer3 import primer3_runner

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

if __name__ == '__main__':
    unittest.main()
