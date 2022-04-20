import unittest
import os

from runner.primer3 import primer3_runner

class TestRunner(unittest.TestCase):
    @classmethod
    def test_runner(self):
        os.environ["PRIMER3_CONFIG"] = "../tests/primer3_test_config.json"
        design_input = {
            "SEQUENCE_ID": "ENSE00000893952",
            "SEQUENCE_TEMPLATE": "TCCACACAGGATGCCAGGCCAAGGTGGAGCAAGCGGTGGAGACAGAGCCGGAGCCCGAGCTGCGCCAGCAGACCGAGTGGCAGAGCGGCCAGCGCTGGGAACTGGCACTGGGTCGCTTTTGGGATTACCTGCGCTGGGTGCAGACACTGTCTGAGCAGGTGCAGGAGGAGCTGCTCAGCTCCCAGGTCACCCAGGAACTGAGGTGAGTGTCC"
        }
        result = primer3_runner(design_input)
        self.assertTrue(result['primer_left_0']['sequence'], 'TCCACACAGGATGCCAGG')
if __name__ == '__main__':
    unittest.main()
