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
        print(result)
        #self.assertTrue(result['PRIMER_LEFT_0_SEQUENCE'], 'GCATCAGTGAGTACAGCATGC')

if __name__ == '__main__':
    unittest.main()
