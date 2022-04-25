import unittest
import os
import json

from runner.exonerate import pair_primers, format_ipcress_input, write_ipcress_input

class TestRunner(unittest.TestCase):
    def test_pair_primers(self):
        test_set = {}
        with open('./tests/test_data/primer3_output.json', "r") as data:
            test_set = json.load(data)
        result = pair_primers(test_set)
        
        self.assertEqual(result['0']['left']['sequence'], 'TCCACACAGGATGCCAGG', 'Left 0 ok')
        self.assertEqual(result['0']['right']['sequence'], 'GGACACTCACCTCAGTTCCTG', 'Right 0 ok')
        
        self.assertEqual(result['1']['left']['sequence'], 'TCCACACAGGATGCCAGGC', 'Left 1 ok')
        self.assertEqual(result['1']['right']['sequence'], 'GGACACTCACCTCAGTTCCTG', 'Right 1 ok')

    def test_format_ipcress(self):
        test_set = {}
        with open('./tests/test_data/primer_paired_data.json', "r") as data:
            test_set = json.load(data)
        result = format_ipcress_input('test', test_set)
        
        self.assertEqual(result[0], 'test_0 TCCACACAGGATGCCAGG GGACACTCACCTCAGTTCCTG 200 400', 'First row ok')
        self.assertEqual(result[1], 'test_1 TCCACACAGGATGCCAGGC GGACACTCACCTCAGTTCCTG 200 400', 'Second row ok')

if __name__ == '__main__':
    unittest.main()
