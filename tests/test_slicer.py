import unittest
import pybedtools
import argparse
from io import StringIO

from designer.slicer import (get_slice_coordinates,
    get_slice_sequences, positive_int, parse_args, main)

class TestSlicer(unittest.TestCase):

    def test_get_slice_coordinates(self):
        expected = [
            ('chr1', 50, 260),
            ('chr1', 55, 265),
            ('chr1', 60, 270),
            ('chr1', 65, 275),
            ('chr1', 70, 280),
            ('chr1', 75, 285),
            ('chr1', 80, 290),
            ('chr1', 85, 295),
            ('chr1', 90, 300)
        ]
        bed = pybedtools.BedTool('chr1\t100\t250', from_string=True)
        params = {
            'flank_5': 50,
            'flank_3': 50,
            'length': 210,
            'offset': 5
        }
        actual = get_slice_coordinates(bed, params)
        self.assertEqual(expected, actual)

    def test_get_slice_sequences(self):
        expected = {
            'chr1:5-10': 'AGTCT',
            'chr1:15-20': 'ATTTT'
        }
        bed = pybedtools.BedTool('chr1\t5\t10\nchr1\t15\t20',
            from_string=True)
        actual = get_slice_sequences(bed,
            pybedtools.example_filename('test.fa'))
        self.assertEqual(expected, actual)

    def test_positive_int(self):
        self.assertEqual(positive_int('1'), 1)
        error = argparse.ArgumentTypeError
        msg = 'Parameter must be above 0'
        self.assertRaisesRegex(error, msg, positive_int, '-1')
        self.assertRaisesRegex(error, msg, positive_int, '0')

    def test_parse_args(self):
        args = parse_args(['bed', 'fasta','-f5', '50', '--length', '200'])
        self.assertEqual(args.bed, 'bed')
        self.assertEqual(args.fasta, 'fasta')
        self.assertEqual(args.flank_5, 50)
        self.assertEqual(args.flank_3, 0)
        self.assertEqual(args.length, 200)
        self.assertEqual(args.offset, 5)

    def test_main(self):
        expected = {
            'chr1:5-10': 'AGTCT',
            'chr1:15-20': 'ATTTT'
        }
        bed = StringIO('chr1\t5\t20')
        fasta = pybedtools.example_filename('test.fa')
        params = {
            'bed': bed,
            'fasta': fasta,
            'flank_5': 0,
            'flank_3': 0,
            'length': 5,
            'offset': 10
        }
        self.assertEqual(expected, main(params))

if __name__ == '__main__':
    unittest.main()
