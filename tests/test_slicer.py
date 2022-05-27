import unittest
import pybedtools
import argparse
from io import StringIO

from designer.slicer import (check_files, get_slice_data,
    positive_int, parse_args, main)

class TestSlicer(unittest.TestCase):

    def test_check_files(self):
        bed = pybedtools.BedTool('chr1\t100\t250\t.\t.\t+',
            from_string=True)
        fasta = pybedtools.example_filename('test.fa')
        check_files(bed, fasta)

    def test_get_slice_data(self):
        expected = [
            ('chr1', 50, 260, 1, '.', '+'),
            ('chr1', 55, 265, 1, '.', '+'),
            ('chr1', 60, 270, 1, '.', '+'),
            ('chr1', 65, 275, 1, '.', '+'),
            ('chr1', 70, 280, 1, '.', '+'),
            ('chr1', 75, 285, 1, '.', '+'),
            ('chr1', 80, 290, 1, '.', '+'),
            ('chr1', 85, 295, 1, '.', '+'),
            ('chr1', 90, 300, 1, '.', '+')
        ]
        bed = pybedtools.BedTool('chr1\t100\t250\t.\t.\t+',
            from_string=True)
        params = {
            'flank_5': 50,
            'flank_3': 50,
            'length': 210,
            'offset': 5
        }
        actual = get_slice_data(bed, params)
        self.assertEqual(expected, actual)
        expected = [
            ('chr1', 50, 260, 'exon1', '.', '-'),
            ('chr1', 55, 265, 'exon1', '.', '-'),
            ('chr1', 60, 270, 'exon1', '.', '-'),
            ('chr1', 65, 275, 'exon1', '.', '-'),
            ('chr1', 70, 280, 'exon1', '.', '-'),
            ('chr1', 75, 285, 'exon1', '.', '-'),
            ('chr1', 80, 290, 'exon1', '.', '-'),
            ('chr1', 85, 295, 'exon1', '.', '-'),
            ('chr1', 90, 300, 'exon1', '.', '-')
        ]
        bed = pybedtools.BedTool('chr1\t100\t250\texon1\t.\t-',
            from_string=True)
        actual = get_slice_data(bed, params)
        self.assertEqual(expected, actual)

    def test_positive_int(self):
        self.assertEqual(positive_int('1'), 1)
        error = argparse.ArgumentTypeError
        msg = 'Parameter must be above 0'
        self.assertRaisesRegex(error, msg, positive_int, '-1')
        self.assertRaisesRegex(error, msg, positive_int, '0')

    def test_parse_args(self):
        args = parse_args(['bed', 'fasta','-f5', '50',
            '--length', '200', '--output_slice_bed'])
        self.assertEqual(args.bed, 'bed')
        self.assertEqual(args.fasta, 'fasta')
        self.assertEqual(args.flank_5, 50)
        self.assertEqual(args.flank_3, 0)
        self.assertEqual(args.length, 200)
        self.assertEqual(args.offset, 5)
        self.assertEqual(args.output_slice_bed, 'slices.bed')

    def test_main(self):
        expected = (
            '>1::chr1:5-10(+)\n'
            'AGTCT\n'
            '>1::chr1:15-20(+)\n'
            'ATTTT'
        )
        bed = StringIO('chr1\t5\t20\t.\t.\t+')
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
        expected = (
            '>exon1::chr1:5-10(-)\n'
            'AGACT\n'
            '>exon1::chr1:15-20(-)\n'
            'AAAAT'
        )
        params['bed'] = StringIO('chr1\t5\t20\texon1\t.\t-')
        self.assertEqual(expected, main(params))

if __name__ == '__main__':
    unittest.main()
