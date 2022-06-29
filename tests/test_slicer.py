import unittest
import pybedtools
import argparse
from io import StringIO

from designer.slicer import (validate_files, get_slice_data,
    positive_int, parse_args, get_slices)

class TestSlicer(unittest.TestCase):

    def test_check_files(self):
        bed = pybedtools.BedTool('chr1\t100\t250\t.\t.\t+',
            from_string=True)
        fasta = pybedtools.example_filename('test.fa')
        validate_files(bed, fasta)

    def test_get_slice_data(self):
        expected = [
            ('chr1', 50, 260, 'region1_1', '.', '+'),
            ('chr1', 55, 265, 'region1_2', '.', '+'),
            ('chr1', 60, 270, 'region1_3', '.', '+'),
            ('chr1', 65, 275, 'region1_4', '.', '+'),
            ('chr1', 70, 280, 'region1_5', '.', '+'),
            ('chr1', 75, 285, 'region1_6', '.', '+'),
            ('chr1', 80, 290, 'region1_7', '.', '+'),
            ('chr1', 85, 295, 'region1_8', '.', '+'),
            ('chr1', 90, 300, 'region1_9', '.', '+')
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
            ('chr1', 50, 260, 'exon1_1', '.', '-'),
            ('chr1', 55, 265, 'exon1_2', '.', '-'),
            ('chr1', 60, 270, 'exon1_3', '.', '-'),
            ('chr1', 65, 275, 'exon1_4', '.', '-'),
            ('chr1', 70, 280, 'exon1_5', '.', '-'),
            ('chr1', 75, 285, 'exon1_6', '.', '-'),
            ('chr1', 80, 290, 'exon1_7', '.', '-'),
            ('chr1', 85, 295, 'exon1_8', '.', '-'),
            ('chr1', 90, 300, 'exon1_9', '.', '-')
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
        args = parse_args(['bed', 'fasta', '-f5', '50', '--length', '200'])
        self.assertEqual(args.input_bed, 'bed')
        self.assertEqual(args.input_fasta, 'fasta')
        self.assertEqual(args.flank_5, 50)
        self.assertEqual(args.flank_3, 50)
        self.assertEqual(args.length, 200)
        self.assertEqual(args.offset, 5)

    def test_get_slices(self):
        expected_bed = (
            'chr1\t5\t10\tregion1_1\t.\t+\n'
            'chr1\t15\t20\tregion1_2\t.\t+\n'
        )
        expected_fasta = (
            '>region1_1::chr1:5-10(+)\n'
            'AGTCT\n'
            '>region1_2::chr1:15-20(+)\n'
            'ATTTT\n'
        )
        in_bed = StringIO('chr1\t5\t20\t.\t.\t+')
        in_fasta = pybedtools.example_filename('test.fa')
        params = {
            'input_bed': in_bed,
            'input_fasta': in_fasta,
            'flank_5': 0,
            'flank_3': 0,
            'length': 5,
            'offset': 10
        }
        slices = get_slices(params)
        self.assertEqual(expected_bed, slices.head(as_string=True))
        self.assertEqual(expected_fasta, slices.print_sequence())
        expected_bed = (
            'chr1\t5\t10\texon1_1\t.\t-\n'
            'chr1\t15\t20\texon1_2\t.\t-\n'
        )
        expected_fasta = (
            '>exon1_1::chr1:5-10(-)\n'
            'AGACT\n'
            '>exon1_2::chr1:15-20(-)\n'
            'AAAAT\n'
        )
        params['input_bed'] = StringIO('chr1\t5\t20\texon1\t.\t-')
        slices = get_slices(params)
        self.assertEqual(expected_bed, slices.head(as_string=True))
        self.assertEqual(expected_fasta, slices.print_sequence())

if __name__ == '__main__':
    unittest.main()
