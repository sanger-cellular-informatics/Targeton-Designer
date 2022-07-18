import unittest
import pybedtools
import argparse

from io import StringIO
from unittest.mock import patch
from pyfakefs.fake_filesystem_unittest import TestCase

from src.slicer.slicer import (get_slice_data, get_slices, decrement_one_based_starts)


class TestSlicer(TestCase):

    def setUp(self):
        self.setUpPyfakefs()
        self.bed_file_data = 'chr1\t100\t250\texon1\t.\t+'
        self.fasta_file_data = '>region1_1::chr1:5-10(+)\nAGTCT\n>region1_2::chr1:15-20(+)\nATTTT\n'

    def create_test_files(self):
        self.fs.create_file('/test.bed', contents=self.bed_file_data)

    def test_get_slice_data_named_exon_success(self):
        # arrange
        expected = [
            ('chr1', 50, 260, 'exon1_1', '.', '+'),
            ('chr1', 90, 300, 'exon1_2', '.', '+')
        ]
        bed = pybedtools.BedTool(self.bed_file_data, from_string=True)
        params = {
            'flank_5': 50,
            'flank_3': 50,
            'length': 210,
            'offset': 40
        }

        # act
        actual = get_slice_data(bed, params)

        # assert
        self.assertEqual(actual, expected)

    def test_get_slice_data_unnamed_exon_success(self):
        # arrange
        expected = [
            ('chr1', 50, 260, 'region1_1', '.', '-'),
            ('chr1', 90, 300, 'region1_2', '.', '-')
        ]
        bed = pybedtools.BedTool('chr1\t100\t250\t.\t.\t-', from_string=True)
        params = {
            'flank_5': 50,
            'flank_3': 50,
            'length': 210,
            'offset': 40
        }

        # act
        actual = get_slice_data(bed, params)

        # assert
        self.assertEqual(actual, expected)

    def test_decrement_one_based_starts(self):
        # arrange
        input_file = [['1', '200', '300', 'name', '0', '+']]
        expected_row = ['1', '199', '300', 'name', '0', '+']

        # act
        result = decrement_one_based_starts(input_file, [])

        # assert
        self.assertEqual(expected_row, result[0])


        #def test_get_slices(self):
###
        # TODO: Rewrite test after Object Orientated refactor to allow for mocking.
        ###
        # expected_bed = (
        #     'chr1\t5\t10\tregion1_1\t.\t+\n'
        #     'chr1\t15\t20\tregion1_2\t.\t+\n'
        # )
        # expected_fasta = (
        #     '>region1_1::chr1:5-10(+)\n'
        #     'AGTCT\n'
        #     '>region1_2::chr1:15-20(+)\n'
        #     'ATTTT\n'
        # )
        # in_bed = StringIO('chr1\t5\t20\t.\t.\t+')
        # in_fasta = pybedtools.example_filename('test.fa')
        # params = {
        #     'input_bed': in_bed,
        #     'input_fasta': in_fasta,
        #     'flank_5': 0,
        #     'flank_3': 0,
        #     'length': 5,
        #     'offset': 10
        # }
        # slices = get_slices(params)
        # self.assertEqual(expected_bed, slices.head(as_string=True))
        # self.assertEqual(expected_fasta, slices.print_sequence())
        # expected_bed = (
        #     'chr1\t5\t10\texon1_1\t.\t-\n'
        #     'chr1\t15\t20\texon1_2\t.\t-\n'
        # )
        # expected_fasta = (
        #     '>exon1_1::chr1:5-10(-)\n'
        #     'AGACT\n'
        #     '>exon1_2::chr1:15-20(-)\n'
        #     'AAAAT\n'
        # )
        # params['input_bed'] = StringIO('chr1\t5\t20\texon1\t.\t-')
        # slices = get_slices(params)
        # self.assertEqual(expected_bed, slices.head(as_string=True))
        # self.assertEqual(expected_fasta, slices.print_sequence())


if __name__ == '__main__':
    unittest.main()
