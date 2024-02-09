from unittest.mock import patch

from pyfakefs.fake_filesystem_unittest import TestCase

from primer.slice_data import SliceData


class TestSliceData(TestCase):

    def setUp(self):
        self.setUpPyfakefs()

    def test_p3_input(self):
        slice_sample = SliceData('slice_name', 'start', 'end', 'strand', 'chrom', 'slice_bases')
        expected_p3_input = {'SEQUENCE_ID': 'slice_name', 'SEQUENCE_TEMPLATE': 'slice_bases'}

        result = slice_sample.p3_input

        self.assertEqual(result, expected_p3_input)

    def test_read_input_fasta_valid_success(self):
        slices_fasta_file = '/fasta.fa'
        self.fs.create_file(slices_fasta_file, contents='>region1_1::chr1:5-10(+)\nGTGATCGAGGAGTTCTAC')

        expected = [SliceData('region1_1', '5', '10', '+', 'chr1', 'GTGATCGAGGAGTTCTAC')]

        result = SliceData.parse_fasta(slices_fasta_file)

        self.assertEqual(result[0].name, expected[0].name)
        self.assertEqual(result[0].start, expected[0].start)
        self.assertEqual(result[0].end, expected[0].end)
        self.assertEqual(result[0].strand, expected[0].strand)
        self.assertEqual(result[0].bases, expected[0].bases)
