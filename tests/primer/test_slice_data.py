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

    def test_get_first_slice(self):
        slices_fasta_file = 'one_slice.fa'
        self.fs.create_file(slices_fasta_file, contents='>region1_1::chr1:5-10(+)\nGTGATCGAGGAGTTCTA')

        expected = SliceData(name='region1_1', start='5', end='10', strand='+', chrom='chr1', bases='GTGATCGAGGAGTTCTA')

        result = SliceData.get_first_pre_targeton(slices_fasta_file)

        self.assertEqual(result.name, expected.name)
        self.assertEqual(result.start, expected.start)
        self.assertEqual(result.end, expected.end)
        self.assertEqual(result.strand, expected.strand)
        self.assertEqual(result.bases, expected.bases)

    @patch('custom_logger.custom_logger.CustomLogger.warning')
    def test_get_first_slice_when_more_than_one_slice(self, logger_warning):
        slices_fasta_file = 'two_slices.fa'
        self.fs.create_file(slices_fasta_file,
                            contents='>region1_1::chr1:5-10(+)\nGTGATCGAGGAGTTCTA\n'
                                     '>region2_1::chr1:5-10(+)\nAAAAGGGCCCTTTAAAA')

        expected = SliceData(name='region1_1', start='5', end='10', strand='+', chrom='chr1', bases='GTGATCGAGGAGTTCTA')

        result = SliceData.get_first_pre_targeton(slices_fasta_file)

        self.assertEqual(result.name, expected.name)
        self.assertEqual(result.start, expected.start)
        self.assertEqual(result.end, expected.end)
        self.assertEqual(result.strand, expected.strand)
        self.assertEqual(result.bases, expected.bases)
        logger_warning.assert_called_once_with("The FASTA file 'two_slices.fa' contains more than one pre-targeton. "
                                               "Only the first pre-targeton is taken.")
