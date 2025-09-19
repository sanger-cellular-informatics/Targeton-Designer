from unittest.mock import patch

from pyfakefs.fake_filesystem_unittest import TestCase

from primer.slice_data import SliceData


class TestSliceData(TestCase):

    def setUp(self):
        self.setUpPyfakefs()

    @patch('primer.slice_data.get_seq_from_ensembl_by_coords')
    def test_p3_input(self, mock_get_seq):

        slice_sample = SliceData(name = 'slice_name',
                                 start = 100,
                                 end = 110,
                                 strand = 'strand',
                                 chromosome = 'chromosome',
                                 bases = 'slice_bases',
                                 region_padding = 5,
                                 region_avoid = 2)

        expected_p3_input = {'SEQUENCE_ID': 'slice_name', 'SEQUENCE_TEMPLATE': 'slice_bases',
                             'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST': [0, 3, 9, 2]}

        result = slice_sample.p3_input

        self.assertEqual(result, expected_p3_input)

    def test_get_first_slice(self):
        slices_fasta_file = 'one_slice.fa'
        self.fs.create_file(slices_fasta_file, contents='>region1_1::chr1:5-10(+)\nGTGATCGAGGAGTTCTA')

        expected = SliceData(name = 'region1_1',
                             start = 5,
                             end = 10,
                             strand = '+',
                             chromosome = '1',
                             bases='GTGATCGAGGAGTTCTA',
                             region_padding = 0,
                             region_avoid = 0)
        result = SliceData.get_first_slice_data(slices_fasta_file, padding = 0, region_avoid = 0)

        self.assertEqual(result.name, expected.name)
        self.assertEqual(result.start, expected.start)
        self.assertEqual(result.end, expected.end)
        self.assertEqual(result.strand, expected.strand)
        self.assertEqual(result.bases, expected.bases)
        self.assertEqual(result.chromosome, expected.chromosome)

    @patch('custom_logger.custom_logger.CustomLogger.warning')
    def test_get_first_slice_when_more_than_one_slice(self, logger_warning):
        slices_fasta_file = 'two_slices.fa'
        self.fs.create_file(slices_fasta_file,
                            contents='>region1_1::chr1:5-10(+)\nGTGATCGAGGAGTTCTA\n'
                                     '>region2_1::chr1:5-10(+)\nAAAAGGGCCCTTTAAAA')

        expected = SliceData(name = 'region1_1',
                             start = 5,
                             end = 10,
                             strand = '+',
                             chromosome = '1',
                             bases = 'GTGATCGAGGAGTTCTA',
                             region_padding = 0,
                             region_avoid = 0)

        result = SliceData.get_first_slice_data(slices_fasta_file, padding = 0, region_avoid = 0)

        self.assertEqual(result.name, expected.name)
        self.assertEqual(result.start, expected.start)
        self.assertEqual(result.end, expected.end)
        self.assertEqual(result.strand, expected.strand)
        self.assertEqual(result.bases, expected.bases)
        self.assertEqual(result.chromosome, expected.chromosome)
        logger_warning.assert_called_once_with(
            f"The FASTA file '{slices_fasta_file}' contains more than one pre-targeton. "
            "Only the first pre-targeton is taken.")

    def test_get_first_slice_when_wrong_sequence_format(self):
        wrong_fasta_file = "wrong.fa"
        self.fs.create_file(wrong_fasta_file, contents='WRONG_SEQUENCE_PATTERN\nGTGATCGAGGAGTTCTA')

        with self.assertRaises(ValueError) as error:
            SliceData.get_first_slice_data(wrong_fasta_file, padding = 0, region_avoid = 0)

        self.assertEqual(str(error.exception), f"Unable to parse the FASTA file '{wrong_fasta_file}'")

    def test_get_first_slice_when_empty_fasta_file(self):
        empty_fasta = "empty.fa"
        self.fs.create_file(empty_fasta, contents='')

        with self.assertRaises(ValueError) as error:
            SliceData.get_first_slice_data(empty_fasta, padding = 0, region_avoid = 0)

        self.assertEqual(str(error.exception), f"Unable to parse the FASTA file '{empty_fasta}'")

    def test_fasta_file_parsing_chromosome_with_characters(self):
        mocked_fasta_1 = 'mocked_fasta_1.fa'
        mocked_fasta_2 = 'mocked_fasta_2.fa'
        mocked_fasta_3 = 'mocked_fasta_3.fa'

        # create fasta with ch
        self.fs.create_file(mocked_fasta_1, contents='>region1_1::ch1:5-10(+)\nGTGATCGAGGAGTTCTA')

        # create fasta with chr
        self.fs.create_file(mocked_fasta_2, contents='>region1_1::chr1:5-10(+)\nGTGATCGAGGAGTTCTA')

        # create fasta with only chromosome number
        self.fs.create_file(mocked_fasta_3, contents='>region1_1::1:5-10(+)\nGTGATCGAGGAGTTCTA')

        expected = SliceData(name = 'region1_1',
                             start = 5,
                             end = 10,
                             strand = '+',
                             chromosome = '1',
                             bases = 'GTGATCGAGGAGTTCTA',
                             region_padding = 0,
                             region_avoid = 0)

        result_1 = SliceData.get_first_slice_data(mocked_fasta_1, padding = 0, region_avoid = 0)

        # check with chr parsed correctly
        self.assertEqual(result_1, expected)

        result_2 = SliceData.get_first_slice_data(mocked_fasta_2, padding = 0, region_avoid = 0)

        # check with ch parsed correctly
        self.assertEqual(result_2, expected)

        # check with numerical chromosome
        result_3 = SliceData.get_first_slice_data(mocked_fasta_3, padding = 0, region_avoid = 0)

        self.assertEqual(result_3, expected)

    def test_fasta_file_parsing_chromosome_with_invalid_characters(self):
        mocked_fasta = 'mocked_fasta.fa'
        self.fs.create_file(mocked_fasta, contents='>region1_1::xyz$#r1:5-10(+)\nGTGATCGAGGAGTTCTA')

        with self.assertRaises(ValueError) as ex:
            _ = SliceData.get_first_slice_data(mocked_fasta, padding = 0, region_avoid = 0)

        self.assertTrue("does not match the expected format" in str(ex.exception))
