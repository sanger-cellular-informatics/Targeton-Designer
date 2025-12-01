import logging
from unittest.mock import patch

from pyfakefs.fake_filesystem_unittest import TestCase
from tests.utils.utils import CapturingStreamHandler

from primer.slice_data import SliceData


class TestSliceData(TestCase):

    def setUp(self):
        self.setUpPyfakefs()
        # Create a custom stream handler to capture logs
        self.handler = CapturingStreamHandler()
        self.logger = self.handler.get_logger(self.handler)

    def tearDown(self):
        # Remove the handler after each test to reset logging
        logger = logging.getLogger()
        logger.removeHandler(self.handler)

    @patch('primer.slice_data.get_seq_from_ensembl_by_coords')
    def test_p3_input(self, mock_get_seq):

        slice_sample = SliceData(name = 'slice_name',
                                 start = 100,
                                 end = 110,
                                 strand = 'strand',
                                 chromosome = 'chromosome',
                                 bases = 'slice_bases',
                                 flanking_region = 5,
                                 exclusion_region = 2)

        expected_p3_input = {'SEQUENCE_ID': 'slice_name', 'SEQUENCE_TEMPLATE': 'slice_bases',
                             'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST': [0, 3, 9, 2]}

        result = slice_sample.p3_input

        self.assertEqual(result, expected_p3_input)

    def test_p3_input_when_flanking_zero(self):

            slice_sample = SliceData(name = 'slice_name',
                                    start = 100,
                                    end = 110,
                                    strand = 'strand',
                                    chromosome = 'chromosome',
                                    bases = 'slice_bases',
                                    flanking_region = 0,
                                    exclusion_region = 200)

            expected_p3_input = {'SEQUENCE_ID': 'slice_name', 'SEQUENCE_TEMPLATE': 'slice_bases',
                                'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST': []}

            result = slice_sample.p3_input

            self.assertEqual(result, expected_p3_input)

    #when flanking_region is = 0
    def test_get_first_slice_no_flank_calls(self):
        slices_fasta_file = 'one_slice.fa'
        self.fs.create_file(slices_fasta_file, contents='>region1_1::chr1:5-10(+)\nGTGATCGAGGAGTTCTA')

        expected = SliceData(name = 'region1_1',
                             start = 5,
                             end = 10,
                             strand = '+',
                             chromosome = '1',
                             bases='GTGATCGAGGAGTTCTA',
                             flanking_region = 0,
                             exclusion_region = 0)
        result = SliceData.get_first_slice_data(slices_fasta_file, flanking = 0, exclusion_region = 0)

        self.assertEqual(result.name, expected.name)
        self.assertEqual(result.start, expected.start)
        self.assertEqual(result.end, expected.end)
        self.assertEqual(result.strand, expected.strand)
        self.assertEqual(result.bases, expected.bases)
        self.assertEqual(result.chromosome, expected.chromosome)

    #when flanking_region is > 0
    @patch('primer.flanking.get_seq_from_ensembl_by_coords')
    def test_get_first_slice_data_with_flanking_calls(self, mock_get_seq):

        slices_fasta_file = 'one_slice_no_flank.fa'
        self.fs.create_file(
            slices_fasta_file,
            contents='>region1_1::chr1:5-10(+)\nGTGATCGAGGAGTTCTA'
        )

        flanking = 2
        exclusion_region = 7

        expected_start = 3
        expected_end = 12

        fake_extended_seq = 'ACTGACTG'
        mock_get_seq.return_value = fake_extended_seq

        result = SliceData.get_first_slice_data(
            slices_fasta_file,
            flanking=flanking,
            exclusion_region=exclusion_region,
        )

        mock_get_seq.assert_called_once_with(
            chromosome='1',
            start=expected_start,
            end=expected_end,
            strand='+',
        )

        self.assertIsInstance(result, SliceData)
        self.assertEqual(result.name, 'region1_1')
        self.assertEqual(result.chromosome, '1')
        self.assertEqual(result.strand, '+')
        self.assertEqual(result.start, expected_start)
        self.assertEqual(result.end, expected_end)
        self.assertEqual(result.bases, fake_extended_seq)
        self.assertEqual(result.flanking_region, flanking)
        self.assertEqual(result.exclusion_region, exclusion_region)

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
                             flanking_region = 0,
                             exclusion_region = 0)

        result = SliceData.get_first_slice_data(slices_fasta_file, flanking = 0, exclusion_region = 0)

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
            SliceData.get_first_slice_data(wrong_fasta_file, flanking = 0, exclusion_region = 0)

        self.assertEqual(str(error.exception), f"Unable to parse the FASTA file '{wrong_fasta_file}'")

    def test_get_first_slice_when_empty_fasta_file(self):
        empty_fasta = "empty.fa"
        self.fs.create_file(empty_fasta, contents='')

        with self.assertRaises(ValueError) as error:
            SliceData.get_first_slice_data(empty_fasta, flanking = 0, exclusion_region = 0)

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
                             flanking_region = 0,
                             exclusion_region = 0)

        result_1 = SliceData.get_first_slice_data(mocked_fasta_1, flanking = 0, exclusion_region = 0)

        # check with chr parsed correctly
        self.assertEqual(result_1, expected)

        result_2 = SliceData.get_first_slice_data(mocked_fasta_2, flanking = 0, exclusion_region = 0)

        # check with ch parsed correctly
        self.assertEqual(result_2, expected)

        # check with numerical chromosome
        result_3 = SliceData.get_first_slice_data(mocked_fasta_3, flanking = 0, exclusion_region = 0)

        self.assertEqual(result_3, expected)

    def test_fasta_file_parsing_chromosome_with_invalid_characters(self):
        mocked_fasta = 'mocked_fasta.fa'
        self.fs.create_file(mocked_fasta, contents='>region1_1::xyz$#r1:5-10(+)\nGTGATCGAGGAGTTCTA')

        with self.assertRaises(ValueError) as ex:
            _ = SliceData.get_first_slice_data(mocked_fasta, flanking = 0, exclusion_region = 0)

        self.assertTrue("does not match the expected format" in str(ex.exception))

    @patch("primer.flanking.logger")
    @patch("primer.flanking.get_seq_from_ensembl_by_coords")
    @patch("primer.flanking._get_chromosome_length")
    def test_fasta_file_flanking_exceeds_chrom_end(
        self,
        mock_get_chr_len,
        mock_get_seq,
        mock_logger,
    ):
        """
        When flanking pushes the extended region beyond the chromosome end,
        get_first_slice_data should clamp extended_end to chrom length and log an warning.
        """

        mock_get_chr_len.return_value = 1000

        fasta_path = "one_slice_exceed.fa"
        self.fs.create_file(
            fasta_path,
            contents=">region1_1::chr1:900-950(+)\n"
                     "ACGTACGTACGTACGTACGT\n"
        )

        fake_extended_seq = "N" * (1000 - 800 + 1)
        mock_get_seq.return_value = fake_extended_seq

        result = SliceData.get_first_slice_data(
            fasta=fasta_path,
            flanking=100,
            exclusion_region=0,
        )

        mock_get_chr_len.assert_called_once_with("1")

        mock_get_seq.assert_called_once_with(
            chromosome="1",
            start=800,
            end=1000,
            strand="+",
        )

        self.assertEqual(result.chromosome, "1")
        self.assertEqual(result.start, 800)
        self.assertEqual(result.end, 1000)
        self.assertEqual(result.bases, fake_extended_seq)

        warning_calls = [c for c in mock_logger.warning.call_args_list]
        self.assertTrue(
            any("Flanking region expands beyond chromosome 1 end" in str(c) for c in warning_calls),
            "Expected an warning log about flanking beyond chromosome end"
        )


class TestGetSliceFromRegion(TestCase):

    @patch('primer.flanking.get_seq_from_ensembl_by_coords')
    def test_get_slice_from_region_valid(self, mock_get_seq):
        """Test successful parsing of region string and creation of SliceData."""
        mock_get_seq.return_value = 'ACTGACTG'

        expected = SliceData(name='ABCD',
                             start=54050,
                             end=54250,
                             strand='+',
                             chromosome='19',
                             bases='ACTGACTG',
                             flanking_region=50,
                             exclusion_region=20)

        result = SliceData.get_slice_from_region(
            targeton_id='ABCD',
            region='chr19:54100-54200',
            strand='+',
            flanking=50,
            exclusion_region=20
        )

        mock_get_seq.assert_called_once_with(
            chromosome='19',
            start=54050,
            end=54250,
            strand='+',
        )

        self.assertIsInstance(result, SliceData)
        self.assertEqual(result.name, expected.name)
        self.assertEqual(result.start, expected.start)
        self.assertEqual(result.end, expected.end)
        self.assertEqual(result.strand, expected.strand)
        self.assertEqual(result.chromosome, expected.chromosome)
        self.assertEqual(result.bases, expected.bases)
        self.assertEqual(result.flanking_region, expected.flanking_region)
        self.assertEqual(result.exclusion_region, expected.exclusion_region)

    @patch('primer.flanking.get_seq_from_ensembl_by_coords')
    def test_get_slice_from_region_valid_neg(self, mock_get_seq):
        """Test that strand '-' is set"""
        mock_get_seq.return_value = 'ACTGACTG'

        expected = SliceData(name='ABCD',
                             start=54050,
                             end=54250,
                             strand='-',
                             chromosome='19',
                             bases='ACTGACTG',
                             flanking_region=50,
                             exclusion_region=20)

        result = SliceData.get_slice_from_region(
            targeton_id='ABCD',
            region='chr19:54100-54200',
            strand='-',
            flanking=50,
            exclusion_region=20
        )

        mock_get_seq.assert_called_once_with(
            chromosome='19',
            start=54050,
            end=54250,
            strand='-',
        )

        self.assertEqual(result.strand, expected.strand)
        self.assertEqual(result.bases, expected.bases)

    def test_get_slice_from_region_invalid_region(self):
        """Test that incorrectly formatted region raises ValueError."""
        with self.assertRaises(ValueError) as ex:
            SliceData.get_slice_from_region(
                targeton_id='ABCD',
                region='chr19_54100-54200',  # invalid format
                strand='+',
                flanking=50,
                exclusion_region=20
            )

        self.assertTrue("does not match the expected format" in str(ex.exception))

    @patch("primer.flanking.logger")
    @patch("primer.flanking.get_seq_from_ensembl_by_coords")
    @patch("primer.flanking._get_chromosome_length")
    def test_get_slice_from_region_flanking_exceeds_chrom_end(
        self,
        mock_get_chr_len,
        mock_get_seq,
        mock_logger,
    ):

        mock_get_chr_len.return_value = 1000

        targeton_id = "REG1"
        region = "chr1:900-950"
        strand = "+"
        flanking = 100
        exclusion_region = 0

        fake_extended_seq = "N" * (1000 - 800 + 1)
        mock_get_seq.return_value = fake_extended_seq

        result = SliceData.get_slice_from_region(
            targeton_id=targeton_id,
            region=region,
            strand=strand,
            flanking=flanking,
            exclusion_region=exclusion_region,
        )

        mock_get_chr_len.assert_called_once_with("1")

        mock_get_seq.assert_called_once_with(
            chromosome="1",
            start=800,
            end=1000,
            strand="+",
        )

        self.assertEqual(result.name, targeton_id)
        self.assertEqual(result.chromosome, "1")
        self.assertEqual(result.start, 800)
        self.assertEqual(result.end, 1000)
        self.assertEqual(result.bases, fake_extended_seq)

        warning_calls = [c for c in mock_logger.warning.call_args_list]
        self.assertTrue(
            any("Flanking region expands beyond chromosome 1 end" in str(c) for c in warning_calls),
            "Expected a warning log about flanking beyond chromosome end"
        )
