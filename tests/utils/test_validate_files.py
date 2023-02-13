import unittest

from unittest.mock import patch
from pyfakefs.fake_filesystem_unittest import TestCase

from utils.exceptions import FileFormatError
from utils.validate_files import validate_bed_content, validate_bed_format, validate_fasta_format, validate_p3_csv, validate_score_tsv


class TestValidateFiles(TestCase):
    def setUp(self):
        self.setUpPyfakefs()
        self.bed_file_data = 'chr1\t100\t250\texon1\t.\t+'
        self.fasta_file_data = '>region1_1::chr1:5-10(+)\nAGTCT\n>region1_2::chr1:15-20(+)\nATTTT\n'

    def create_slicer_test_files(self):
        self.fs.create_file('/test.bed', contents=self.bed_file_data)
        self.fs.create_file('/test.fa', contents=self.fasta_file_data)

    def test_validate_bed_format_valid_bed_success(self):
        # arrange
        test_arg = '/test.bed'
        self.create_slicer_test_files()

        # act
        validate_bed_format(test_arg)

    def test_validate_bed_format_valid_long_bed_success(self):
        # arrange
        test_arg = '/longtest.bed'
        self.fs.create_file('/longtest.bed', contents='chr1\t100\t250\texon1\t.\t+\tfoo\tbar')

        # act
        validate_bed_format(test_arg)

    def test_validate_bed_format_invalid_csv_fail(self):
        # arrange
        test_arg = '/badtest.bed'
        self.fs.create_file('/badtest.bed', contents='chr1,100,250,exon1,.,+')
        expected = 'Unable to read in BED file correctly. Check file format on line 1.'

        # act
        with self.assertRaises(FileFormatError) as exception_context:
            validate_bed_format(test_arg)

        # assert
        self.assertEqual(str(exception_context.exception), expected)

    def test_validate_bed_format_invalid_len_fail(self):
        # arrange
        test_arg = '/badtest.bed'
        self.fs.create_file('/badtest.bed', contents='chr1\t100\t250\texon1\t.')
        expected = 'Unable to read in BED file correctly. Check file format on line 1.'

        # act
        with self.assertRaises(FileFormatError) as exception_context:
            validate_bed_format(test_arg)

        # assert
        self.assertEqual(str(exception_context.exception), expected)

    def test_validate_fasta_format_valid_fasta_success(self):
        # arrange
        test_arg = '/test.fa'
        self.create_slicer_test_files()

        # act
        validate_fasta_format(test_arg)

    def test_validate_fasta_format_invalid_bed_fail(self):
        # arrange
        test_arg = '/test.bed'
        self.create_slicer_test_files()
        expected = 'Unable to read in FastA file correctly. Check file format.'

        # act
        with self.assertRaises(FileFormatError) as exception_context:
            validate_fasta_format(test_arg)

        # assert
        self.assertEqual(str(exception_context.exception), expected)

    def test_validate_bed_content_valid_bed_success(self):
        # arrange
        test_arg = '/test.bed'
        self.create_slicer_test_files()

        # act
        validate_bed_content(test_arg)

    def test_validate_bed_content_valid_bed_alt_success(self):
        # arrange
        test_arg = '/goodtest.bed'
        self.fs.create_file('/goodtest.bed', contents='1\t100\t250\texon1\t.\t+')

        # act
        validate_bed_content(test_arg)

    def test_validate_bed_content_invalid_chr_fail(self):
        # arrange
        test_arg = '/badtest.bed'
        self.fs.create_file('/badtest.bed', contents='NOCHR\t100\t250\texon1\t.\t+')
        expected = 'Chromosome format incorrect on line 1: NOCHR'

        # act
        with self.assertRaises(ValueError) as exception_context:
            validate_bed_content(test_arg)

        # assert
        self.assertEqual(str(exception_context.exception), expected)

    def test_validate_bed_content_null_chr_fail(self):
        # arrange
        test_arg = '/badtest.bed'
        self.fs.create_file('/badtest.bed', contents='\t100\t250\texon1\t.\t+')
        expected = 'Chromosome format incorrect on line 1: '

        # act
        with self.assertRaises(ValueError) as exception_context:
            validate_bed_content(test_arg)

        # assert
        self.assertEqual(str(exception_context.exception), expected)

    def test_validate_bed_content_alpha_start_fail(self):
        # arrange
        test_arg = '/badtest.bed'
        self.fs.create_file('/badtest.bed', contents='chr1\ta\t250\texon1\t.\t+')
        expected = 'Start coordinate format incorrect on line 1: a'

        # act
        with self.assertRaises(ValueError) as exception_context:
            validate_bed_content(test_arg)

        # assert
        self.assertEqual(str(exception_context.exception), expected)

    def test_validate_bed_content_null_start_fail(self):
        # arrange
        test_arg = '/badtest.bed'
        self.fs.create_file('/badtest.bed', contents='chr1\t\t250\texon1\t.\t+')
        expected = 'Start coordinate format incorrect on line 1: '

        # act
        with self.assertRaises(ValueError) as exception_context:
            validate_bed_content(test_arg)

        # assert
        self.assertEqual(str(exception_context.exception), expected)

    def test_validate_bed_content_negative_start_fail(self):
        # arrange
        test_arg = '/badtest.bed'
        self.fs.create_file('/badtest.bed', contents='chr1\t-100\t250\texon1\t.\t+')
        expected = 'Start coordinate format incorrect on line 1: -100'

        # act
        with self.assertRaises(ValueError) as exception_context:
            validate_bed_content(test_arg)

        # assert
        self.assertEqual(str(exception_context.exception), expected)

    def test_validate_bed_content_alpha_end_fail(self):
        # arrange
        test_arg = '/badtest.bed'
        self.fs.create_file('/badtest.bed', contents='chr1\t100\ta\texon1\t.\t+')
        expected = 'End coordinate format incorrect on line 1: a'

        # act
        with self.assertRaises(ValueError) as exception_context:
            validate_bed_content(test_arg)

        # assert
        self.assertEqual(str(exception_context.exception), expected)

    def test_validate_bed_content_null_end_fail(self):
        # arrange
        test_arg = '/badtest.bed'
        self.fs.create_file('/badtest.bed', contents='chr1\t100\t\texon1\t.\t+')
        expected = 'End coordinate format incorrect on line 1: '

        # act
        with self.assertRaises(ValueError) as exception_context:
            validate_bed_content(test_arg)

        # assert
        self.assertEqual(str(exception_context.exception), expected)

    def test_validate_bed_content_negative_end_fail(self):
        # arrange
        test_arg = '/badtest.bed'
        self.fs.create_file('/badtest.bed', contents='chr1\t100\t-250\texon1\t.\t+')
        expected = 'End coordinate format incorrect on line 1: -250'

        # act
        with self.assertRaises(ValueError) as exception_context:
            validate_bed_content(test_arg)

        # assert
        self.assertEqual(str(exception_context.exception), expected)

    def test_validate_bed_content_larger_start_fail(self):
        # arrange
        test_arg = '/badtest.bed'
        self.fs.create_file('/badtest.bed', contents='chr1\t1000\t250\texon1\t.\t+')
        expected = 'End coordinate must be greater than start coordinate on line 1. Start: 1000 End: 250'

        # act
        with self.assertRaises(ValueError) as exception_context:
            validate_bed_content(test_arg)

        # assert
        self.assertEqual(str(exception_context.exception), expected)

    def test_validate_bed_content_big_diff_cood_fail(self):
        # arrange
        test_arg = '/badtest.bed'
        self.fs.create_file('/badtest.bed', contents='chr1\t100\t20100\texon1\t.\t+')
        expected = 'Difference between start coordinate and end coordinate must be less than 10000. On line 1 Difference: 20000'

        # act
        with self.assertRaises(ValueError) as exception_context:
            validate_bed_content(test_arg)

        # assert
        self.assertEqual(str(exception_context.exception), expected)

    def test_validate_bed_content_null_name_fail(self):
        # arrange
        test_arg = '/badtest.bed'
        self.fs.create_file('/badtest.bed', contents='chr1\t100\t250\t\t.\t+')
        expected = 'Error with name field, if no name is supplied please mark with a \'.\' on line 1: '

        # act
        with self.assertRaises(ValueError) as exception_context:
            validate_bed_content(test_arg)

        # assert
        self.assertEqual(str(exception_context.exception), expected)

    def test_validate_bed_content_null_score_fail(self):
        # arrange
        test_arg = '/badtest.bed'
        self.fs.create_file('/badtest.bed', contents='chr1\t100\t250\texon1\t\t+')
        expected = 'Error with score field, if no score is supplied please mark with a \'.\' on line 1: '

        # act
        with self.assertRaises(ValueError) as exception_context:
            validate_bed_content(test_arg)

        # assert
        self.assertEqual(str(exception_context.exception), expected)

    def test_validate_bed_content_invalid_strand_fail(self):
        # arrange
        test_arg = '/badtest.bed'
        self.fs.create_file('/badtest.bed', contents='chr1\t100\t250\texon1\t.\ta')
        expected = 'Strand format incorrect on line 1: a'

        # act
        with self.assertRaises(ValueError) as exception_context:
            validate_bed_content(test_arg)

        # assert
        self.assertEqual(str(exception_context.exception), expected)

    def test_validate_bed_content_both_strand_fail(self):
        # arrange
        test_arg = '/badtest.bed'
        self.fs.create_file('/badtest.bed', contents='chr1\t100\t250\texon1\t.\t-+')
        expected = 'Strand format incorrect on line 1: -+'

        # act
        with self.assertRaises(ValueError) as exception_context:
            validate_bed_content(test_arg)

        # assert
        self.assertEqual(str(exception_context.exception), expected)

    def test_validate_p3_csv_wrong_headers_fail(self):
        # arrange
        test_arg = '/not_p3_output.csv'
        self.fs.create_file('/not_p3_output.csv', contents='Incorrect,Headers')
        expected = 'Unexpected columns in Primer3 CSV'

        # act
        with self.assertRaises(FileFormatError) as exception_context:
            validate_p3_csv(test_arg)

        # assert
        self.assertIn(expected, str(exception_context.exception))

    def test_validate_p3_csv_correct_headers_success(self):
        # arrange
        test_arg = '/p3_output.csv'
        headers = [
            'primer', 'sequence', 'chr', 'primer_start', 'primer_end',
            'tm', 'gc_percent', 'penalty', 'self_any_th', 'self_end_th',
            'hairpin_th', 'end_stability'
        ]
        self.fs.create_file('/p3_output.csv', contents=','.join(headers))

        # act
        validate_p3_csv(test_arg)

    def test_validate_score_tsv_wrong_headers_fail(self):
        # arrange
        test_arg = '/not_score_output.tsv'
        self.fs.create_file('/not_score_output.tsv', contents='Incorrect\tHeaders')
        expected = 'Unexpected columns in Scoring TSV'

        # act
        with self.assertRaises(FileFormatError) as exception_context:
            validate_score_tsv(test_arg)

        # assert
        self.assertIn(expected, str(exception_context.exception))

    def test_validate_score_tsv_correct_headers_success(self):
        # arrange
        test_arg = '/scoring_output.tsv'
        headers = [
            'Targeton', 'Primer pair', 'A/B/Total', '0', '1', '2',
            '3', '4', '5', '6', '7', '8', '9', '10', 'WGE format', 'Score'
        ]
        self.fs.create_file('/scoring_output.tsv', contents='\t'.join(headers))

        # act
        validate_score_tsv(test_arg)


if __name__ == '__main__':
    unittest.main()
