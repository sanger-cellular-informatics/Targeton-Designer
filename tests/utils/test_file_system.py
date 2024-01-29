import unittest
from unittest.mock import patch

from pyfakefs.fake_filesystem_unittest import TestCase

from src.utils.file_system import check_file_exists


class TestFileSystem(TestCase):
    def setUp(self):
        self.setUpPyfakefs()
        self.bed_file_data = 'chr1\t100\t250\texon1\t.\t+'
        self.fasta_file_data = '>region1_1::chr1:5-10(+)\nAGTCT\n>region1_2::chr1:15-20(+)\nATTTT\n'

    def create_slicer_test_files(self):
        self.fs.create_file('/test.bed', contents=self.bed_file_data)
        self.fs.create_file('/test.fa', contents=self.fasta_file_data)

    def test_check_file_exists_valid_bed_arg_success(self):
        # arrange
        test_arg = '/test.bed'
        self.create_slicer_test_files()

        # act
        check_file_exists(test_arg)

    def test_check_file_exists_valid_fasta_arg_success(self):
        # arrange
        test_arg = '/test.fa'
        self.create_slicer_test_files()

        # act
        check_file_exists(test_arg)

    def test_check_file_exists_invalid_args_fail(self):
        # arrange
        test_arg = '/test.bed'
        expected = "Unable to find file: /test.bed"

        # act
        with self.assertRaises(FileNotFoundError) as exception_context:
            check_file_exists(test_arg)

        # assert
        self.assertEqual(str(exception_context.exception), expected)


if __name__ == '__main__':
    unittest.main()
