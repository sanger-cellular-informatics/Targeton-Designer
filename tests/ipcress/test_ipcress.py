import sys
import unittest
from unittest import TestCase
from unittest.mock import patch

from ipcress.ipcress import Ipcress
from utils.arguments_parser import ParsedInputArguments


class TestIpcress(TestCase):
    def setUp(self):
        with patch.object(sys, 'argv', ["./designer.sh", "ipcress"]):
            parsed_input = ParsedInputArguments()
            args = parsed_input.get_args()
            self.ipcress = Ipcress(args)

    def test_prettify_output_defined(self):
        expected = 'test_cmd --pretty true'
        prettify_param = True
        cmd = 'test_cmd'

        actual = self.ipcress.prettify_output(prettify_param, cmd)

        self.assertEqual(actual, expected)

    def test_prettify_output_undef(self):
        expected = 'test_cmd --pretty false'
        prettify_param = ''
        cmd = 'test_cmd'

        actual = self.ipcress.prettify_output(prettify_param, cmd)

        self.assertEqual(actual, expected)

    @patch('utils.logger.Logger.log')
    def test_validate_primers_pretty_prints_skipping_message(self, mock_print):
        ipcress_output = ''
        primer_data = {}
        params = {"pretty": True, "quiet": False}

        self.ipcress.validate_primers(ipcress_output, primer_data, params)

        mock_print.assert_called_with('Output is pretty, skipping validation')

    @patch('utils.logger.Logger.log')
    def test_validate_primers_wrong_coord_prints_warning(self, mock_print):
        ipcress_output = 'ipcress: 1:filter(unmasked) test-primer-pair 200 A 122 0 B 456 0 forward\n'
        primer_data = {'test-primer-pair': {'F': {'start': 123, 'seq': 'ATCG'}, 'R': {'start': 456, 'seq': 'GCTA'}}}
        params = {"pretty": False, "quiet": False}

        self.ipcress.validate_primers(ipcress_output, primer_data, params, validate_coords=True)

        mock_print.assert_called_with('No valid primer pair found for test-primer-pair')

    @patch('utils.logger.Logger.log')
    def test_validate_primers_mismatch_prints_warning(self, mock_print):
        ipcress_output = 'ipcress: 1:filter(unmasked) test-primer-pair 200 A 123 0 B 456 1 forward\n'
        primer_data = {'test-primer-pair': {'F': {'start': 123, 'seq': 'ATCG'}, 'R': {'start': 456, 'seq': 'GCTA'}}}
        params = {"pretty": False, "quiet": False}

        self.ipcress.validate_primers(ipcress_output, primer_data, params)

        mock_print.assert_called_with('No valid primer pair found for test-primer-pair')

    @patch('utils.logger.Logger.log')
    def test_validate_primers_not_forward_prints_warning(self, mock_print):
        ipcress_output = 'ipcress: 1:filter(unmasked) test-primer-pair 200 A 123 0 B 456 0 revcomp\n'
        primer_data = {'test-primer-pair': {'F': {'start': 123, 'seq': 'ATCG'}, 'R': {'start': 456, 'seq': 'GCTA'}}}
        params = {"pretty": False, "quiet": False}

        self.ipcress.validate_primers(ipcress_output, primer_data, params)

        mock_print.assert_called_with('No valid primer pair found for test-primer-pair')

    @patch('utils.logger.Logger.log')
    def test_validate_primers_match_does_not_print_warning(self, mock_print):
        ipcress_output = 'ipcress: 1:filter(unmasked) test-primer-pair 200 A 123 0 B 456 0 forward\n'
        primer_data = {'test-primer-pair': {'F': {'start': 123, 'seq': 'ATCG'}, 'R': {'start': 456, 'seq': 'GCTA'}}}
        params = {"pretty": False, "quiet": False}

        self.ipcress.validate_primers(ipcress_output, primer_data, params)

        mock_print.assert_called_once_with('Validating primers...')


if __name__ == '__main__':
    unittest.main()
