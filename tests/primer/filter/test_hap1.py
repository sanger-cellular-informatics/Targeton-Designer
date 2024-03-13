from unittest import TestCase
from unittest.mock import patch, Mock

import requests

from primer.filter.hap1 import contain_variant


class TestPrimer3(TestCase):

    @patch('requests.get')
    def test_contain_variant_with_variants(self, mock_get):
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.text = '{"variants": ["variant1", "variant2"]}'
        mock_get.return_value = mock_response

        chromosome = "chr1"
        start = 100
        end = 200

        self.assertTrue(contain_variant(chromosome, start, end))

    @patch('requests.get')
    def test_contain_variant_without_variants(self, mock_get):
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.text = '{"variants": []}'
        mock_get.return_value = mock_response

        chromosome = "chr1"
        start = 100
        end = 200

        self.assertFalse(contain_variant(chromosome, start, end))

    @patch('requests.get')
    def test_contain_variant_request_exception(self, mock_get):
        mock_get.side_effect = requests.exceptions.RequestException

        chromosome = "chr1"
        start = 100
        end = 200

        with self.assertRaises(requests.exceptions.RequestException):
            contain_variant(chromosome, start, end)
