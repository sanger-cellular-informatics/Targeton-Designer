from unittest import TestCase
from unittest.mock import patch, Mock

import requests

from utils.get_data.hap1 import contain_variant


class TestHap1(TestCase):

    @patch('requests.get')
    def test_contain_variant_with_variants(self, request_get):
        http_response = Mock(status_code=200, text='{"variants": ["variant1", "variant2"]}')

        request_get.return_value = http_response

        self.assertTrue(contain_variant(chromosome="1", start=100, end=200))

    @patch('requests.get')
    def test_contain_variant_without_variants(self, request_get):
        http_response = Mock(status_code=200, text='{"variants": []}')

        request_get.return_value = http_response

        self.assertFalse(contain_variant(chromosome="1", start=100, end=200))

    @patch('requests.get')
    def test_contain_variant_request_exception(self, mock_get):
        mock_get.side_effect = requests.exceptions.RequestException

        with self.assertRaises(requests.exceptions.RequestException):
            contain_variant(chromosome="1", start=100, end=200)
