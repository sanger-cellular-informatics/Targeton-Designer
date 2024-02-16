from unittest import TestCase
from unittest.mock import patch, MagicMock

import requests

from primer.ensembl import get_seq_from_ensembl_by_coords


class TestEnsemble(TestCase):
    @patch('requests.get')
    @patch('time.sleep', return_value=None)
    def test_get_seq_from_ensembl_by_coords_success(self, mock_sleep, mock_get):
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.text = "Mocked response text"
        mock_get.return_value = mock_response

        result = get_seq_from_ensembl_by_coords("X", 1000, 2000)

        mock_get.assert_called_once_with('https://rest.ensembl.org/sequence/region/human/X:1000..2000:1',
                                         headers={'Content-type': 'text/plain'})

        self.assertEqual(result, "Mocked response text")

    @patch('requests.get')
    @patch('time.sleep', return_value=None)
    def test_get_seq_from_ensembl_by_coords_failure(self, mock_sleep, mock_get):
        mock_response = MagicMock()
        mock_response.status_code = 404
        mock_get.return_value = mock_response

        with self.assertRaises(requests.exceptions.RequestException):
            get_seq_from_ensembl_by_coords("X", 1000, 2000)

        mock_get.assert_called_once_with('https://rest.ensembl.org/sequence/region/human/X:1000..2000:1',
                                         headers={'Content-type': 'text/plain'})
