import json
import unittest
from unittest.mock import MagicMock, patch
from pyfakefs.fake_filesystem_unittest import TestCase

from primer.primer_pair import PrimerPair
from primer.primer3 import Primer3
from primer.slice_data import SliceData


class TestPrimer3(TestCase):

    def setUp(self):
        self._stringency_vector = [1, 0.5, 0.1]
        designer_config = {'stringency_vector': self._stringency_vector}
        p3_config = {'p3_config': 'p3_config'}
        self.primer3_test_instance = Primer3(designer_config, p3_config)

    @patch('primer.primer3.Primer3._get_primer_pairs')
    @patch.object(SliceData, 'parse_fasta')
    def test_get_primers(self, parse_fasta_file, mock_get_primers):
        slice1 = MagicMock(spec=SliceData)
        slice2 = MagicMock(spec=SliceData)
        parse_fasta_file.return_value = [slice1, slice2]

        primer_pair1 = MagicMock(spec=PrimerPair)
        primer_pair2 = MagicMock(spec=PrimerPair)
        primer_pair3 = MagicMock(spec=PrimerPair)
        pairs_slice1 = [primer_pair1, primer_pair2]
        pairs_slice2 = [primer_pair3]
        mock_get_primers.side_effect = [pairs_slice1, pairs_slice2]

        expected_primer_pairs = [primer_pair1, primer_pair2, primer_pair3]

        result_primer_pairs = self.primer3_test_instance.get_primers("FASTA_FILENAME")

        self.assertEqual(result_primer_pairs, expected_primer_pairs)

    @patch('primer.primer3.create_primer_pairs')
    @patch('primer.primer3.Primer3._get_primer3_designs')
    def test_get_primer_pairs(self, _get_primer3_designs, create_primer_pairs_mock):
        pair1 = MagicMock(spec=PrimerPair)
        pair2 = MagicMock(spec=PrimerPair)
        pair3 = MagicMock(spec=PrimerPair)
        create_primer_pairs_mock.side_effect = [[pair1], [pair2], [pair3]]

        slice = MagicMock(spec=SliceData)
        _get_primer3_designs.side_effect = [{"design1": "value"}, {"design2": "value"}, {"design2": "value"}]

        result_primer_pairs = self.primer3_test_instance._get_primer_pairs(slice)

        self.assertEqual(result_primer_pairs, [pair1, pair2, pair3])

    @patch('primer3.bindings.design_primers')
    def test_get_primer3_designs(self, design_primers):
        expected = {'designs': ["design1", "design2"]}
        design_primers.return_value = expected

        result = self.primer3_test_instance._get_primer3_designs({}, {})

        self.assertEqual(result, expected)


if __name__ == '__main__':
    unittest.main()
