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
    def test_get_primer_pairs_for_all_slices_in_fasta_file(self, parse_fasta_file, mock_get_primer_pairs):
        slice1 = MagicMock(spec=SliceData)
        slice2 = MagicMock(spec=SliceData)
        parse_fasta_file.return_value = [slice1, slice2]

        primer_pair1 = MagicMock(spec=PrimerPair)
        primer_pair2 = MagicMock(spec=PrimerPair)
        primer_pair3 = MagicMock(spec=PrimerPair)
        pairs_slice1 = [primer_pair1, primer_pair2]
        pairs_slice2 = [primer_pair3]
        mock_get_primer_pairs.side_effect = [pairs_slice1, pairs_slice2]

        expected_primer_pairs = [primer_pair1, primer_pair2, primer_pair3]

        result_primer_pairs = self.primer3_test_instance.get_primers("FASTA_FILENAME")

        self.assertEqual(result_primer_pairs, expected_primer_pairs)

    @patch('primer.primer3.build_primer_pairs')
    @patch('primer.primer3.Primer3._get_primer3_designs')
    def test_get_slice_primer_pairs_through_all_stringencies(self, _get_primer3_designs, build_primer_pairs_mock):
        stringency_vector = [1, 2, 3]
        self.primer3_test_instance._stringency_vector = stringency_vector
        
        pair1 = MagicMock(spec=PrimerPair)
        pair2 = MagicMock(spec=PrimerPair)
        pairs_stringency_1 = [pair1, pair2]
        
        pair3 = MagicMock(spec=PrimerPair)
        pairs_stringency_2 = [pair3]
        
        pair4 = MagicMock(spec=PrimerPair)
        pairs_stringency_3 = [pair4]
        
        build_primer_pairs_mock.side_effect = [pairs_stringency_1, pairs_stringency_2, pairs_stringency_3]

        slice = MagicMock(spec=SliceData)
        _get_primer3_designs.side_effect = \
            [{"designs_stringency_1": "A"}, {"designs_stringency_2": "B"}, {"designs_stringency_3": "C"}]

        result_primer_pairs = self.primer3_test_instance._get_primer_pairs(slice)

        self.assertEqual(result_primer_pairs, [pair1, pair2, pair3, pair4])

    @patch('primer3.bindings.design_primers')
    def test_get_primer3_designs(self, design_primers):
        expected = {"designs": ["design1", "design2"]}
        design_primers.return_value = expected

        result = self.primer3_test_instance._get_primer3_designs({}, {})

        self.assertEqual(result, expected)


if __name__ == '__main__':
    unittest.main()
