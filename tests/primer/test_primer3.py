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

    def test_get_primer3_designs(self):
        config_test = {'PRIMER_TASK': 'generic', 'PRIMER_PICK_LEFT_PRIMER': 1, 'PRIMER_PICK_RIGHT_PRIMER': 1,
                       'PRIMER_OPT_SIZE': 20, 'PRIMER_MIN_SIZE': 18, 'PRIMER_MAX_SIZE': 30, 'P3_FILE_FLAG': 1,
                       'SEQUENCE_INCLUDED_REGION': [0, 212], 'PRIMER_EXPLAIN_FLAG': 1, 'PRIMER_MASK_TEMPLATE': 1,
                       'PRIMER_MASK_FAILURE_RATE': 1, 'PRIMER_MASK_5P_DIRECTION': 1, 'PRIMER_MASK_3P_DIRECTION': 0,
                       'PRIMER_MASK_KMERLIST_PATH': 'kmer/', 'PRIMER_WT_MASK_FAILURE_RATE': 1.0,
                       'PRIMER_NUM_RETURN': 20}
        slice_info = {'SEQUENCE_ID': 'mask_mask_1',
                      'SEQUENCE_TEMPLATE': 'CTTTTTTCTCTTTCCTTCTGCTTTTGTTTAAAGCGACAAGATGTTGCTCTTTTCCCAGGCTGGAATACAGTGGCATGATCATAGCTCAAGCTCCTGGGCTCAAGTGATCCTCCCGCCTCAGCCTCTCAAGTAGCTAGGACTACAGGCATATCACCACACCAGCGTTTTCTTTGTAGAGGCAGAGTCTCACTCTGTTGCTCAGGCAGGTGTTGAACTCCTGCCTCAAGCAATCCTCCCACCTCAGCCTCCCAGAGCCCTCAAATTATAAGCCACTGTGCTCGGGGCATCCTTTTTGGGGGGTAATCAGCAAACTGAAAAACCTCTTCTTACAACTCCCTATACATTCTCATTCCCAGTATAGAGGAGACTTTTTGTTTTTAAACACTTCCAAAGAATGCAAATTTATAATCCAGAGTATATACATTCTCACTGAATTATTGTACTGTTTCAG'}

        result = self.primer3_test_instance._get_primer3_designs(slice_info, config_test)

        expected = _get_file_content('tests/primer/primer3_output.json')

        self.assertDictEqual(result, json.loads(expected))


def _get_file_content(filename: str) -> str:
    with open(filename, 'r') as file:
        content = file.read()
    return content


if __name__ == '__main__':
    unittest.main()
