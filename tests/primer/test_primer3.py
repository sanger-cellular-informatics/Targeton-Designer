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


class IntegrationTestPrimer3(TestCase):

    def setUp(self):
        self.setUpPyfakefs()

    def test_get_primer_pairs_from_fasta_file(self):
        stringency = 1
        chromosome = "1"
        pre_targeton_name = "ENSE00000769557_HG8_11"
        pre_targeton_start = 42929593
        pre_targeton_end = 42929803

        # arrange
        slices_fasta_file = self.fs.create_file(
            'fasta.fa',
            contents=f'>{pre_targeton_name}::{chromosome}:{pre_targeton_start}-{pre_targeton_end}(-)\nCACCTTCCCTCCGGTCCCCCCAGTGCTAAAGAAGCTGCGCGGGACAGCTGACGTGACCCATGACCTGCAGGAGATGAAGGAAGAGAGTCGGCAGATGATGCGGGAGAAGAAGGTCACCATCCTGGAGCTGTTCCGCTCCCCCGCCTACCGCCAGCCCATCCTCATCGCTGTGGTGCTGCAGCTGTCCCAGCAGCTGTCTGGCATCAACGC')

        designer_config = {"stringency_vector": [stringency]}

        p3_config = {
            "PRIMER_TASK": "pick_cloning_primers",
            "PRIMER_PICK_LEFT_PRIMER": 1,
            "PRIMER_PICK_RIGHT_PRIMER": 1,
            "PRIMER_OPT_SIZE": 20,
            "PRIMER_MIN_SIZE": 18,
            "PRIMER_MAX_SIZE": 23,
            "P3_FILE_FLAG": 1,
            "SEQUENCE_INCLUDED_REGION": [0, 200],
            "PRIMER_EXPLAIN_FLAG": 1
        }

        # act
        result = Primer3(designer_config, p3_config).get_primers(slices_fasta_file.name)

        # assert
        expected_primer_pair = PrimerPair(
            pair_id=f'{pre_targeton_name}_0_str1',
            chromosome=chromosome,
            pre_targeton_start=f'{pre_targeton_start}',
            pre_targeton_end=f'{pre_targeton_end}'
        )
        expected_forward = {'primer': f'{pre_targeton_name}_LibAmpF_0', 'penalty': 2.7456977357412597,
                            'side': 'right', 'stringency': f'{stringency}',
                            'pair_id': f'{pre_targeton_name}_0_str{stringency}', 'sequence': 'CAGACAGCTGCTGGGACA',
                            'coords': [199, 18],
                            'primer_start': 42929775, 'primer_end': 42929793, 'strand': '+', 'tm': 59.25430226425874,
                            'gc_percent': 61.111111111111114, 'self_any_th': 30.996860910464648, 'self_end_th': 0.0,
                            'hairpin_th': 35.513327628973116, 'end_stability': 4.02}
        expected_reverse = {'primer': f'{pre_targeton_name}_LibAmpR_0', 'penalty': 3.400054355094312, 'side': 'left',
                            'stringency': f'{stringency}', 'pair_id': f'{pre_targeton_name}_0_str{stringency}',
                            'sequence': 'CACCTTCCCTCCGGTCCC', 'coords': [0, 18],
                            'primer_start': 42929593, 'primer_end': 42929611, 'strand': '-', 'tm': 61.40005435509431,
                            'gc_percent': 72.22222222222223, 'self_any_th': 0.0, 'self_end_th': 0.0, 'hairpin_th': 0.0,
                            'end_stability': 4.46}
        expected_primer_pair.forward = expected_forward
        expected_primer_pair.reverse = expected_reverse

        self.assertEqual(result, [expected_primer_pair])


if __name__ == '__main__':
    unittest.main()
