from unittest import TestCase
from unittest.mock import patch, MagicMock, Mock

from primer.filter.designed_primer import DesignedPrimer
from primer.filter.duplicates_filter import DuplicatesFilter, _get_max_stringency_pair, \
    _take_pair_with_max_stringency_and_others, _group_duplicates_pairs
from primer.filter.filter_response import PrimerPairDiscarded
from primer.primer_pair import PrimerPair
import pandas as pd


class TestDuplicatesFilter(TestCase):

    def setUp(self) -> None:
        self.test_instance = DuplicatesFilter()

    @patch('primer.filter.duplicates_filter._group_duplicates_pairs')
    @patch('primer.filter.duplicates_filter._take_pair_with_max_stringency_and_others')
    def test_apply_duplicates_filter(self, pair_with_max_stringency_and_others, pair_groups):
        pair_max_stringency_group1 = MagicMock(spec=PrimerPair)
        pair_to_discard_group1 = MagicMock(spec=PrimerPair)

        pair_max_stringency_group2 = MagicMock(spec=PrimerPair)

        pair_groups.return_value = [[pair_max_stringency_group1, pair_to_discard_group1], [pair_max_stringency_group2]]

        pair_with_max_stringency_and_others.side_effect = [(pair_max_stringency_group1, [pair_to_discard_group1]),
                                                           (pair_max_stringency_group2, [])]

        result = self.test_instance.apply(
            [pair_max_stringency_group1, pair_to_discard_group1, pair_max_stringency_group2])

        self.assertEqual(result.primer_pairs_to_keep, [pair_max_stringency_group1, pair_max_stringency_group2])

        discarded_pair1 = PrimerPairDiscarded(pair_to_discard_group1, DuplicatesFilter.reason_discarded)
        self.assertEqual(result.primer_pairs_to_discard, [discarded_pair1])

    def test_apply_test_when_no_primer_pairs(self):
        response = self.test_instance.apply([])

        self.assertEqual(response.primer_pairs_to_keep, [])
        self.assertEqual(response.primer_pairs_to_discard, [])


class TestDuplicatesFilterAuxFunctions(TestCase):

    def test_get_max_stringency_pair(self):
        max_stringency_pair = MagicMock(spec=PrimerPair, forward=Mock(stringency=0.75))
        pair1 = MagicMock(spec=PrimerPair, forward=Mock(stringency=0.5))
        pair2 = MagicMock(spec=PrimerPair, forward=Mock(stringency=0.25))

        result = _get_max_stringency_pair([pair1, pair2, max_stringency_pair])

        self.assertEqual(result, max_stringency_pair)

    def test_take_pair_with_max_stringency_and_others(self):
        max_stringency_pair = MagicMock(spec=PrimerPair, forward=Mock(stringency=0.75))
        pair1 = MagicMock(spec=PrimerPair, forward=Mock(stringency=0.5))
        pair2 = MagicMock(spec=PrimerPair, forward=Mock(stringency=0.25))

        result_max_stringency_pair, other_pairs = _take_pair_with_max_stringency_and_others(
            [pair1, pair2, max_stringency_pair])

        self.assertEqual(result_max_stringency_pair, max_stringency_pair)
        self.assertEqual(other_pairs, [pair1, pair2])

    # TODO
    # def test_group_duplicates_pairs(self):
