from unittest import TestCase
from unittest.mock import patch, MagicMock

from primer.filter.duplicates_filter import DuplicatesFilter
from primer.filter.filter_response import PrimerPairDiscarded
from primer.primer_pair import PrimerPair


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
