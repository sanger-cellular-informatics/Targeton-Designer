from unittest import TestCase
from unittest.mock import patch, MagicMock, Mock

from primer.designed_primer import DesignedPrimer, Interval
from primer.filter.duplicates_filter import DuplicatesFilter, _get_max_stringency_pair, \
    _take_pair_with_max_stringency_and_others, _group_duplicates_pairs
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


class TestDuplicatesFilterAuxFunctions(TestCase):

    def test_get_max_stringency_pair(self):
        max_stringency_pair = MagicMock(spec=PrimerPair, stringency=0.75, forward=Mock())
        pair1 = MagicMock(spec=PrimerPair, stringency=0.5, forward=Mock())
        pair2 = MagicMock(spec=PrimerPair, stringency=0.25, forward=Mock())

        result = _get_max_stringency_pair([pair1, pair2, max_stringency_pair])

        self.assertEqual(result, max_stringency_pair)

    def test_take_pair_with_max_stringency_and_others(self):
        max_stringency_pair = MagicMock(spec=PrimerPair, stringency=0.75, forward=Mock())
        pair1 = MagicMock(spec=PrimerPair, stringency=0.5, forward=Mock())
        pair2 = MagicMock(spec=PrimerPair, stringency=0.25, forward=Mock())

        result_max_stringency_pair, other_pairs = _take_pair_with_max_stringency_and_others(
            [pair1, pair2, max_stringency_pair])

        self.assertEqual(result_max_stringency_pair, max_stringency_pair)
        self.assertEqual(other_pairs, [pair1, pair2])

    def test_group_duplicates_pairs(self):
        # Note: The 'name' attribute is excluded from the comparison in the DesignedPrimer
        # The 'pair_id' and 'stringency' attributes are excluded from the comparison of PrimerPair

        # pair1_group1 and pair2_group1 belongs to the same group because they only differ in stringency and name
        pair1_group1 = PrimerPair(pair_id="1",
                                  chromosome="chr1",
                                  pre_targeton_start=42930996,
                                  pre_targeton_end=42931206,
                                  product_size=131,
                                  stringency=0.5,
                                  targeton_id="ENSE")
        pair1_group1.forward = _designed_primer(name='forward', penalty=1)
        pair1_group1.reverse = _designed_primer(name='reverse', penalty=1)
        pair2_group1 = PrimerPair(pair_id="2",
                                  chromosome="chr1",
                                  pre_targeton_start=42930996,
                                  pre_targeton_end=42931206,
                                  product_size=131,
                                  stringency=0.75,
                                  targeton_id="ENSE")
        pair2_group1.forward = _designed_primer(name='forward', penalty=1)
        pair2_group1.reverse = _designed_primer(name='reverse', penalty=1)

        # pair_group2 is not part of group1 because it differs from the other pairs,
        # having one primer with a different penalty and penalty is included in the comparison in PrimerPair
        pair_group2 = PrimerPair(pair_id="3",
                                 chromosome="chr1",
                                 pre_targeton_start=42930996,
                                 pre_targeton_end=42931206,
                                 product_size=131,
                                 stringency=0.5,
                                 targeton_id="ENSE")
        pair_group2.forward = _designed_primer(name='forward', penalty=1)
        pair_group2.reverse = _designed_primer(name='reverse', penalty=0.7)

        result = _group_duplicates_pairs([pair1_group1, pair2_group1, pair_group2])

        self.assertEqual(result, [[pair1_group1, pair2_group1], [pair_group2]])


def _designed_primer(name: str, penalty: float):
    return DesignedPrimer(
        name=name,
        penalty=penalty,
        pair_id="pair_id",
        sequence="ATCGATCGATCG",
        coords=Interval(start=100, end=200),
        primer_start=100,
        primer_end=110,
        strand="+",
        tm=60.5,
        gc_percent=50.0,
        self_any_th=0.1,
        self_end_th=0.2,
        hairpin_th=0.3,
        end_stability=0.4
    )
