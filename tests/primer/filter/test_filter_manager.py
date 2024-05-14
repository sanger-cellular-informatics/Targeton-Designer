from unittest import TestCase

from primer.designed_primer import DesignedPrimer, Interval
from primer.filter.duplicates_filter import DuplicatesFilter
from primer.filter.hap1_variant_filter import HAP1VariantFilter
from primer.filter.filter_manager import FilterManager
from primer.primer_pair import PrimerPair
from primer.primer_pair_discarded import PrimerPairDiscarded


class TestFilterManager(TestCase):
    def setUp(self) -> None:
        # There is a HAP1 Variant at 11542 in chromosome 1
        self.primer_with_variant = DesignedPrimer(
            name="primer_with_variant",
            penalty=0.5,
            pair_id="pair_id",
            sequence="ATCGATCG",
            coords=Interval(start=199, end=18),
            primer_start=11540,
            primer_end=11545,
            strand="+",
            tm=60.0,
            gc_percent=50.0,
            self_any_th=30.0,
            self_end_th=10.0,
            hairpin_th=20.0,
            end_stability=25.0
        )

        self.primer_with_no_variant = DesignedPrimer(
            name="primer_with_no_variant",
            penalty=0.5,
            pair_id="pair_id",
            sequence="ATCGATCG",
            coords=Interval(start=199, end=18),
            primer_start=10,
            primer_end=20,
            strand="+",
            tm=60.0,
            gc_percent=50.0,
            self_any_th=30.0,
            self_end_th=10.0,
            hairpin_th=20.0,
            end_stability=25.0
        )

    def test_apply_filters(self):
        # Arrange
        pair_with_variant = PrimerPair(
            pair_id="pair_with_hap1_variant",
            chromosome="1",
            pre_targeton_start=11540,
            pre_targeton_end=11545,
            product_size="200",
            stringency=0.1,
            targeton_id="targeton_id",
            uid="uid")
        pair_with_variant.forward = self.primer_with_variant
        pair_with_variant.reverse = self.primer_with_no_variant

        pair_with_no_variant = PrimerPair(
            pair_id="pair_with_no_variant",
            chromosome="1",
            pre_targeton_start=11540,
            pre_targeton_end=11545,
            product_size="200",
            stringency=1,
            targeton_id="targeton_id",
            uid="uid")
        pair_with_no_variant.forward = self.primer_with_no_variant
        pair_with_no_variant.reverse = self.primer_with_no_variant

        pair_max_stringency = PrimerPair(
            pair_id="pair_max_stringency",
            chromosome="2",
            pre_targeton_start=11540,
            pre_targeton_end=11545,
            product_size="200",
            stringency=1,
            targeton_id="targeton_id",
            uid="uid")
        pair_max_stringency.forward = self.primer_with_no_variant
        pair_max_stringency.reverse = self.primer_with_no_variant

        pair_min_stringency = PrimerPair(
            pair_id="pair_min_stringency",
            chromosome="2",
            pre_targeton_start=11540,
            pre_targeton_end=11545,
            product_size="200",
            stringency=0.1,
            targeton_id="targeton_id",
            uid="uid")
        pair_min_stringency.forward = self.primer_with_no_variant
        pair_min_stringency.reverse = self.primer_with_no_variant

        # Act
        pairs_to_filter = [pair_with_variant, pair_with_no_variant, pair_max_stringency, pair_min_stringency]
        filter_response = FilterManager().apply_filters(pairs_to_filter)

        # Assertion
        self.assertEqual(len(filter_response.primer_pairs_to_keep), 2)
        self.assertIn(pair_with_no_variant, filter_response.primer_pairs_to_keep)
        self.assertTrue(pair_max_stringency, filter_response.primer_pairs_to_keep)

        self.assertEqual(len(filter_response.primer_pairs_to_discard), 2)
        self.assertIn(PrimerPairDiscarded(pair_with_variant, reason_discarded=HAP1VariantFilter.reason_discarded),
                      filter_response.primer_pairs_to_discard)
        self.assertIn(PrimerPairDiscarded(pair_min_stringency, reason_discarded=DuplicatesFilter.reason_discarded),
                      filter_response.primer_pairs_to_discard)

    def test_apply_filters_when_all_primer_pairs_are_kept(self):
        # Arrange
        pair_with_no_duplicates_nor_variant1 = PrimerPair(
            pair_id="pair_with_no_duplicates_nor_variant1",
            chromosome="1",
            pre_targeton_start=11540,
            pre_targeton_end=11545,
            product_size="200",
            stringency=1,
            targeton_id="targeton_id",
            uid="uid")
        pair_with_no_duplicates_nor_variant1.forward = self.primer_with_no_variant
        pair_with_no_duplicates_nor_variant1.reverse = self.primer_with_no_variant

        pair_with_no_duplicates_nor_variant2 = PrimerPair(
            pair_id="pair_with_no_duplicates_nor_variant2",
            chromosome="2",
            pre_targeton_start=11540,
            pre_targeton_end=11545,
            product_size="200",
            stringency=1,
            targeton_id="targeton_id",
            uid="uid")
        pair_with_no_duplicates_nor_variant2.forward = self.primer_with_no_variant
        pair_with_no_duplicates_nor_variant2.reverse = self.primer_with_no_variant

        # Act
        pairs_to_filter = [pair_with_no_duplicates_nor_variant1, pair_with_no_duplicates_nor_variant2]
        filter_response = FilterManager().apply_filters(pairs_to_filter)

        # Assertion
        self.assertEqual(len(filter_response.primer_pairs_to_keep), 2)
        self.assertIn(pair_with_no_duplicates_nor_variant1, filter_response.primer_pairs_to_keep)
        self.assertTrue(pair_with_no_duplicates_nor_variant2, filter_response.primer_pairs_to_keep)

        self.assertEqual(len(filter_response.primer_pairs_to_discard), 0)

    def test_apply_filters_when_all_primer_pairs_are_discarded(self):
        # Arrange
        pair_with_variant = PrimerPair(
            pair_id="pair_with_variant",
            chromosome="1",
            pre_targeton_start=11540,
            pre_targeton_end=11545,
            product_size="200",
            stringency=1,
            targeton_id="targeton_id",
            uid="uid")
        pair_with_variant.forward = self.primer_with_variant
        pair_with_variant.reverse = self.primer_with_no_variant

        pair_duplicate = PrimerPair(
            pair_id="pair_duplicate",
            chromosome="1",
            pre_targeton_start=11540,
            pre_targeton_end=11545,
            product_size="200",
            stringency=0.1,
            targeton_id="targeton_id",
            uid="uid")
        pair_duplicate.forward = self.primer_with_variant
        pair_duplicate.reverse = self.primer_with_no_variant

        # Act
        pairs_to_filter = [pair_with_variant, pair_duplicate]
        filter_response = FilterManager().apply_filters(pairs_to_filter)

        # Assertion
        self.assertEqual(len(filter_response.primer_pairs_to_keep), 0)

        self.assertEqual(len(filter_response.primer_pairs_to_discard), 2)
        self.assertIn(PrimerPairDiscarded(pair_with_variant, reason_discarded=HAP1VariantFilter.reason_discarded),
                      filter_response.primer_pairs_to_discard)
        self.assertIn(PrimerPairDiscarded(pair_duplicate, reason_discarded=DuplicatesFilter.reason_discarded),
                      filter_response.primer_pairs_to_discard)

    def test_apply_filters_when_no_primer_pairs(self):
        # Act
        filter_response = FilterManager().apply_filters([])

        # Assertion
        self.assertEqual(len(filter_response.primer_pairs_to_keep), 0)
        self.assertEqual(len(filter_response.primer_pairs_to_discard), 0)
