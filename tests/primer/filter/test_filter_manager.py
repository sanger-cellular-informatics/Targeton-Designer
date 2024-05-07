from unittest import TestCase

from primer.designed_primer import DesignedPrimer, Interval
from primer.filter import HAP1VariantFilter, DuplicatesFilter
from primer.filter.filter_manager import FilterManager
from primer.filter.filter_response import PrimerPairDiscarded
from primer.primer_pair import PrimerPair


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

    def test_filter_manager_init_when_valid_filter_names(self):
        test_instance = FilterManager(filter_names=["HAP1_variant", "duplicates"])

        self.assertEqual(len(test_instance._filters_to_apply), 2)

        has_hap1_variant_filter = any(isinstance(f, HAP1VariantFilter) for f in test_instance._filters_to_apply)
        has_duplicates_filter = any(isinstance(f, DuplicatesFilter) for f in test_instance._filters_to_apply)
        self.assertTrue(has_hap1_variant_filter)
        self.assertTrue(has_duplicates_filter)

    def test_filter_manager_init_when_wrong_filter_name(self):
        with self.assertRaises(ValueError) as exception_context:
            FilterManager(filter_names=["NON_EXISTING_FILTER"])

        self.assertEqual(str(exception_context.exception), "Filter name 'NON_EXISTING_FILTER' does not exist")

    def test_filter_manager_init_when_repeated_filter_name(self):
        test_instance = FilterManager(filter_names=["HAP1_variant", "HAP1_variant"])

        self.assertEqual(len(test_instance._filters_to_apply), 1)
        self.assertIsInstance(test_instance._filters_to_apply[0], HAP1VariantFilter)

    def test_filter_manager_init_when_no_filter_names(self):
        test_instance = FilterManager(filter_names=[])

        self.assertEqual(test_instance._filters_to_apply, [])

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
        filter_response = FilterManager(filter_names=["HAP1_variant", "duplicates"]).apply_filters(pairs_to_filter)

        # Assertion
        self.assertEqual(len(filter_response.primer_pairs_to_keep), 2)
        self.assertIn(pair_with_no_variant, filter_response.primer_pairs_to_keep)
        self.assertTrue(pair_max_stringency, filter_response.primer_pairs_to_keep)

        self.assertEqual(len(filter_response.primer_pairs_to_discard), 2)
        self.assertIn(PrimerPairDiscarded(pair_with_variant, reason_discarded=HAP1VariantFilter.reason_discarded),
                      filter_response.primer_pairs_to_discard)
        self.assertIn(PrimerPairDiscarded(pair_min_stringency, reason_discarded=DuplicatesFilter.reason_discarded),
                      filter_response.primer_pairs_to_discard)
