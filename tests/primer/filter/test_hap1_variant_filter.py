import unittest
from unittest.mock import MagicMock
from primer.primer_pair import PrimerPair
from primer.filter.hap1_variant_filter import HAP1VariantFilter


class TestHAP1VariantFilter(unittest.TestCase):

    def setUp(self) -> None:
        self.test_instance = HAP1VariantFilter()

    def test_apply_hap_one_filter(self):
        pair1_with_variant = MagicMock(spec=PrimerPair)
        pair1_with_variant.contain_hap_one_variant = True
        pair2_with_variant = MagicMock(spec=PrimerPair)
        pair2_with_variant.contain_hap_one_variant = True

        pair_without_variant = MagicMock(spec=PrimerPair)
        pair_without_variant.contain_hap_one_variant = False

        response = self.test_instance.apply([pair1_with_variant, pair_without_variant, pair2_with_variant])

        self.assertEqual(len(response.primer_pairs_to_discard), 2)
        self.assertEqual(len(response.primer_pairs_to_keep), 1)
        self.assertTrue(pair_without_variant in response.primer_pairs_to_keep)
        self.assertTrue(
            pair1_with_variant in [pair_discarded.primer_pair for pair_discarded in response.primer_pairs_to_discard])
        self.assertEqual(response.primer_pairs_to_discard[0].filter_applied, HAP1VariantFilter.reason_discarded)
        self.assertTrue(
            pair2_with_variant in [pair_discarded.primer_pair for pair_discarded in response.primer_pairs_to_discard])
        self.assertEqual(response.primer_pairs_to_discard[1].filter_applied, HAP1VariantFilter.reason_discarded)

    def test_apply_test_when_no_primer_pairs(self):
        pair_without_variant = MagicMock(spec=PrimerPair)
        pair_without_variant.contain_hap_one_variant = False

        response = self.test_instance.apply([pair_without_variant])

        self.assertEqual(len(response.primer_pairs_to_keep), 1)
        self.assertEqual(len(response.primer_pairs_to_discard), 0)
        self.assertTrue(pair_without_variant in response.primer_pairs_to_keep)


if __name__ == '__main__':
    unittest.main()
