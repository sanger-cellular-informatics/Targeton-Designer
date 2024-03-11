from unittest import TestCase

from primer.hap1 import contain_variant


class TestPrimer3(TestCase):

    def test_check_variant_in_region_when_in_region(self):
        # variant position 42930743
        result = contain_variant(chromosome='chr1', start=42930700, end=42930800)

        self.assertEqual(result, True)

    def test_check_variant_in_region_when_not_in_region(self):
        result = contain_variant(chromosome='chr1', start=100, end=200)

        self.assertEqual(result, False)
