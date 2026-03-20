from unittest import TestCase, main
from parameterized import parameterized
import unittest
from primer.designed_primer import DesignedPrimer, Interval, Orientation, Strand


class TestDesignedPrimers(TestCase):

    def test_designed_primers_only_with_different_name_are_equal(self):
        # ONLY the 'name' attribute is EXCLUDED from the comparison of DesignedPrimer
        primer = DesignedPrimer(
            name="Primer1",
            penalty=0.5,
            pair_id="Pair1",
            sequence="ATCGATCG",
            coords=Interval(start=100, end=200),
            primer_start=100,
            primer_end=108,
            strand=Strand.POSITIVE,
            orientation=Orientation.FORWARD,
            tm=60.0,
            gc_percent=50.0,
            self_any_th=30.0,
            self_end_th=10.0,
            hairpin_th=20.0,
            end_stability=25.0
        )

        primer_with_different_name = DesignedPrimer(
            name="Primer2",
            penalty=0.5,
            pair_id="Pair1",
            sequence="ATCGATCG",
            coords=Interval(start=100, end=200),
            primer_start=100,
            primer_end=108,
            strand=Strand.POSITIVE,
            orientation=Orientation.FORWARD,
            tm=60.0,
            gc_percent=50.0,
            self_any_th=30.0,
            self_end_th=10.0,
            hairpin_th=20.0,
            end_stability=25.0
        )

        self.assertEqual(primer, primer_with_different_name)

    def test_compare_different_designed_primers(self):
        # ONLY the 'name' attribute is EXCLUDED from the comparison of DesignedPrimer
        primer = DesignedPrimer(
            name="Primer1",
            penalty=0.5,
            pair_id="Pair1",
            sequence="ATCGATCG",
            coords=Interval(start=100, end=200),
            primer_start=100,
            primer_end=108,
            strand=Strand.POSITIVE,
            orientation=Orientation.FORWARD,
            tm=60.0,
            gc_percent=50.0,
            self_any_th=30.0,
            self_end_th=10.0,
            hairpin_th=20.0,
            end_stability=25.0,
        )

        primer_with_different_penalty = DesignedPrimer(
            name="Primer1",
            penalty=1.5,
            pair_id="Pair1",
            sequence="ATCGATCG",
            coords=Interval(start=100, end=200),
            primer_start=100,
            primer_end=108,
            strand=Strand.POSITIVE,
            orientation=Orientation.FORWARD,
            tm=60.0,
            gc_percent=50.0,
            self_any_th=30.0,
            self_end_th=10.0,
            hairpin_th=20.0,
            end_stability=25.0
        )

        self.assertNotEqual(primer, primer_with_different_penalty)


class TestOrientation(TestCase):

    @parameterized.expand([
        ("+", "left", Orientation.FORWARD),
        ("+", "right", Orientation.REVERSE),
        ("-", "left", Orientation.REVERSE),
        ("-", "right", Orientation.FORWARD),
    ])
    def test_valid_combinations(self, strand, side, expected):
        orientation = Orientation.from_side_and_strand(side, strand)

        self.assertEqual(orientation, expected)

    @parameterized.expand([
        ("invalid_strand", "*", "left"),
        ("invalid_side", "+", "center"),
        ("empty_side", "-", ""),
        ("both_invalid", "x", "y"),
    ])
    def test_invalid_combinations(self, name, strand, side):
        with self.assertRaises(ValueError):
            Orientation.from_side_and_strand(side, strand)


class TestPrimerStrands(TestCase):
    @parameterized.expand([
        ('left', '+', Strand.POSITIVE),
        ('left', '-', Strand.NEGATIVE),
        ('right', '+', Strand.NEGATIVE),
        ('right', '-', Strand.POSITIVE),
    ])
    def test_strands_primers(self, side, strand, expected):
        strand = Strand.from_side_and_slice_strand(side, strand)

        # assert
        self.assertEqual(strand, expected)


if __name__ == '__main__':
    unittest.main()
