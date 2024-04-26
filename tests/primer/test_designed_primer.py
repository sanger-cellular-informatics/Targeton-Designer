from unittest import TestCase

from primer.designed_primer import DesignedPrimer, Interval, map_to_designed_primer


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
            strand="+",
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
            strand="+",
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
            strand="+",
            tm=60.0,
            gc_percent=50.0,
            self_any_th=30.0,
            self_end_th=10.0,
            hairpin_th=20.0,
            end_stability=25.0
        )

        primer_with_different_penalty = DesignedPrimer(
            name="Primer1",
            penalty=1.5,
            pair_id="Pair1",
            sequence="ATCGATCG",
            coords=Interval(start=100, end=200),
            primer_start=100,
            primer_end=108,
            strand="+",
            tm=60.0,
            gc_percent=50.0,
            self_any_th=30.0,
            self_end_th=10.0,
            hairpin_th=20.0,
            end_stability=25.0
        )

        self.assertNotEquals(primer, primer_with_different_penalty)

    def test_map_to_designed_primer(self):
        primer_dict = {
            'primer': 'Primer1',
            'penalty': 0.5,
            'side': 'right',
            'pair_id': 'Pair1',
            'sequence': 'ATCGATCG',
            'coords': [199, 18],
            'primer_start': 119,
            'primer_end': 18,
            'strand': '+',
            'tm': 60.0,
            'gc_percent': 50.0,
            'self_any_th': 30.0,
            'self_end_th': 10.0,
            'hairpin_th': 20.0,
            'end_stability': 25.0
        }

        result = map_to_designed_primer(primer_dict)

        expected = DesignedPrimer(
            name="Primer1",
            penalty=0.5,
            pair_id="Pair1",
            sequence="ATCGATCG",
            coords=Interval(start=199, end=18),
            primer_start=119,
            primer_end=18,
            strand="+",
            tm=60.0,
            gc_percent=50.0,
            self_any_th=30.0,
            self_end_th=10.0,
            hairpin_th=20.0,
            end_stability=25.0
        )

        self.assertEqual(result, expected)
