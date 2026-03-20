import unittest

from pyfakefs.fake_filesystem_unittest import TestCase
from parameterized import parameterized
from unittest.mock import patch, Mock

from primer.designed_primer import Interval, DesignedPrimer, Orientation, Strand
from primer.primer_pair import PrimerPair
from primer.slice_data import SliceData
from primer.build_primer_pairs import _calculate_primer_coords, build_primer_pairs


class TestBuildPrimerPair(TestCase):

    def test_build_primer_pairs_when_no_primer_pairs_found(self):
        designs = {
            'PRIMER_LEFT': [],
            'PRIMER_PAIR': [],
            'PRIMER_RIGHT': []
        }

        result = build_primer_pairs(design=designs, slice_data=Mock(), stringency=0.5, primer_type="LibAmp")

        self.assertEqual(result, [])

    def test_build_primer_pairs_when_primer_pairs_found(self):
        slice_data = SliceData(
            name="AABB",
            start=50398701,
            end=50399203,
            strand="-",
            chromosome="19",
            bases="CAGGCTGCACGGGGTGGAAGGGCAGTCCTTTGTAGGATTCAAAATTGGCCAGAGTGGGGTCACAGGTCAGCTTTGCAGCTGACTCTTG"
                  "GACCAAAGGCAAGACCAGGGGCAATGGGTTGGGCAGCAGGAACCTCCAACTCCAAGCCTTACCGTCTGCAACCCCCTCCAGGACTGAC"
                  "TGCAGCTCCTCCTCCTCCTGCTCCTGCAGCCTGTGTTCTGCCTCCATCTCCTCCATCAGTGCCAGGTCCTCCTCGAATTGGGATGGCC"
                  "GAGGTGCATCATCATCATCCCAGAGGCCCCCACGGGCCCGCTTTGGGGGCACCCCGGGCCCTGGGCCTGGCCGCCGCTTGCCATCCAT"
                  "CCTGCTGGGCAAGTTGGAGCTTGGTGGAGGTTCTGACCGGAGACACCTCTTACCCTGTCTCCATGGCATCCCATGCACCCACCATGAC"
                  "CACGAGCCTTTTCCTACTCTGGACTGGCTGTACTAACATAGCTCTCTGCACCCTGCAGGCAGT",
            flanking_region=150,
            exclusion_region=5
        )
        design = {
            "PRIMER_PAIR": [
                {"PENALTY": 0.3122182485828068, "COMPL_ANY_TH": 6.098089414466301, "COMPL_END_TH": 2.0852502864723874,
                 "PRODUCT_SIZE": 279, "PRODUCT_TM": 92.61358067328065},
                {"PENALTY": 0.3125104180142119, "COMPL_ANY_TH": 6.098089414466301, "COMPL_END_TH": 2.0852502864723874,
                 "PRODUCT_SIZE": 280, "PRODUCT_TM": 92.52731595233853}],
            "PRIMER_LEFT": [
                {"PENALTY": 0.16915001782036365, "SEQUENCE": "GCAGGAACCTCCAACTCCAA", "COORDS": [123, 20],
                 "TM": 59.88951398263703, "GC_PERCENT": 55.0, "SELF_ANY_TH": 0.0, "SELF_END_TH": 0.0,
                 "HAIRPIN_TH": 39.610015792775926, "END_STABILITY": 3.53},
                {"PENALTY": 0.16915001782036365, "SEQUENCE": "GCAGGAACCTCCAACTCCAA", "COORDS": [123, 20],
                 "TM": 59.88951398263703, "GC_PERCENT": 55.0, "SELF_ANY_TH": 0.0, "SELF_END_TH": 0.0,
                 "HAIRPIN_TH": 39.610015792775926, "END_STABILITY": 3.53}],
            "PRIMER_RIGHT": [
                {"PENALTY": 0.14306823076244313, "SEQUENCE": "AGAGGTGTCTCCGGTCAGAA", "COORDS": [401, 20],
                 "TM": 59.887393730886686, "GC_PERCENT": 55.0, "SELF_ANY_TH": 0.0, "SELF_END_TH": 0.0,
                 "HAIRPIN_TH": 40.83258802708207, "END_STABILITY": 3.02},
                {"PENALTY": 0.14336040019384827, "SEQUENCE": "AAGAGGTGTCTCCGGTCAGA", "COORDS": [402, 20],
                 "TM": 59.887393730886686, "GC_PERCENT": 55.0, "SELF_ANY_TH": 0.0, "SELF_END_TH": 0.0,
                 "HAIRPIN_TH": 40.83258802708207, "END_STABILITY": 3.27}]
        }

        result = build_primer_pairs(design=design, slice_data=slice_data, stringency=0.1, primer_type="LibAmp")

        primer1 = PrimerPair(pair_id='AABB_LibAmp_0_str01', uid='uid', chromosome='19',
                             pre_targeton_start=50398701, pre_targeton_end=50399203, product_size=279,
                             stringency=0.1, targeton_id='AABB'
                             )
        forward1 = DesignedPrimer(name='AABB_LibAmpF_0_str01', penalty=0.14306823076244313,
                                  pair_id='AABB_LibAmp_0_str01', sequence='AGAGGTGTCTCCGGTCAGAA',
                                  coords=Interval(start=401, end=20), primer_start=50398802,
                                  primer_end=50398821, strand=Strand.POSITIVE, tm=59.887393730886686,
                                  gc_percent=55.0, self_any_th=0.0, self_end_th=0.0,
                                  hairpin_th=40.83258802708207, end_stability=3.02, orientation=Orientation.FORWARD,)
        reverse1 = DesignedPrimer(name='AABB_LibAmpR_0_str01', penalty=0.16915001782036365,
                                  pair_id='AABB_LibAmp_0_str01', sequence='GCAGGAACCTCCAACTCCAA',
                                  coords=Interval(start=123, end=20), primer_start=50399061,
                                  primer_end=50399080, strand=Strand.NEGATIVE, tm=59.88951398263703,
                                  gc_percent=55.0, self_any_th=0.0, self_end_th=0.0,
                                  hairpin_th=39.610015792775926, end_stability=3.53, orientation=Orientation.REVERSE,)

        primer2 = PrimerPair(pair_id='AABB_LibAmp_1_str01', uid='uid', chromosome='19',
                             pre_targeton_start=50398701, pre_targeton_end=50399203, product_size=280, stringency=0.1,
                             targeton_id='AABB',
                             )
        forward2 = DesignedPrimer(name='AABB_LibAmpF_1_str01', penalty=0.14336040019384827,
                                  pair_id='AABB_LibAmp_1_str01',
                                  sequence='AAGAGGTGTCTCCGGTCAGA',
                                  coords=Interval(start=402, end=20),
                                  primer_start=50398801, primer_end=50398820,
                                  strand=Strand.POSITIVE, tm=59.887393730886686, gc_percent=55.0,
                                  self_any_th=0.0, self_end_th=0.0,
                                  hairpin_th=40.83258802708207, end_stability=3.27, orientation=Orientation.FORWARD,)
        reverse2 = DesignedPrimer(name='AABB_LibAmpR_1_str01', penalty=0.16915001782036365,
                                  pair_id='AABB_LibAmp_1_str01', sequence='GCAGGAACCTCCAACTCCAA',
                                  coords=Interval(start=123, end=20), primer_start=50399061,
                                  primer_end=50399080, strand=Strand.NEGATIVE, tm=59.88951398263703, gc_percent=55.0,
                                  self_any_th=0.0, self_end_th=0.0, hairpin_th=39.610015792775926,
                                  end_stability=3.53, orientation=Orientation.REVERSE,)
        primer1.forward = forward1
        primer1.reverse = reverse1
        primer2.forward = forward2
        primer2.reverse = reverse2

        self.assertEqual(len(result), 2)
        self.assertEqual(result[0].id, primer1.id)
        self.assertEqual(result[0].chromosome, primer1.chromosome)
        self.assertEqual(result[0].pre_targeton_start, primer1.pre_targeton_start)
        self.assertEqual(result[0].pre_targeton_end, primer1.pre_targeton_end)
        self.assertEqual(result[0].product_size, primer1.product_size)
        self.assertEqual(result[0].stringency, primer1.stringency)
        self.assertEqual(result[0].targeton_id, primer1.targeton_id)

        self.assertEqual(result[0].forward.name, forward1.name)
        self.assertEqual(result[0].forward.pair_id, forward1.pair_id)
        self.assertEqual(result[0].forward, forward1)
        self.assertEqual(result[0].reverse.name, reverse1.name)
        self.assertEqual(result[0].reverse.pair_id, reverse1.pair_id)
        self.assertEqual(result[0].reverse, reverse1)

        self.assertEqual(result[1].id, primer2.id)
        self.assertEqual(result[1].chromosome, primer2.chromosome)
        self.assertEqual(result[1].pre_targeton_start, primer2.pre_targeton_start)
        self.assertEqual(result[1].pre_targeton_end, primer2.pre_targeton_end)
        self.assertEqual(result[1].product_size, primer2.product_size)
        self.assertEqual(result[1].stringency, primer2.stringency)
        self.assertEqual(result[1].targeton_id, primer2.targeton_id)

        self.assertEqual(result[1].forward.name, forward2.name)
        self.assertEqual(result[1].forward.pair_id, forward2.pair_id)
        self.assertEqual(result[1].forward, forward2)
        self.assertEqual(result[1].reverse.name, reverse2.name)
        self.assertEqual(result[1].reverse.pair_id, reverse2.pair_id)
        self.assertEqual(result[1].reverse, reverse2)


class TestCalculatePrimerCoords(unittest.TestCase):

    @parameterized.expand(
        [
            ('left', [55, 20], 44490254, 44490755, '+', (44490309, 44490328)),
            ('right', [168, 19], 44490254, 44490755, '+', (44490404, 44490422)),
            ('left', [39, 20], 44490254, 44490755, '-', (44490697, 44490716)),
            ('right', [167, 19], 44490254, 44490755, '-', (44490588, 44490606)),
        ])
    def test_calculate_primer_coords(
            self,
            side,
            coords,
            slice_start,
            slice_end,
            strand,
            expected_result
    ):
        result = _calculate_primer_coords(side, coords, slice_start, slice_end, strand)

        self.assertEqual(result, expected_result)


if __name__ == '__main__':
    unittest.main()
