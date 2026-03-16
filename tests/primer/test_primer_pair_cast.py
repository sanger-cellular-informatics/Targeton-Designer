import unittest

from pyfakefs.fake_filesystem_unittest import TestCase
from parameterized import parameterized
from unittest.mock import patch, Mock

from primer.designed_primer import Interval, DesignedPrimer
from primer.primer_pair import build_primer_pairs, PrimerPair
from primer.slice_data import SliceData


class TestPrimerPairCast(TestCase):

    def test_simple(self):
        slice_data = SliceData(
            name="AABB",
            start=50398701,
            end=50399203,
            strand="-",
            chromosome="19",
            bases="CAGGCTGCACGGGGTGGAAGGGCAGTCCTTTGTAGGATTCAAAATTGGCCAGAGTGGGGTCACAGGTCAGCTTTGCAGCTGACTCTTGGACCAAAGGCAAGACCAGGGGCAATGGGTTGGGCAGCAGGAACCTCCAACTCCAAGCCTTACCGTCTGCAACCCCCTCCAGGACTGACTGCAGCTCCTCCTCCTCCTGCTCCTGCAGCCTGTGTTCTGCCTCCATCTCCTCCATCAGTGCCAGGTCCTCCTCGAATTGGGATGGCCGAGGTGCATCATCATCATCCCAGAGGCCCCCACGGGCCCGCTTTGGGGGCACCCCGGGCCCTGGGCCTGGCCGCCGCTTGCCATCCATCCTGCTGGGCAAGTTGGAGCTTGGTGGAGGTTCTGACCGGAGACACCTCTTACCCTGTCTCCATGGCATCCCATGCACCCACCATGACCACGAGCCTTTTCCTACTCTGGACTGGCTGTACTAACATAGCTCTCTGCACCCTGCAGGCAGT",
            flanking_region=150,
            exclusion_region=5
        )

        design = {
            "PRIMER_LEFT_EXPLAIN": "considered 3744, low tm 18, high tm 529, not in any ok left region 3112, ok 85",
            "PRIMER_RIGHT_EXPLAIN": "considered 3744, low tm 43, high tm 440, not in any ok right region 3112, ok 149",
            "PRIMER_PAIR_EXPLAIN": "considered 79, unacceptable product size 55, ok 24", "PRIMER_LEFT_NUM_RETURNED": 3,
            "PRIMER_RIGHT_NUM_RETURNED": 3, "PRIMER_INTERNAL_NUM_RETURNED": 0, "PRIMER_PAIR_NUM_RETURNED": 3,
            "PRIMER_PAIR": [
                {"PENALTY": 0.3122182485828068, "COMPL_ANY_TH": 6.098089414466301, "COMPL_END_TH": 2.0852502864723874,
                 "PRODUCT_SIZE": 279, "PRODUCT_TM": 92.61358067328065},
                {"PENALTY": 0.3125104180142119, "COMPL_ANY_TH": 6.098089414466301, "COMPL_END_TH": 2.0852502864723874,
                 "PRODUCT_SIZE": 280, "PRODUCT_TM": 92.52731595233853},
                {"PENALTY": 0.44316111700200334, "COMPL_ANY_TH": 1.9912483019189722, "COMPL_END_TH": 1.2218210672552345,
                 "PRODUCT_SIZE": 278, "PRODUCT_TM": 92.70046600372598}], "PRIMER_LEFT": [
                {"PENALTY": 0.16915001782036365, "SEQUENCE": "GCAGGAACCTCCAACTCCAA", "COORDS": [123, 20],
                 "TM": 59.88951398263703, "GC_PERCENT": 55.0, "SELF_ANY_TH": 0.0, "SELF_END_TH": 0.0,
                 "HAIRPIN_TH": 39.610015792775926, "END_STABILITY": 3.53},
                {"PENALTY": 0.16915001782036365, "SEQUENCE": "GCAGGAACCTCCAACTCCAA", "COORDS": [123, 20],
                 "TM": 59.88951398263703, "GC_PERCENT": 55.0, "SELF_ANY_TH": 0.0, "SELF_END_TH": 0.0,
                 "HAIRPIN_TH": 39.610015792775926, "END_STABILITY": 3.53},
                {"PENALTY": 0.16915001782036365, "SEQUENCE": "GCAGGAACCTCCAACTCCAA", "COORDS": [123, 20],
                 "TM": 59.88951398263703, "GC_PERCENT": 55.0, "SELF_ANY_TH": 0.0, "SELF_END_TH": 0.0,
                 "HAIRPIN_TH": 39.610015792775926, "END_STABILITY": 3.53}], "PRIMER_RIGHT": [
                {"PENALTY": 0.14306823076244313, "SEQUENCE": "AGAGGTGTCTCCGGTCAGAA", "COORDS": [401, 20],
                 "TM": 59.887393730886686, "GC_PERCENT": 55.0, "SELF_ANY_TH": 0.0, "SELF_END_TH": 0.0,
                 "HAIRPIN_TH": 40.83258802708207, "END_STABILITY": 3.02},
                {"PENALTY": 0.14336040019384827, "SEQUENCE": "AAGAGGTGTCTCCGGTCAGA", "COORDS": [402, 20],
                 "TM": 59.887393730886686, "GC_PERCENT": 55.0, "SELF_ANY_TH": 0.0, "SELF_END_TH": 0.0,
                 "HAIRPIN_TH": 40.83258802708207, "END_STABILITY": 3.27},
                {"PENALTY": 0.2740110991816397, "SEQUENCE": "GAGGTGTCTCCGGTCAGAAC", "COORDS": [400, 20],
                 "TM": 59.755822102021625, "GC_PERCENT": 60.0, "SELF_ANY_TH": 0.0, "SELF_END_TH": 0.0,
                 "HAIRPIN_TH": 41.99310719958544, "END_STABILITY": 3.01}], "PRIMER_INTERNAL": [],
            "PRIMER_PAIR_0_PENALTY": 0.3122182485828068, "PRIMER_LEFT_0_PENALTY": 0.16915001782036365,
            "PRIMER_RIGHT_0_PENALTY": 0.14306823076244313, "PRIMER_LEFT_0_SEQUENCE": "GCAGGAACCTCCAACTCCAA",
            "PRIMER_RIGHT_0_SEQUENCE": "AGAGGTGTCTCCGGTCAGAA", "PRIMER_LEFT_0": [123, 20], "PRIMER_RIGHT_0": [401, 20],
            "PRIMER_LEFT_0_TM": 59.88951398263703, "PRIMER_RIGHT_0_TM": 59.887393730886686,
            "PRIMER_LEFT_0_GC_PERCENT": 55.0, "PRIMER_RIGHT_0_GC_PERCENT": 55.0, "PRIMER_LEFT_0_SELF_ANY_TH": 0.0,
            "PRIMER_RIGHT_0_SELF_ANY_TH": 0.0, "PRIMER_LEFT_0_SELF_END_TH": 0.0, "PRIMER_RIGHT_0_SELF_END_TH": 0.0,
            "PRIMER_LEFT_0_HAIRPIN_TH": 39.610015792775926, "PRIMER_RIGHT_0_HAIRPIN_TH": 40.83258802708207,
            "PRIMER_LEFT_0_END_STABILITY": 3.53, "PRIMER_RIGHT_0_END_STABILITY": 3.02,
            "PRIMER_PAIR_0_COMPL_ANY_TH": 6.098089414466301, "PRIMER_PAIR_0_COMPL_END_TH": 2.0852502864723874,
            "PRIMER_PAIR_0_PRODUCT_SIZE": 279, "PRIMER_PAIR_0_PRODUCT_TM": 92.61358067328065,
            "PRIMER_PAIR_1_PENALTY": 0.3125104180142119, "PRIMER_LEFT_1_PENALTY": 0.16915001782036365,
            "PRIMER_RIGHT_1_PENALTY": 0.14336040019384827, "PRIMER_LEFT_1_SEQUENCE": "GCAGGAACCTCCAACTCCAA",
            "PRIMER_RIGHT_1_SEQUENCE": "AAGAGGTGTCTCCGGTCAGA", "PRIMER_LEFT_1": [123, 20], "PRIMER_RIGHT_1": [402, 20],
            "PRIMER_LEFT_1_TM": 59.88951398263703, "PRIMER_RIGHT_1_TM": 59.887393730886686,
            "PRIMER_LEFT_1_GC_PERCENT": 55.0, "PRIMER_RIGHT_1_GC_PERCENT": 55.0, "PRIMER_LEFT_1_SELF_ANY_TH": 0.0,
            "PRIMER_RIGHT_1_SELF_ANY_TH": 0.0, "PRIMER_LEFT_1_SELF_END_TH": 0.0, "PRIMER_RIGHT_1_SELF_END_TH": 0.0,
            "PRIMER_LEFT_1_HAIRPIN_TH": 39.610015792775926, "PRIMER_RIGHT_1_HAIRPIN_TH": 40.83258802708207,
            "PRIMER_LEFT_1_END_STABILITY": 3.53, "PRIMER_RIGHT_1_END_STABILITY": 3.27,
            "PRIMER_PAIR_1_COMPL_ANY_TH": 6.098089414466301, "PRIMER_PAIR_1_COMPL_END_TH": 2.0852502864723874,
            "PRIMER_PAIR_1_PRODUCT_SIZE": 280, "PRIMER_PAIR_1_PRODUCT_TM": 92.52731595233853,
            "PRIMER_PAIR_2_PENALTY": 0.44316111700200334, "PRIMER_LEFT_2_PENALTY": 0.16915001782036365,
            "PRIMER_RIGHT_2_PENALTY": 0.2740110991816397, "PRIMER_LEFT_2_SEQUENCE": "GCAGGAACCTCCAACTCCAA",
            "PRIMER_RIGHT_2_SEQUENCE": "GAGGTGTCTCCGGTCAGAAC", "PRIMER_LEFT_2": [123, 20], "PRIMER_RIGHT_2": [400, 20],
            "PRIMER_LEFT_2_TM": 59.88951398263703, "PRIMER_RIGHT_2_TM": 59.755822102021625,
            "PRIMER_LEFT_2_GC_PERCENT": 55.0, "PRIMER_RIGHT_2_GC_PERCENT": 60.0, "PRIMER_LEFT_2_SELF_ANY_TH": 0.0,
            "PRIMER_RIGHT_2_SELF_ANY_TH": 0.0, "PRIMER_LEFT_2_SELF_END_TH": 0.0, "PRIMER_RIGHT_2_SELF_END_TH": 0.0,
            "PRIMER_LEFT_2_HAIRPIN_TH": 39.610015792775926, "PRIMER_RIGHT_2_HAIRPIN_TH": 41.99310719958544,
            "PRIMER_LEFT_2_END_STABILITY": 3.53, "PRIMER_RIGHT_2_END_STABILITY": 3.01,
            "PRIMER_PAIR_2_COMPL_ANY_TH": 1.9912483019189722, "PRIMER_PAIR_2_COMPL_END_TH": 1.2218210672552345,
            "PRIMER_PAIR_2_PRODUCT_SIZE": 278, "PRIMER_PAIR_2_PRODUCT_TM": 92.70046600372598}

        result = build_primer_pairs(design=design,
                                    slice_data=slice_data,
                                    stringency=1)

        primer1 = PrimerPair(pair_id='AABB_LibAmp_0_str1', uid='9f6b0fea-1e71-11f1-8d25-9aac72ca3ecb', chromosome='19',
                             pre_targeton_start=50398701, pre_targeton_end=50399203, product_size=279,
                             stringency=1, targeton_id='AABB'
                             )
        forward1 = DesignedPrimer(name='AABB_LibAmpF_0', penalty=0.14306823076244313,
                                  pair_id='AABB_LibAmp_0_str1', sequence='AGAGGTGTCTCCGGTCAGAA',
                                  coords=Interval(start=401, end=20), primer_start=50398802,
                                  primer_end=50398821, strand='+', tm=59.887393730886686,
                                  gc_percent=55.0, self_any_th=0.0, self_end_th=0.0,
                                  hairpin_th=40.83258802708207, end_stability=3.02)
        reverse1 = DesignedPrimer(name='AABB_LibAmpR_0', penalty=0.16915001782036365,
                                  pair_id='AABB_LibAmp_0_str1', sequence='GCAGGAACCTCCAACTCCAA',
                                  coords=Interval(start=123, end=20), primer_start=50399061,
                                  primer_end=50399080, strand='-', tm=59.88951398263703,
                                  gc_percent=55.0, self_any_th=0.0, self_end_th=0.0,
                                  hairpin_th=39.610015792775926, end_stability=3.53)

        primer2 = PrimerPair(pair_id='AABB_LibAmp_1_str1', uid='9f6b1986-1e71-11f1-8d25-9aac72ca3ecb', chromosome='19',
                             pre_targeton_start=50398701, pre_targeton_end=50399203, product_size=280, stringency=1,
                             targeton_id='AABB',
                             )
        forward2 = DesignedPrimer(name='AABB_LibAmpF_1', penalty=0.14336040019384827,
                                  pair_id='AABB_LibAmp_1_str1',
                                  sequence='AAGAGGTGTCTCCGGTCAGA',
                                  coords=Interval(start=402, end=20),
                                  primer_start=50398801, primer_end=50398820,
                                  strand='+', tm=59.887393730886686, gc_percent=55.0,
                                  self_any_th=0.0, self_end_th=0.0,
                                  hairpin_th=40.83258802708207, end_stability=3.27)
        reverse2 = DesignedPrimer(name='AABB_LibAmpR_1', penalty=0.16915001782036365,
                                  pair_id='AABB_LibAmp_1_str1', sequence='GCAGGAACCTCCAACTCCAA',
                                  coords=Interval(start=123, end=20), primer_start=50399061,
                                  primer_end=50399080, strand='-', tm=59.88951398263703, gc_percent=55.0,
                                  self_any_th=0.0, self_end_th=0.0, hairpin_th=39.610015792775926,
                                  end_stability=3.53)

        primer3 = PrimerPair(pair_id='AABB_LibAmp_2_str1', uid='9f6b2156-1e71-11f1-8d25-9aac72ca3ecb', chromosome='19',
                             pre_targeton_start=50398701, pre_targeton_end=50399203, product_size=278, stringency=1,
                             targeton_id='AABB'
                             )
        forward3 = DesignedPrimer(name='AABB_LibAmpF_2', penalty=0.2740110991816397,
                                  pair_id='AABB_LibAmp_2_str1',
                                  sequence='GAGGTGTCTCCGGTCAGAAC',
                                  coords=Interval(start=400, end=20),
                                  primer_start=50398803, primer_end=50398822,
                                  strand='+', tm=59.755822102021625, gc_percent=60.0,
                                  self_any_th=0.0, self_end_th=0.0,
                                  hairpin_th=41.99310719958544, end_stability=3.01)
        reverse3 = DesignedPrimer(name='AABB_LibAmpR_2', penalty=0.16915001782036365,
                                  pair_id='AABB_LibAmp_2_str1', sequence='GCAGGAACCTCCAACTCCAA',
                                  coords=Interval(start=123, end=20), primer_start=50399061,
                                  primer_end=50399080, strand='-', tm=59.88951398263703, gc_percent=55.0,
                                  self_any_th=0.0, self_end_th=0.0, hairpin_th=39.610015792775926,
                                  end_stability=3.53)

        self.assertEqual(len(result), 3)
        self.assertEqual(result[0].id, primer1.id)
        self.assertEqual(result[0].chromosome, primer1.chromosome)
        self.assertEqual(result[0].pre_targeton_start, primer1.pre_targeton_start)
        self.assertEqual(result[0].pre_targeton_end, primer1.pre_targeton_end)
        self.assertEqual(result[0].product_size, primer1.product_size)
        self.assertEqual(result[0].stringency, primer1.stringency)
        self.assertEqual(result[0].targeton_id, primer1.targeton_id)
        self.assertEqual(result[0].forward, forward1)
        self.assertEqual(result[0].reverse, reverse1)

        self.assertEqual(result[1].id, primer2.id)
        self.assertEqual(result[1].chromosome, primer2.chromosome)
        self.assertEqual(result[1].pre_targeton_start, primer2.pre_targeton_start)
        self.assertEqual(result[1].pre_targeton_end, primer2.pre_targeton_end)
        self.assertEqual(result[1].product_size, primer2.product_size)
        self.assertEqual(result[1].stringency, primer2.stringency)
        self.assertEqual(result[1].targeton_id, primer2.targeton_id)
        self.assertEqual(result[1].forward, forward2)
        self.assertEqual(result[1].reverse, reverse2)

        self.assertEqual(result[2].id, primer3.id)
        self.assertEqual(result[2].chromosome, primer3.chromosome)
        self.assertEqual(result[2].pre_targeton_start, primer3.pre_targeton_start)
        self.assertEqual(result[2].pre_targeton_end, primer3.pre_targeton_end)
        self.assertEqual(result[2].product_size, primer3.product_size)
        self.assertEqual(result[2].stringency, primer3.stringency)
        self.assertEqual(result[2].targeton_id, primer3.targeton_id)
        self.assertEqual(result[2].forward, forward3)
        self.assertEqual(result[2].reverse, reverse3)

if __name__ == '__main__':
    unittest.main()
