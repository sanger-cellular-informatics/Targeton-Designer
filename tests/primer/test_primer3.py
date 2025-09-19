import unittest
from pyfakefs.fake_filesystem_unittest import TestCase
from unittest.mock import patch

from primer.designed_primer import DesignedPrimer, Interval
from primer.primer_pair import PrimerPair
from primer.primer3 import Primer3
from primer.slice_data import SliceData
from utils.exceptions import Primer3Error


class IntegrationTestPrimer3(TestCase):

    def setUp(self):
        self.setUpPyfakefs()

    def test_get_primer_pairs_from_fasta_file(self):
        stringency = 1
        chromosome = "1"
        pre_targeton_name = "ENSE00000769557_HG8_11"
        pre_targeton_start = 42929593
        pre_targeton_end = 42929803
        bases = "CACCTTCCCTCCGGTCCCCCCAGTGCTAAAGAAGCTGCGCGGGACAGCTGACGTGACCCATGACCTGCAGGAGATGAAGGAAGAGAGTCGGCAGATGATGCGGGAGAAGAAGGTCACCATCCTGGAGCTGTTCCGCTCCCCCGCCTACCGCCAGCCCATCCTCATCGCTGTGGTGCTGCAGCTGTCCCAGCAGCTGTCTGGCATCAACGC"

        # arrange
        slice = SliceData(
            name=pre_targeton_name,
            chromosome=chromosome,
            start=pre_targeton_start,
            end=pre_targeton_end,
            bases=bases,
            strand="-",
            region_padding=210,
            region_avoid=0
        )

        p3_config = {
             "PRIMER_TASK": "pick_cloning_primers",
             "PRIMER_PICK_LEFT_PRIMER": 1,
             "PRIMER_PICK_RIGHT_PRIMER": 1,
             "PRIMER_OPT_SIZE": 20,
             "PRIMER_MIN_SIZE": 18,
             "PRIMER_MAX_SIZE": 30,
             "P3_FILE_FLAG": 1,
             "PRIMER_PRODUCT_SIZE_RANGE": "100-210",
             "PRIMER_EXPLAIN_FLAG": 1,
             "PRIMER_MASK_TEMPLATE": 0,
             "PRIMER_NUM_RETURN": 1
             }

        # act
        result = Primer3(stringency_vector=[stringency], p3_config=p3_config).get_primers(slice)

        # assert
        expected_primer_pair = PrimerPair(
            uid="24ef4f66-9491-11f0-8314-8773ed291a1b",
            pair_id=f"{pre_targeton_name}_0_str1",
            chromosome=chromosome,
            pre_targeton_start=pre_targeton_start,
            pre_targeton_end=pre_targeton_end,
            product_size=210,
            stringency=stringency,
            targeton_id="ENSE"
        )

        expected_forward = DesignedPrimer(
            name=f"{pre_targeton_name}_LibAmpF_0",
            penalty=2.0371361439404154,
            pair_id=f"{pre_targeton_name}_0_str1",
            sequence="GCGTTGATGCCAGACAGCT",
            coords=Interval(start=209, end=19),
            primer_start=42929594,
            primer_end=42929612,
            strand="+",
            tm=61.037136143940415,
            gc_percent=57.89473684210526,
            self_any_th=0.0,
            self_end_th=0.0,
            hairpin_th=35.057108503623226,
            end_stability=4.24
        )

        expected_reverse = DesignedPrimer(
            name=f"{pre_targeton_name}_LibAmpR_0",
            penalty=3.400054355094312,
            pair_id=f"{pre_targeton_name}_0_str1",
            sequence="CACCTTCCCTCCGGTCCC",
            coords=Interval(start=0, end=18),
            primer_start=42929786,
            primer_end=42929803,
            strand="-",
            tm=61.40005435509431,
            gc_percent=72.22222222222223,
            self_any_th=0.0,
            self_end_th=0.0,
            hairpin_th=0.0,
            end_stability=4.46
        )

        expected_primer_pair.forward = expected_forward
        expected_primer_pair.reverse = expected_reverse
        self.assertEqual(result, [expected_primer_pair])

    @patch('primer3.bindings.design_primers')
    def test_get_primer_pairs_when_primer3_error(self, mock_design_primers):
        mock_design_primers.return_value = {
            "PRIMER_LEFT_EXPLAIN": "considered 1469, GC content failed 769, low tm 1, high tm 657, high hairpin stability 2, ok 40",
            "PRIMER_RIGHT_EXPLAIN": "considered 1469, GC content failed 235, low tm 1, high tm 1159, ok 74",
            "PRIMER_PAIR_EXPLAIN": "considered 2960, unacceptable product size 2960, ok 0",
            "PRIMER_PAIR_NUM_RETURNED": 0, "PRIMER_PAIR": []}

        # arrange
        bases = "GCTCGGGACCCGCACCGAGCCAGGCTCGGAGAGGCGCGCGGCCCGCCCCGGGCGCACAGCGCAGCGGGGCGGCGGGGGAGGCCCTGGCCGGCGTAAGGCGGGCAGGAGTCTGCGCCTTTGTTCCTGGCGGGAGGGCCCGCGGGCGCGCGACTCACCTTGCTGCTGGGCTCCATGGCAGCGCTGCGCTGGTGGCTCTGGCTGCGCCGGGTACGCGGGTGGCGACGGGCGTGCGAGCGGCGCTCTCCCGCTCAGGCTCGTGCTCCGGTCCGGGGACTCCCACTGCGACTCTGACTCCGACCCCCGTCGTTTGGTCTCCTGCTCCCTGGCG"
        slice = SliceData(
            name="ARTY",
            chromosome="1",
            start=42958479,
            end=42958806,
            bases=bases,
            strand="+",
            region_padding=328,
            region_avoid=0
        )

        p3_config = {
             "PRIMER_TASK": "generic",
             "PRIMER_PICK_LEFT_PRIMER": 1,
             "PRIMER_PICK_RIGHT_PRIMER": 1,
             "PRIMER_OPT_SIZE": 20,
             "PRIMER_MIN_SIZE": 18,
             "PRIMER_MAX_SIZE": 30,
             "P3_FILE_FLAG": 1,
             "PRIMER_PRODUCT_SIZE_RANGE": "250-328",
             "PRIMER_EXPLAIN_FLAG": 1,
             "PRIMER_MASK_TEMPLATE": 0
             }

        # act
        with self.assertRaises(Primer3Error) as primer_error:
            Primer3(stringency_vector=[1], p3_config=p3_config).get_primers(slice)

        # assert
        expected_msg = """\
            NO PRIMER PAIRS BUILT BY PRIMER3: \
            Stringency level 1 -- PRIMER_LEFT_EXPLAIN: considered 1469, GC content failed 769, low tm 1, \
            high tm 657, high hairpin stability 2, ok 40; PRIMER_RIGHT_EXPLAIN: considered 1469, GC \
            content failed 235, low tm 1, high tm 1159, ok 74; PRIMER_PAIR_EXPLAIN: considered 2960, \
            unacceptable product size 2960, ok 0"""
        expected_msg = ''.join(expected_msg.strip().split())

        error_msg = str(primer_error.exception)
        error_msg = ''.join(error_msg.strip().split())

        self.assertEqual(error_msg, expected_msg)

    @patch('primer.primer3.build_primer_pairs')
    def test_get_primer_pairs_empty_list(self, mock_build_primer_pairs):
        # arrange
        bases = "CACCTTCCCTCCGGTCCCCCCAGTGCTAAAGAAGCTGCGCGGGACAGCTGACGTGACCCATGACCTGCAGGAGATGAAGGAAGAGAGTCGGCAGATGATGCGGGAGAAGAAGGTCACCATCCTGGAGCTGTTCCGCTCCCCCGCCTACCGCCAGCCCATCCTCATCGCTGTGGTGCTGCAGCTGTCCCAGCAGCTGTCTGGCATCAACGC"
        slice = SliceData(
            name="ENSE00000769557_HG8_11",
            chromosome="1",
            start=42929593,
            end=42929803,
            bases=bases,
            strand="-",
            region_padding=210,
            region_avoid=0
        )

        mock_build_primer_pairs.return_value = []

        p3_config = {
             "PRIMER_TASK": "pick_cloning_primers",
             "PRIMER_PICK_LEFT_PRIMER": 1,
             "PRIMER_PICK_RIGHT_PRIMER": 1,
             "PRIMER_OPT_SIZE": 20,
             "PRIMER_MIN_SIZE": 18,
             "PRIMER_MAX_SIZE": 30,
             "P3_FILE_FLAG": 1,
             "PRIMER_PRODUCT_SIZE_RANGE": "100-210",
             "PRIMER_EXPLAIN_FLAG": 1,
             "PRIMER_MASK_TEMPLATE": 0
             }

        # act
        with self.assertRaises(ValueError) as value_error:
            Primer3(stringency_vector=[1], p3_config=p3_config).get_primers(slice)

        # assert
        self.assertEqual(str(value_error.exception), "No primer pairs returned")

    def test_kmer_lists_exist_success(self):
        # arrange
        p3_config = {
            "PRIMER_MASK_TEMPLATE": 1,
            "PRIMER_MASK_KMERLIST_PATH": "kmers/",
            }

        self.fs.create_dir("kmers/")
        self.fs.create_file("kmers/homo_sapiens_11.list", contents="list of 11-mers")
        self.fs.create_file("kmers/homo_sapiens_16.list", contents="list of 16-mers")

        # act
        result = Primer3(stringency_vector=[1], p3_config=p3_config)._kmer_lists_exist()

        # assert
        self.assertIsNone(result)

    def test_kmer_lists_exist_no_path(self):
        # arrange
        p3_config = {
            "PRIMER_MASK_TEMPLATE": 1,
            "PRIMER_MASK_KMERLIST_PATH": "kmers/",
            }

        # act
        with self.assertRaises(ValueError) as kmer_dir_error:
            Primer3(stringency_vector=[1], p3_config=p3_config)._kmer_lists_exist()

        # assert
        expected_result = "Missing directory with kmer lists required for masking. Expected path: 'kmers/'"
        result = str(kmer_dir_error.exception)
        self.assertEqual(result, expected_result)

    def test_kmer_lists_exist_no_lists(self):
        # arrange
        p3_config = {
            "PRIMER_MASK_TEMPLATE": 1,
            "PRIMER_MASK_KMERLIST_PATH": "kmers/",
            }

        self.fs.create_dir("kmers/")

        # act
        with self.assertRaises(ValueError) as kmer_lists_error:
            Primer3(stringency_vector=[1], p3_config=p3_config)._kmer_lists_exist()

        # assert
        expected_result = "Missing kmer list file(s) required for masking: 'kmers/homo_sapiens_11.list', 'kmers/homo_sapiens_16.list'"
        result = str(kmer_lists_error.exception)
        self.assertEqual(result, expected_result)


if __name__ == '__main__':
    unittest.main()
