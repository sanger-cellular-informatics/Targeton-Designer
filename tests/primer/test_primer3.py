import unittest
from pyfakefs.fake_filesystem_unittest import TestCase

from primer.designed_primer import DesignedPrimer, Interval
from primer.primer_pair import PrimerPair
from primer.primer3 import Primer3
from utils.exceptions import Primer3Error
from src import primer
from io import StringIO
import sys


class IntegrationTestPrimer3(TestCase):

    def setUp(self):
        self.setUpPyfakefs()

    def test_get_primer_pairs_from_fasta_file(self):
        stringency = 1
        chromosome = "1"
        pre_targeton_name = "ENSE00000769557_HG8_11"
        pre_targeton_start = 42929593
        pre_targeton_end = 42929803
        targeton_id = "ENSE"

        # arrange
        slices_fasta_file = self.fs.create_file(
            'fasta.fa',
            contents=f'>{pre_targeton_name}::{chromosome}:{pre_targeton_start}-{pre_targeton_end}(-)\nCACCTTCCCTCCGGTCCCCCCAGTGCTAAAGAAGCTGCGCGGGACAGCTGACGTGACCCATGACCTGCAGGAGATGAAGGAAGAGAGTCGGCAGATGATGCGGGAGAAGAAGGTCACCATCCTGGAGCTGTTCCGCTCCCCCGCCTACCGCCAGCCCATCCTCATCGCTGTGGTGCTGCAGCTGTCCCAGCAGCTGTCTGGCATCAACGC')

        designer_config = {"stringency_vector": [stringency]}

        p3_config = {
             "PRIMER_TASK": "pick_cloning_primers",
             "PRIMER_PICK_LEFT_PRIMER": 1,
             "PRIMER_PICK_RIGHT_PRIMER": 1,
             "PRIMER_OPT_SIZE": 20,
             "PRIMER_MIN_SIZE": 18,
             "PRIMER_MAX_SIZE": 23,
             "P3_FILE_FLAG": 1,
             "SEQUENCE_INCLUDED_REGION": [0, 200],
             "PRIMER_EXPLAIN_FLAG": 1
             }

        # act
        result = Primer3(designer_config, p3_config).get_primers(slices_fasta_file.name)

        uid = 'bc09fcac-07c0-11ef-b244-fa163e9abfe1'
        # assert
        expected_primer_pair = PrimerPair(
            uid={uid},
            pair_id=f'{pre_targeton_name}_0_str1',
            chromosome=chromosome,
            pre_targeton_start=f'{pre_targeton_start}',
            pre_targeton_end=f'{pre_targeton_end}',
            product_size=200,
            stringency=1,
            targeton_id="ENSE"
        )

        expected_forward = DesignedPrimer(
            name='ENSE00000769557_HG8_11_LibAmpF_0',
            penalty=2.7456977357412597,
            pair_id='ENSE00000769557_HG8_11_0_str1',
            sequence='CAGACAGCTGCTGGGACA',
            coords=Interval(start=199, end=18),
            primer_start=42929775,
            primer_end=42929793,
            strand='+',
            tm=59.25430226425874,
            gc_percent=61.111111111111114,
            self_any_th=30.996860910464648,
            self_end_th=0.0,
            hairpin_th=35.513327628973116,
            end_stability=4.02
        )

        expected_reverse = DesignedPrimer(
            name='ENSE00000769557_HG8_11_LibAmpR_0',
            penalty=3.400054355094312,
            pair_id='ENSE00000769557_HG8_11_0_str1',
            sequence='CACCTTCCCTCCGGTCCC',
            coords=Interval(start=0, end=18),
            primer_start=42929593,
            primer_end=42929611,
            strand='-',
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

    def test_get_primer_pairs_when_primer3_error(self):
        stringency = 1
        chromosome = "1"
        pre_targeton_name = "ARTY"
        pre_targeton_start = 42958479
        pre_targeton_end = 42958806

        # arrange
        slices_fasta_file = self.fs.create_file(
            'fasta.fa',
            contents=f'>{pre_targeton_name}::{chromosome}:{pre_targeton_start}-{pre_targeton_end}(+)\nGCTCGGGACCCGCACCGAGCCAGGCTCGGAGAGGCGCGCGGCCCGCCCCGGGCGCACAGCGCAGCGGGGCGGCGGGGGAGGCCCTGGCCGGCGTAAGGCGGGCAGGAGTCTGCGCCTTTGTTCCTGGCGGGAGGGCCCGCGGGCGCGCGACTCACCTTGCTGCTGGGCTCCATGGCAGCGCTGCGCTGGTGGCTCTGGCTGCGCCGGGTACGCGGGTGGCGACGGGCGTGCGAGCGGCGCTCTCCCGCTCAGGCTCGTGCTCCGGTCCGGGGACTCCCACTGCGACTCTGACTCCGACCCCCGTCGTTTGGTCTCCTGCTCCCTGGCG')

        designer_config = {"stringency_vector": [stringency]}

        p3_config = {
             "PRIMER_TASK": "generic",
             "PRIMER_PICK_LEFT_PRIMER": 1,
             "PRIMER_PICK_RIGHT_PRIMER": 1,
             "PRIMER_OPT_SIZE": 20,
             "PRIMER_MIN_SIZE": 18,
             "PRIMER_MAX_SIZE": 30,
             "P3_FILE_FLAG": 1,
             "SEQUENCE_INCLUDED_REGION": [0, 212],
             "PRIMER_EXPLAIN_FLAG": 1
             }

        # act
        with self.assertRaises(Primer3Error) as primer_error:
            Primer3(designer_config, p3_config).get_primers(slices_fasta_file.name)

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



if __name__ == '__main__':
    unittest.main()
