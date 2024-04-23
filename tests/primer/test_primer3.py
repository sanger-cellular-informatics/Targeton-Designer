import unittest
from pyfakefs.fake_filesystem_unittest import TestCase

from primer.primer_pair import PrimerPair
from primer.primer3 import Primer3


class IntegrationTestPrimer3(TestCase):

    def setUp(self):
        self.setUpPyfakefs()

    def test_get_primer_pairs_from_fasta_file(self):
        stringency = 1
        chromosome = "1"
        pre_targeton_name = "ENSE00000769557_HG8_11"
        pre_targeton_start = 42929593
        pre_targeton_end = 42929803

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

        # assert
        expected_primer_pair = PrimerPair(
            pair_id=f'{pre_targeton_name}_0_str1',
            chromosome=chromosome,
            pre_targeton_start=f'{pre_targeton_start}',
            pre_targeton_end=f'{pre_targeton_end}',
            product_size=200
        )
        expected_forward = {'primer': f'{pre_targeton_name}_LibAmpF_0', 'penalty': 2.7456977357412597,
                            'side': 'right', 'stringency': f'{stringency}',
                            'pair_id': f'{pre_targeton_name}_0_str{stringency}', 'sequence': 'CAGACAGCTGCTGGGACA',
                            'coords': [199, 18],
                            'primer_start': 42929775, 'primer_end': 42929793, 'strand': '+', 'tm': 59.25430226425874,
                            'gc_percent': 61.111111111111114, 'self_any_th': 30.996860910464648, 'self_end_th': 0.0,
                            'hairpin_th': 35.513327628973116, 'end_stability': 4.02}
        expected_reverse = {'primer': f'{pre_targeton_name}_LibAmpR_0', 'penalty': 3.400054355094312, 'side': 'left',
                            'stringency': f'{stringency}', 'pair_id': f'{pre_targeton_name}_0_str{stringency}',
                            'sequence': 'CACCTTCCCTCCGGTCCC', 'coords': [0, 18],
                            'primer_start': 42929593, 'primer_end': 42929611, 'strand': '-', 'tm': 61.40005435509431,
                            'gc_percent': 72.22222222222223, 'self_any_th': 0.0, 'self_end_th': 0.0, 'hairpin_th': 0.0,
                            'end_stability': 4.46}
        expected_primer_pair.forward = expected_forward
        expected_primer_pair.reverse = expected_reverse

        self.assertEqual(result, [expected_primer_pair])


if __name__ == '__main__':
    unittest.main()
