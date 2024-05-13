import primer3

from typing import List

from primer.slice_data import SliceData
from primer.primer3_prepare_config import prepare_p3_config
from primer.primer_pair import PrimerPair, build_primer_pairs
from utils.exceptions import Primer3Error


class Primer3:
    def __init__(
            self,
            designer_config: dict,
            p3_config: dict
    ) -> None:

        self._p3_config = p3_config
        self._stringency_vector = designer_config.get('stringency_vector', [""])

    def get_primers(self, fasta: str) -> List[PrimerPair]:
        primer_pairs = []

        print('Reading Fasta file')
        slices = SliceData.parse_fasta(fasta)

        for slice in slices:
            slice_primer_pairs = self._get_primer_pairs(slice)
            primer_pairs.extend(slice_primer_pairs)

        return primer_pairs

    def _get_primer_pairs(self, slice_data: SliceData) -> List[PrimerPair]:
        primer_pairs = []
        primer_explain = []

        for stringency in self._stringency_vector:
            designs = self._get_primer3_designs(slice_data.p3_input, stringency)

            number_pairs = designs['PRIMER_PAIR_NUM_RETURNED']
            primer_explain_flag = self._p3_config['PRIMER_EXPLAIN_FLAG']
             
            # If Primer3 does not return any primer pairs for this stringency
            if not number_pairs:
                if primer_explain_flag:
                    msg = 'Stringency level ' + str(stringency) + " -- "
                    msg += self._get_primer3_explain(designs, stringency)
                else:
                    msg = 'Stringency level ' + str(stringency) + " -- "
                    msg += 'No primer pairs returned; add PRIMER_EXPLAIN_FLAG == 1 to config file for more details'
                
                primer_explain.append(msg)

            else:
                built_primer_pairs = build_primer_pairs(designs, slice_data, stringency)
                primer_pairs.extend(built_primer_pairs)

        # If Primer3 did not return any primer pairs for at least one stringency
        if primer_explain:
            message = '\n'.join([msg for msg in primer_explain])

            if len(primer_explain) == len(self._stringency_vector):
                message = 'NO PRIMER PAIRS BUILT BY PRIMER3: \n' + message
                raise Primer3Error(message)

            else:
                message = 'Warning: No primer pairs built by Primer3 with the following stringencies: \n' + message
                print(message)

        return primer_pairs

    def _get_primer3_designs(self, slice_info: dict, stringency) -> dict:
        config_data = prepare_p3_config(self._p3_config, stringency)
        return primer3.bindings.design_primers(slice_info, config_data)

    def _get_primer3_explain(self, designs, stringency) -> dict:
        keys = ["PRIMER_LEFT_EXPLAIN",
                "PRIMER_RIGHT_EXPLAIN",
                "PRIMER_PAIR_EXPLAIN"]
        msg = dict((key, designs[key]) for key in keys)
        msg_formatted = '; '.join([f"{key}: {value}" for key, value in msg.items()])
        return msg_formatted
