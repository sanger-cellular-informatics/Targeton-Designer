import primer3

from typing import List

from primer.slice_data import SliceData
from primer.primer3_prepare_config import prepare_config
from primer.primer_pair import parse_designs_to_primer_pairs, PrimerPair


class Primer3:
    def __init__(
            self,
            designer_params: dict,
            p3_params: dict
    ) -> None:

        self._p3_params = p3_params
        self._stringency_vector = designer_params['stringency_vector']

    def get_primers(self, fasta: str) -> List[PrimerPair]:
        primer_pairs = []

        print('Reading Fasta file')
        slices = SliceData.parse_fasta(fasta)

        for stringency in self._stringency_vector:
            print('Designing primers, stringency: ', stringency)
            designs = self._primer3_run(slices, str(stringency))

            print('Parsing primer pairs: ', stringency)
            pairs = parse_designs_to_primer_pairs(designs)
            primer_pairs.extend(pairs)

        return primer_pairs

    def _primer3_run(self, slices: List[SliceData], stringency: str) -> List[SliceData]:

        config_data = prepare_config(self._p3_params, stringency)
        result = self._primer3_design(slices, config_data, stringency)

        return result

    def _primer3_design(self, slices: List[SliceData], config: dict, stringency: str) -> List[SliceData]:
        designs = []
        for slice in slices:
            primer3_input = slice.p3_input
            design = primer3.bindings.design_primers(primer3_input, config)

            design['stringency'] = stringency

            slice.designs.append(design)
            designs.append(slice)

        return designs

