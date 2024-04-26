import primer3

from typing import List

from primer.slice_data import SliceData
from primer.primer3_prepare_config import prepare_p3_config
from primer.primer_pair import PrimerPair, build_primer_pairs


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

        for stringency in self._stringency_vector:
            designs = self._get_primer3_designs(slice_data.p3_input, stringency)
            built_primer_pairs = build_primer_pairs(designs, slice_data, str(stringency))

            primer_pairs.extend(built_primer_pairs)

        return primer_pairs


    def _get_primer3_designs(self, slice_info: dict, stringency) -> dict:
        config_data = prepare_p3_config(self._p3_config, stringency)
        return primer3.bindings.design_primers(slice_info, config_data)
