import primer3

from Bio.Seq import Seq
from typing import List

from utils.exceptions import InvalidConfigError
from utils.file_system import parse_json

from primer.slice_data import SliceData
from primer.primer3_prepare_config import prepare_config
from primer.primer_class import parse_designs_to_primers


class Primer3:
    def __init__(self, user_config: str = None) -> None:
        default_config = './src/primer/primer3.config.json'
        self._config = user_config or default_config

    def get_primers(self, fasta: str) -> List[dict]:
        print('Reading Fasta file')
        slices = SliceData.parse_fasta(fasta)

        print('Designing primers')
        designs = self._primer3_run(slices)

        print('Naming primers')
        slices_dict = parse_designs_to_primers(designs)

        return slices_dict

    def _primer3_run(self, slices: List[SliceData]) -> List[SliceData]:
       # STRINGENCY_VECTOR = ["0.1", "0.25", "0.5", "0.75", "1.0"]
        STRINGENCY_VECTOR = ["0.1", "1.0"]

        base_config = self._get_config_data()

        vector = STRINGENCY_VECTOR
        if base_config.get('PRIMER_MASK_FAILURE_RATE') is None:
            vector = [""]

        for stringency in vector:
            config_data = prepare_config(base_config, stringency)
            result = self._primer3_design(slices, config_data)

        return result

    def _primer3_design(self, slices: List[SliceData], config: dict) -> List[SliceData]:
        designs = []
        for slice in slices:
            primer3_input = slice.p3_input
            design = primer3.bindings.design_primers(primer3_input, config)

            stringency = ""
            if "PRIMER_MASK_FAILURE_RATE" in config:
                stringency = config["PRIMER_MASK_FAILURE_RATE"]

            design["stringency"] = stringency

            slice.designs.append(design)
            designs.append(slice)

        return designs

    def _get_config_data(self) -> dict:
        try:
            config_data = parse_json(self._config)
        except FileNotFoundError:
            raise
        except Exception:
            raise InvalidConfigError('Primer3 config file is not a correct JSON')

        return config_data
