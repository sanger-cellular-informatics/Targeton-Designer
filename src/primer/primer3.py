import primer3

from Bio.Seq import Seq
from typing import List

from utils.exceptions import InvalidConfigError
from utils.file_system import parse_json

from primer.slice_data import SliceData
from primer.primer3_prepare_config import prepare_config
from primer.primer_pair import parse_designs_to_primer_pairs
from config.config import Config


DEFAULT_P3_CONFIG = './src/primer/primer3.config.json'

class Primer3:
    def __init__(
            self,
            designer_config: Config,
            user_p3_config: str = DEFAULT_P3_CONFIG
    ) -> None:
        print('DESIGNER_CONFIG::::', designer_config)

        self._p3_config = user_p3_config
        self._stringency_vector = designer_config.stringency_vector

    def get_primers(self, fasta: str) -> List[dict]:
        # STRINGENCY_VECTOR = ["0.1", "0.25", "0.5", "0.75", "1.0"]
        STRINGENCY_VECTOR = ["1.0", "0.1"]
        primer_pairs = []

        print('Reading Fasta file')
        slices = SliceData.parse_fasta(fasta)

        for stringency in STRINGENCY_VECTOR:
            print('Designing primers, stringency: ', stringency)
            designs = self._primer3_run(slices, stringency)

            print('Parsing primer pairs: ', stringency)
            pairs = parse_designs_to_primer_pairs(designs)
            primer_pairs.extend(pairs)

        return primer_pairs

    def _primer3_run(self, slices: List[SliceData], stringency: str) -> List[SliceData]:
        base_p3_config = self._get_config_data()

        config_data = prepare_config(base_p3_config, stringency)
        result = self._primer3_design(slices, config_data, stringency)

        return result

    def _primer3_design(self, slices: List[SliceData], config: dict, stringency: str) -> List[SliceData]:
        designs = []
        for slice in slices:
            primer3_input = slice.p3_input
            design = primer3.bindings.design_primers(primer3_input, config)

            design["stringency"] = stringency

            slice.designs.append(design)
            designs.append(slice)

        return designs

    def _get_config_data(self) -> dict:
        try:
            config_data = parse_json(self._p3_config)
        except FileNotFoundError:
            raise
        except Exception:
            raise InvalidConfigError('Primer3 config file is not a correct JSON')

        return config_data
