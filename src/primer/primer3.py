import primer3
import re
from collections import defaultdict
from _collections_abc import dict_keys
from typing import Tuple, List

from Bio.Seq import Seq

from utils.exceptions import InvalidConfigError
from utils.file_system import parse_json

from primer.slice_data import SliceData
from primer.primer3_prepare_config import prepare_config


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
        slices_dict = self._locate_primers(designs)

        return slices_dict

    def _primer3_run(self, slices: List[SliceData]) -> List[SliceData]:
       # STRINGENCY_VECTOR = ["0.1", "0.25", "0.5", "0.75", "1.0"]
        STRINGENCY_VECTOR = ["0.1", "1.0"]

        base_config = self._get_config_data()

        for stringency in STRINGENCY_VECTOR:
            config_data = prepare_config(base_config, stringency)

            print(config_data)
            result = self._primer3_design(slices, config_data)

        return result

    def _primer3_design(self, slices: List[SliceData], config: dict) -> List[SliceData]:
        designs = []
        for slice in slices:
            primer3_input = slice.p3_input
            design = primer3.bindings.design_primers(primer3_input, config)

            stringency = config["PRIMER_MASK_FAILURE_RATE"]
            design["stringency"] = stringency

            slice.designs.append(design)
            designs.append(slice)

        return designs

    def _locate_primers(self, slices: List[SliceData]) -> List[dict]:
        slice_designs = []

        for slice_data in slices:
            slice_data.primers = {}
            for design in slice_data.designs:

                primer_keys = design.keys()

                primers = self._build_primers_dict(design, primer_keys, slice_data, design['stringency'])

                for primer in primers:
                    slice_data.primers[primer] = primers[primer]
                    slice_designs.append(slice_data)

            del slice_data.designs

        return slice_designs

    def _build_primers_dict(
            self,
            design,
            primer_keys: dict_keys,
            slice_data: dict,
            stringency: str,
    ) -> defaultdict(dict):
        primers = defaultdict(dict)

        for key in primer_keys:
            primer_details = self.capture_primer_details(key)

            if primer_details:
                libamp_name = self.name_primers(primer_details, slice_data.strand)
                primer_name = slice_data.name + "_" + libamp_name + "_" + primer_details['pair']

                primer_name_with_stringency = primer_name + "_str" + stringency.replace(".", "_")
                primer_pair_id = slice_data.name + "_" + primer_details['pair'] + "_str" + stringency.replace(".", "_")

                primers[primer_name_with_stringency] = self._build_primer_loci(
                    primers[primer_name_with_stringency],
                    key,
                    design,
                    primer_details,
                    slice_data,
                    stringency,
                    primer_name,
                    primer_pair_id,
                )

        return primers

    def _build_primer_loci(
        self,
        primer,
        key,
        design,
        primer_details,
        slice_data: SliceData,
        stringency: str,
        primer_name: str,
        primer_pair_id: str,
    ) -> dict:

        primer_field = primer_details['field']

        primer['name'] = primer_name
        primer[primer_field] = design[key]
        
        primer['side'] = primer_details['side']
        primer['stringency'] = stringency
        primer['pair_id'] = primer_pair_id

        if primer_field == 'coords':
            primer_coords = self.calculate_primer_coords(
                primer_details['side'], design[key], slice_data.start)

            primer['primer_start'] = primer_coords[0]
            primer['primer_end'] = primer_coords[1]
            primer['strand'] = self.determine_primer_strands(
                primer_details['side'], slice_data.strand
            )

        return primer

    def _get_config_data(self) -> dict:
        try:
            config_data = parse_json(self._config)
        except FileNotFoundError:
            raise
        except Exception:
            raise InvalidConfigError('Primer3 config file is not a correct JSON')

        return config_data

    @staticmethod
    def name_primers(primer_details, strand) -> str:
        fwd_primers = {
            'left': 'LibAmpF',
            'right': 'LibAmpR',
        }
        rev_primers = {
            'left': 'LibAmpR',
            'right': 'LibAmpF',
        }
        names = {
            '+': fwd_primers,
            '-': rev_primers,
        }

        primer_name = names[strand][primer_details['side']]

        return primer_name

    @staticmethod
    def capture_primer_details(primer_name) -> dict:
        match = re.search(r'^(primer_(left|right)_(\d+))(\_(\S+))?$', primer_name.lower())
        result = {}
        if match:
            primer_id = match.group(1)
            primer_side = match.group(2)
            pair_number = match.group(3)
            primer_field = match.group(5)
            if primer_field is None:
                primer_field = 'coords'
            result = {
                'id': primer_id,
                'side': primer_side,
                'field': primer_field,
                'pair': pair_number
            }

        return result

    @staticmethod
    def calculate_primer_coords(side, coords, slice_start) -> Tuple[int, int]:
        slice_start = int(slice_start)
        left_flank = {
            'start': slice_start,
            'end': slice_start + int(coords[1])
        }

        slice_end = slice_start + int(coords[0])
        right_flank = {
            'start': 1 + slice_end - int(coords[1]),
            'end': 1 + slice_end,
        }

        slice_coords = {
            'left': left_flank,
            'right': right_flank
        }

        start = slice_coords[side]['start']
        end = slice_coords[side]['end']

        return start, end

    @staticmethod
    def determine_primer_strands(side, slice_strand) -> str:
        positive = {
            'left': '+',
            'right': '-',
        }

        negative = {
            'left': '-',
            'right': '+',
        }

        strands = {
            '+': positive,
            '-': negative,
        }

        return strands[slice_strand][side]

    @staticmethod
    def revcom_reverse_primer(seq, strand) -> Seq:
        seq_obj = Seq(seq)

        if strand == '-':
            seq_obj = seq_obj.reverse_complement()

        return seq_obj
