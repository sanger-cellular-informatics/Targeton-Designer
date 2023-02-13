from __future__ import annotations

from typing import Any, List
from collections import defaultdict
import re
import json
import numpy as np
from utils.file_system import read_csv_to_list_dict
from utils.exceptions import InputTypeError
from utils.write_output_files import DesignOutputData


class PrimerDesigner():
    def __init__(self, data=DesignOutputData('')) -> None:
        self.primer_pairs = []
        if self.validate_input(data):
            self.prepare_primer_designer(data)

    def get_primer_pairs(self) -> List[PrimerPair]:
        return self.primer_pairs

    def append_pair(self, pair) -> None:
        self.primer_pairs.append(pair)

    def get_fields(self) -> List[str]:
        keys = []
        for sub_dict in self.to_list_dicts():
            keys.extend(list(sub_dict.keys()))
        return list(set(keys))

    def dump_json(self, *args, **kwargs) -> None:
        json.dump(self.to_list_dicts(), *args, **kwargs)

    def to_list_dicts(self) -> List[defaultdict(dict)]:
        pair_dict_list = []
        for pair in self.get_primer_pairs():
            return_dict = pair._asdict()
            pair_dict_list.append(return_dict)
        return pair_dict_list

    def from_dict(self, data_dict) -> None:
        if isinstance(data_dict, list):
            for element in data_dict:
                self.append_pair(PrimerPair(element))
        elif isinstance(data_dict, dict):
            self.append_pair(PrimerPair(data_dict))
        else:
            raise InputTypeError("PrimerDesigner.from_dict expects a list of dicts or dict input.")

    def flatten(self) -> List[dict]:
        flat_dict_list = []
        for pair in self.get_primer_pairs():
            return_dict = pair._asdict()
            side_list = []
            for side in ('left', 'right'):
                side_return_dict = return_dict.copy()
                side_return_dict['side'] = side
                side_return_dict.update(side_return_dict[side])
                side_return_dict.pop('left')
                side_return_dict.pop('right')
                side_list.append(side_return_dict)
            flat_dict_list.extend(side_list)
        return flat_dict_list

    def copy(self) -> PrimerDesigner:
        new_primer_designer = PrimerDesigner()
        new_primer_designer.primer_pairs = [pair.copy() for pair in self.get_primer_pairs()]
        return new_primer_designer

    def prepare_primer_designer(self, design_output_data: DesignOutputData) -> None:
        if not self.validate_input(design_output_data):
            raise FileNotFoundError("Primer designer design output data not found in input.")

        primers = read_csv_to_list_dict(design_output_data.p3_csv)
        scoring = read_csv_to_list_dict(design_output_data.scoring_tsv, delimiter='\t')
        scoring = [score for score in scoring if score['A/B/Total'] == 'Total']
        pairs = iterate_design(primers, scoring)
        self.build_pair_classes(pairs)

    def build_pair_classes(self, pairs: defaultdict(dict)) -> None:
        for pair_key, pair in pairs.items():
            left = extract_primer_data(pair['F'])
            right = extract_primer_data(pair['R'])

            pair_class = PrimerPair({
                'pair' : pair_key,
                'score' : pair['score'],
                'left' : Primer(left),
                'right' : Primer(right),
            })

            self.append_pair(pair_class)

    @staticmethod
    def validate_input(data : DesignOutputData, needed_fields=["p3_csv", "scoring_tsv"]) -> bool:
        data_check = True
        for field in needed_fields:
            if not getattr(data, field):
                data_check = False
        return data_check


class PrimerPair():
    def __init__(self, data : dict) -> None:
        self.pair = data['pair']
        self.score = data['score']
        self.left = Primer(data['left'])
        self.right = Primer(data['right'])
        self.product_size = self.get_product_size()

    def get_paired_dict(self) -> dict:
        return vars(self)

    def get_fields(self) -> List[str]:
        return list(self._asdict().keys())

    def _asdict(self) -> dict:
        return_dict = vars(self)
        for k, v in return_dict.items():
            if isinstance(v, Primer):
                return_dict[k] = v._asdict()
        return return_dict

    def get_product_size(self) -> int:
        starts = [int(self.left.chr_start), int(self.right.chr_start)]
        ends = [int(self.left.chr_end), int(self.right.chr_end)]
        product_size = int(np.max(ends) - np.min(starts))
        return product_size

    def copy(self) -> PrimerPair:
        new_primer_pair = PrimerPair(self._asdict().copy())
        return new_primer_pair


class Primer():
    def __init__(self, primer_data: dict) -> None:
        self.chromosome = primer_data['chromosome']
        self.chr_start = primer_data['chr_start']
        self.chr_end = primer_data['chr_end']
        self.seq = primer_data['seq']
        self.melting_temp = primer_data['melting_temp']

    def __getitem__(self, item) -> Any:
        return getattr(self, item)

    def get_fields(self) -> List[str]:
        return list(self._asdict().keys())

    def _asdict(self) -> dict:
        return vars(self)


def extract_primer_data(primer_dict: dict) -> dict:
    record = {
        'chromosome' : primer_dict['chr'],
        'chr_start' : primer_dict['primer_start'],
        'chr_end' : primer_dict['primer_end'],
        'seq' : primer_dict['sequence'],
        'melting_temp' : primer_dict['tm'],
    }
    # For ipcress compatibility.
    del primer_dict['chr']

    return record


def iterate_design(primers: list, scoring: list) -> defaultdict:
    pairs = defaultdict(dict)
    for primer in primers:
        pairs = map_primer_data(primer, scoring, pairs)

    return pairs


def map_primer_data(primer: dict, scoring: list, pairs: defaultdict) -> defaultdict:
    pair = primer['primer']
    if (pair):
        chrom = primer['chr']
        # ENSE00000769557_HG8_16_LibAmpR_0
        # (ENSE00000769557_HG8_16_LibAmp)(R)_(0)
        match = re.search("^(\w+LibAmp)([F|R])\_(\d+)$", pair)
        pair_key = match.group(1) + '_' + match.group(3)
        side = match.group(2)
        pairs[pair_key][side] = primer
        pairs[pair_key][side]['chr'] = chrom
        pairs[pair_key]['score'] = [score['Score'] for score in scoring if score['Primer pair'] == pair_key][0]

    return pairs
