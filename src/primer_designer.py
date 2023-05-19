from __future__ import annotations

from typing import Any, List, Union
from collections import defaultdict
import re
import json
import numpy as np
from utils.file_system import read_csv_to_list_dict
from utils.exceptions import InputTypeError
from utils.write_output_files import DesignOutputData

VERSION = '01'


class PrimerDesigner():
    def __init__(self, data: Union[DesignOutputData, dict, List[dict]] = dict()) -> None:
        self.primer_pairs = []
        if data:
            self.process_input(data)

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

    def from_design_output(self, design_output_data: DesignOutputData) -> None:
        primers = read_csv_to_list_dict(design_output_data.p3_csv)
        scoring = read_csv_to_list_dict(design_output_data.scoring_tsv, delimiter='\t')
        scoring = [score for score in scoring if score['A/B/Total'] == 'Total']
        pairs = iterate_design(primers, scoring)
        self.build_pair_classes(pairs)

    def build_pair_classes(self, pairs: defaultdict(dict)) -> None:
        for _, pair in pairs.items():
            pair_class = PrimerPair(pair)
            self.append_pair(pair_class)

    def process_input(self, data: Union[DesignOutputData, dict, List[dict]]) -> bool:
        if isinstance(data, (dict, list)):
            data_check = self.validate_dict_input(data)
            if data_check:
                self.from_dict(data)
            else:
                raise InputTypeError("PrimerDesigner class input dict (or list of dicts) doesn't contain needed fields.")
        elif data.__class__.__name__ == DesignOutputData.__name__:
            data_check = self.validate_file_input(data)
            if data_check:
                self.from_design_output(data)
            else:
                raise InputTypeError("PrimerDesigner class input DesignOutputData doesn't contain needed fields.")
        else:
            raise InputTypeError(f"PrimerDesigner class expects data input to be DesignOutputData or dict (or list of dicts), not: {type(data)}")

    @staticmethod
    def validate_file_input(data: DesignOutputData, needed_fields=["p3_csv", "scoring_tsv"]) -> bool:
        data_check = True
        for field in needed_fields:
            if not getattr(data, field):
                data_check = False

        return data_check

    def validate_dict_input(
            self,
            data : Union[dict,List[dict]],
            expected_structure={
                'pair'  : str,
                'score' : str,
                'left' : dict,
                'right' : dict,
                'product_size' : int,
                'targeton' : str,
                'version': str,
            }) -> bool:
        data_check = True
        if isinstance(data, list):
            for element in data:
                self.validate_dict_input(element)
        elif isinstance(data, dict):
            for key, type in expected_structure.items():
                if key in data:
                    if isinstance(data[key], type):
                        if key == 'left' or key == 'right':
                            expected_primer_structure = {
                                'chromosome' : str,
                                'chr_start' : str,
                                'chr_end' : str,
                                'seq' : str,
                                'melting_temp' : str,
                                'gc_content': str,
                            }
                            data_check = self.validate_dict_input(
                                data[key],
                                expected_structure=expected_primer_structure
                                )
                    else:
                        return False
                else:
                    return False
        else:
            return False

        return data_check


class PrimerPair():
    def __init__(self, data : dict) -> None:
        fields = [
            'left',
            'right',
            'score',
            'targeton',
            'version',
            'pair',
        ]
        if data:
            self.assign_data(data, fields=fields)

    def assign_data(self, data: dict, fields=[]) -> None:
        translation_dict = {
            'left' : 'F',
            'right' : 'R',
        }
        data = translate_dict(data, translation_dict=translation_dict, fields=fields)
        for k, v in data.items():
            if k in ['left', 'right']:
                self.__setattr__(k, Primer(v))
            else:
                self.__setattr__(k, v)
        self.product_size = self.get_product_size()

    def get_paired_dict(self) -> dict:
        return self._asdict()

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
        fields = [
            'chromosome',
            'chr_start',
            'chr_end',
            'seq',
            'melting_temp',
            'gc_content',
        ]
        if primer_data:
            self.assign_data(primer_data, fields=fields)

    def __getitem__(self, item) -> Any:
        return getattr(self, item)

    def get_fields(self) -> List[str]:
        return list(self._asdict().keys())

    def _asdict(self) -> dict:
        return vars(self)

    def assign_data(self, primer_dict: dict, fields=[]):
        translation_dict = {
            'chromosome' : 'chr',
            'chr_start' : 'primer_start',
            'chr_end' : 'primer_end',
            'seq' : 'sequence',
            'melting_temp' : 'tm',
            'gc_content' : 'gc_percent',
        }
        primer_dict = translate_dict(primer_dict, translation_dict=translation_dict, fields=fields)
        for k, v in primer_dict.items():
            self.__setattr__(k, v)


def translate_dict(data_dict: dict, translation_dict={}, fields=[]):
    translated_dict = {}
    keys = list(data_dict.keys())
    for new_key, old_key in translation_dict.items():
        if old_key in data_dict:
            translated_dict[new_key] = data_dict[old_key]
            keys.remove(old_key)
        elif new_key in data_dict:
            translated_dict[new_key] = data_dict[new_key]
            keys.remove(new_key)
    if not fields:
        fields = keys
    leftover_keys_dict = {k: data_dict[k] for k in keys if k in fields}
    translated_dict.update(leftover_keys_dict)
    return translated_dict


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
        match = re.search(r"^(\w+LibAmp)([F|R])\_(\d+)$", pair)
        pair_key = match.group(1) + '_' + match.group(3)
        side = match.group(2)
        pairs[pair_key][side] = primer
        pairs[pair_key][side]['chr'] = chrom
        pairs[pair_key]['version'] = VERSION
        pairs[pair_key]['pair'] = pair_key
        for score in scoring:
            if score['Primer pair'] == pair_key:
                pairs[pair_key]['score'] = score['Score']
                pairs[pair_key]['targeton'] = score['Targeton']

    return pairs
