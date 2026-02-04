from typing import Tuple, List, Optional
from collections import defaultdict
import re
import uuid

from primer.designed_primer import map_to_designed_primer
from utils.get_data.hap1 import contain_variant
from primer.slice_data import SliceData


class PrimerPair:
    def __init__(self, pair_id: str, 
                       chromosome: str,
                       pre_targeton_start: int,
                       pre_targeton_end: int,
                       product_size: int,
                       stringency: float,
                       targeton_id: str,
                       uid: str):
        self.id = pair_id
        self.uid = uid
        self.chromosome = chromosome
        self.pre_targeton_start = pre_targeton_start
        self.pre_targeton_end = pre_targeton_end
        self.product_size = product_size
        self.stringency = stringency
        self.targeton_id = targeton_id
        self.forward_primer_data = {}
        self.reverse_primer_data = {}
        self.reverse = None
        self.forward = None


    def __repr__(self):
        return (f"PrimerPair(pair_id='{self.id}', "
                f"uid='{self.uid}', "
                f"chromosome='{self.chromosome}', "
                f"pre_targeton_start='{self.pre_targeton_start}', "
                f"pre_targeton_end='{self.pre_targeton_end}', "
                f"product_size='{self.product_size}', "
                f"stringency='{self.stringency}',"
                f"targeton_id='{self.targeton_id}', "
                f"forward={self.forward}, "
                f"reverse={self.reverse})"
                )

    def __eq__(self, other):
        if isinstance(other, PrimerPair):
            return (
                    self.chromosome == other.chromosome and
                    self.forward == other.forward and
                    self.reverse == other.reverse
            )
        return False

    def __hash__(self):
        return hash((self.chromosome, self.forward, self.reverse))

    @property
    def contain_hap_one_variant(self) -> bool:
        forward_start, forward_end = self.forward.primer_start, self.forward.primer_end
        reverse_start, reverse_end = self.reverse.primer_start, self.reverse.primer_end

        return (contain_variant(self.chromosome, forward_start, forward_end) or
                contain_variant(self.chromosome, reverse_start, reverse_end))


def build_primer_loci(
        primer,
        key,
        design,
        primer_details,
        slice_data: SliceData,
        primer_name: str,
        primer_pair_id: str,
) -> dict:
    primer_field = primer_details['field']

    primer['primer'] = primer_name
    primer[primer_field] = design[key]

    primer['side'] = primer_details['side']
    primer['pair_id'] = primer_pair_id

    if primer_field == 'coords':
        primer_coords = calculate_primer_coords(
            primer_details['side'],
            design[key],
            slice_data.start,
            slice_data.end,
            slice_data.strand
        )

        primer['primer_start'] = primer_coords[0]
        primer['primer_end'] = primer_coords[1]
        primer['strand'] = determine_primer_strands(
            primer_details['side'], slice_data.strand)

    return primer


def name_primers(side: str, strand: str) -> str:
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

    primer_name = names[strand][side]

    return primer_name


def capture_primer_details(primer_name: str) -> dict:
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


def calculate_primer_coords(side: str, coords: list,
                            slice_start: int, slice_end: int,
                            strand: str) -> Tuple[int, int]:
    if strand == "+":
        left_flank = {
            'start': slice_start + int(coords[0]),
            'end': slice_start + int(coords[0]) + int(coords[1]) - 1
        }

        right_end = slice_start + int(coords[0])
        right_flank = {
            'start': 1 + right_end - int(coords[1]),
            'end': right_end,
        }

    if strand == "-":
        left_flank = {
            'start': slice_end - int(coords[0]) - int(coords[1]) + 1,
            'end': slice_end - int(coords[0])
        }

        right_start = slice_end - int(coords[0])
        right_flank = {
            'start': right_start,
            'end': right_start + coords[1] - 1,
        }

    slice_coords = {
        'left': left_flank,
        'right': right_flank
    }

    start = slice_coords[side]['start']
    end = slice_coords[side]['end']

    return start, end


def determine_primer_strands(side: str, slice_strand: str) -> str:
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


def build_primer_pairs(
        design,
        slice_data: SliceData,
        stringency: float,
) -> List[PrimerPair]:
    primer_pairs = []
    primers = defaultdict(dict)

    for key in design.keys():
        primer_details = capture_primer_details(key)

        if primer_details:
            libamp_name = name_primers(primer_details['side'], slice_data.strand)
            primer_name = slice_data.name + "_" + libamp_name + "_" + \
                          primer_details['pair']

            stringency_string = "_str" + str(stringency).replace(".", "")

            primer_name_with_stringency = primer_name + stringency_string
            primer_pair_id = slice_data.name + "_LibAmp_" + primer_details['pair'] + stringency_string

            primer_pair_product_size = design['PRIMER_PAIR_' + primer_details['pair'] + '_PRODUCT_SIZE']

            primer = build_primer_loci(
                primers[primer_name_with_stringency],
                key,
                design,
                primer_details,
                slice_data,
                primer_name,
                primer_pair_id,
            )

            pair = _find_pair_by_id(primer_pairs, primer_pair_id)
            if pair is None:
                uid = str(uuid.uuid1())
                pair = PrimerPair(
                    pair_id=primer_pair_id,
                    chromosome=slice_data.chromosome,
                    pre_targeton_start=slice_data.start,
                    pre_targeton_end=slice_data.end,
                    product_size=primer_pair_product_size,
                    stringency=stringency,
                    targeton_id=slice_data.targeton_id,
                    uid=uid
                )
                primer_pairs.append(pair)

            if libamp_name == "LibAmpF":
                pair.forward_primer_data = primer
            if libamp_name == "LibAmpR":
                pair.reverse_primer_data = primer

    primer_pairs = _map_primers_into_designed_primers_objects(primer_pairs)
    return primer_pairs


def _map_primers_into_designed_primers_objects(primer_pairs: List[PrimerPair]) -> List[PrimerPair]:
    for pair in primer_pairs:
        pair.forward = map_to_designed_primer(pair.forward_primer_data)
        pair.reverse = map_to_designed_primer(pair.reverse_primer_data)

    return primer_pairs


def _find_pair_by_id(pairs: List[PrimerPair], pair_id: str) -> Optional[PrimerPair]:
    for pair in pairs:
        if pair.id == pair_id:
            return pair
    return None
