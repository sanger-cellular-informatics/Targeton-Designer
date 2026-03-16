from typing import Tuple, List, Optional
from collections import defaultdict
import re
import uuid

from primer.designed_primer import DesignedPrimer, Interval
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
        design: dict,
        slice_data: SliceData,
        stringency: float,
) -> List[PrimerPair]:
    stringency_tag = str(stringency).replace(".", "")

    primer_pairs = []
    zipped = zip(design["PRIMER_PAIR"], design["PRIMER_LEFT"], design["PRIMER_RIGHT"], )

    for index, (pair_dict, left_dict, right_dict) in enumerate(zipped):
        primer_pair = PrimerPair(
            pair_id=f"{slice_data.name}_LibAmp_{index}_str{stringency_tag}",
            chromosome=slice_data.chromosome,
            pre_targeton_start=slice_data.start,
            pre_targeton_end=slice_data.end,
            product_size=pair_dict["PRODUCT_SIZE"],
            stringency=stringency,
            targeton_id=slice_data.targeton_id,
            uid=str(uuid.uuid1()),
        )

        left_primer = _build_designed_primer("left", left_dict, slice_data, primer_pair.id, index)
        right_primer = _build_designed_primer("right", right_dict, slice_data, primer_pair.id, index)

        _assign_forward_reverse(primer_pair, left_primer, right_primer)
        primer_pairs.append(primer_pair)

    return primer_pairs

def _build_designed_primer(
        side: str,
        primer_dict: dict,
        slice: SliceData,
        pair_id: str,
        index: int,
) -> DesignedPrimer:
    start, end = calculate_primer_coords(side, primer_dict["COORDS"], slice.start, slice.end, slice.strand)

    return DesignedPrimer(
        name=f"{slice.name}_{name_primers(side, slice.strand)}_{index}",
        penalty=primer_dict["PENALTY"],
        pair_id=pair_id,
        sequence=primer_dict["SEQUENCE"],
        coords=Interval(start=primer_dict["COORDS"][0], end=primer_dict["COORDS"][1]),
        primer_start=start,
        primer_end=end,
        strand=determine_primer_strands(side, slice.strand),
        tm=primer_dict["TM"],
        gc_percent=primer_dict["GC_PERCENT"],
        self_any_th=primer_dict["SELF_ANY_TH"],
        self_end_th=primer_dict["SELF_END_TH"],
        hairpin_th=primer_dict["HAIRPIN_TH"],
        end_stability=primer_dict["END_STABILITY"],
    )

def _assign_forward_reverse(pp: PrimerPair, left: DesignedPrimer, right: DesignedPrimer) -> None:
    if "LibAmpF" in left.name and "LibAmpR" in right.name:
        pp.forward, pp.reverse = left, right
    else:
        pp.forward, pp.reverse = right, left
