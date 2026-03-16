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
        design,
        slice_data: SliceData,
        stringency: float,
) -> List[PrimerPair]:
    primer_pairs = []

    primer_pairs_dict_list = design['PRIMER_PAIR']
    primer_left_dict_list = design['PRIMER_LEFT']
    primer_right_dict_list = design['PRIMER_RIGHT']

    for index,primer_pair2 in enumerate(primer_pairs_dict_list):
        pp = PrimerPair(
            pair_id= slice_data.name + "_LibAmp_" + str(index) + "_str" + str(stringency).replace(".", ""),
            chromosome=slice_data.chromosome,
            pre_targeton_start=slice_data.start,
            pre_targeton_end=slice_data.end,
            product_size=primer_pair2["PRODUCT_SIZE"],
            stringency=stringency,
            targeton_id=slice_data.targeton_id,
            uid = str(uuid.uuid1())
        )

        lefty_dict = primer_left_dict_list[index]

        start_left, end_left = calculate_primer_coords("left",lefty_dict["COORDS"],slice_data.start,slice_data.end,slice_data.strand)
        lefty = DesignedPrimer(
            name=f"{slice_data.name}_{name_primers('left', slice_data.strand)}_{index}",
            penalty=lefty_dict["PENALTY"],
            pair_id=pp.id,
            sequence=lefty_dict["SEQUENCE"],
            coords=Interval(start=lefty_dict["COORDS"][0], end=lefty_dict["COORDS"][1]),
            primer_start=start_left,
            primer_end=end_left,
            strand=determine_primer_strands("left", slice_data.strand),
            tm=lefty_dict["TM"],
            gc_percent=lefty_dict["GC_PERCENT"],
            self_any_th=lefty_dict["SELF_ANY_TH"],
            self_end_th=lefty_dict["SELF_END_TH"],
            hairpin_th=lefty_dict["HAIRPIN_TH"],
            end_stability=lefty_dict["END_STABILITY"]

        )
        righty_dict = primer_right_dict_list[index]

        start_right, end_right = calculate_primer_coords("right",righty_dict["COORDS"],slice_data.start,slice_data.end,slice_data.strand)

        righty = DesignedPrimer(
            name=f"{slice_data.name}_{name_primers('right', slice_data.strand)}_{index}",
            penalty=righty_dict["PENALTY"],
            pair_id=pp.id,
            sequence=righty_dict["SEQUENCE"],
            coords=Interval(start=righty_dict["COORDS"][0], end=righty_dict["COORDS"][1]),
            primer_start=start_right,
            primer_end=end_right,
            strand=determine_primer_strands("right", slice_data.strand),
            tm=righty_dict["TM"],
            gc_percent=righty_dict["GC_PERCENT"],
            self_any_th=righty_dict["SELF_ANY_TH"],
            self_end_th=righty_dict["SELF_END_TH"],
            hairpin_th=righty_dict["HAIRPIN_TH"],
            end_stability=righty_dict["END_STABILITY"]
        )


        if "LibAmpF" in lefty.name and "LibAmpR" in righty.name:
            pp.forward = lefty
            pp.reverse = righty
        else:
            pp.forward = righty
            pp.reverse = lefty

        primer_pairs.append(pp)

    return primer_pairs
