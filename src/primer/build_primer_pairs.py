from primer.designed_primer import DesignedPrimer, Interval, Orientation, Strand
from primer.primer_pair import PrimerPair
from primer.slice_data import SliceData

from typing import Tuple, List, Optional
import uuid


def build_primer_pairs(
        design: dict,
        slice_data: SliceData,
        stringency: float,
        primer_type: str
) -> List[PrimerPair]:
    stringency_tag = str(stringency).replace(".", "")

    primer_pairs = []
    zipped = zip(design["PRIMER_PAIR"], design["PRIMER_LEFT"], design["PRIMER_RIGHT"], )

    for index, (pair_dict, left_dict, right_dict) in enumerate(zipped):
        primer_pair = PrimerPair(
            pair_id=f"{slice_data.name}_{primer_type}_{index}_str{stringency_tag}",
            chromosome=slice_data.chromosome,
            pre_targeton_start=slice_data.start,
            pre_targeton_end=slice_data.end,
            product_size=pair_dict["PRODUCT_SIZE"],
            stringency=stringency,
            targeton_id=slice_data.targeton_id,
            uid=str(uuid.uuid1()),
        )

        left_primer = _build_designed_primer("left", left_dict, slice_data, primer_pair.id, index, primer_type)
        right_primer = _build_designed_primer ("right", right_dict, slice_data, primer_pair.id, index, primer_type)

        _assign_forward_reverse2(primer_pair, left_primer, right_primer)
        primer_pairs.append(primer_pair)

    return primer_pairs


def _build_designed_primer(
        side: str,
        primer_dict: dict,
        slice: SliceData,
        pair_id: str,
        index: int,
        primer_type: str
) -> DesignedPrimer:
    start, end = _calculate_primer_coords(side, primer_dict["COORDS"], slice.start, slice.end, slice.strand)
    orientation = Orientation.from_side_and_strand(side=side, strand=slice.strand)

    return DesignedPrimer(
        name=f"{slice.name}_{primer_type}{orientation.value}_{index}",
        penalty=primer_dict["PENALTY"],
        pair_id=pair_id,
        sequence=primer_dict["SEQUENCE"],
        coords=Interval(start=primer_dict["COORDS"][0], end=primer_dict["COORDS"][1]),
        primer_start=start,
        primer_end=end,
        strand=Strand.from_side_and_slice_strand(side=side, slice_strand=slice.strand),
        tm=primer_dict["TM"],
        gc_percent=primer_dict["GC_PERCENT"],
        self_any_th=primer_dict["SELF_ANY_TH"],
        self_end_th=primer_dict["SELF_END_TH"],
        hairpin_th=primer_dict["HAIRPIN_TH"],
        end_stability=primer_dict["END_STABILITY"],
        orientation=orientation,
    )


def _name_primers(side: str, strand: str) -> str:
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


def _calculate_primer_coords(side: str, coords: list,
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


def _determine_primer_strands(side: str, slice_strand: str) -> str:
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


def _assign_forward_reverse(pp: PrimerPair, left: DesignedPrimer, right: DesignedPrimer) -> None:
    if "LibAmpF" in left.name and "LibAmpR" in right.name:
        pp.forward, pp.reverse = left, right
    else:
        pp.forward, pp.reverse = right, left

def _assign_forward_reverse2(pp: PrimerPair, left: DesignedPrimer, right: DesignedPrimer) -> None:
    if left.orientation == Orientation.FORWARD and right.orientation == Orientation.REVERSE:
        pp.forward, pp.reverse = left, right
    elif left.orientation == Orientation.REVERSE and right.orientation == Orientation.FORWARD:
        pp.forward, pp.reverse = right, left
    else:
        raise ValueError(
            f"Cannot assign forward/reverse: "
            f"left.orientation={left.orientation}, right.orientation={right.orientation}"
        )
