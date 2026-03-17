from dataclasses import dataclass
from enum import Enum


@dataclass
class Interval:
    start: int
    end: int

class Orientation(Enum):
    FORWARD = "F"
    REVERSE = "R"

    @classmethod
    def from_side_and_strand(cls, side: str, strand: str) -> "Orientation":
        mapping = {
            ("+", "left"): cls.FORWARD,
            ("+", "right"): cls.REVERSE,
            ("-", "left"): cls.REVERSE,
            ("-", "right"): cls.FORWARD,
        }

        try:
            return mapping[(strand, side)]
        except KeyError:
            raise ValueError(f"Invalid combination: strand={strand}, side={side}")

class Strand(Enum):
    NEGATIVE = "-"
    POSITIVE = "+"

    @classmethod
    def from_side_and_slice_strand(cls, side: str, slice_strand: str) -> "Strand":
        mapping = {
            ("+", "left"): cls.POSITIVE,
            ("+", "right"): cls.NEGATIVE,
            ("-", "left"): cls.NEGATIVE,
            ("-", "right"): cls.POSITIVE,
        }

        try:
            return mapping[(slice_strand, side)]
        except KeyError:
            raise ValueError(f"Invalid combination: slice_strand={slice_strand}, side={side}")

@dataclass
class DesignedPrimer:
    name: str
    penalty: float
    pair_id: str
    sequence: str
    coords: Interval
    primer_start: int
    primer_end: int
    strand: Strand
    tm: float
    gc_percent: float
    self_any_th: float
    self_end_th: float
    hairpin_th: float
    end_stability: float
    orientation: Orientation

    def __eq__(self, other):
        if isinstance(other, DesignedPrimer):
            # Exclude the 'name' and 'pair_id' attributes from the comparison
            return (
                    self.penalty == other.penalty and
                    self.sequence == other.sequence and
                    self.coords == other.coords and
                    self.primer_start == other.primer_start and
                    self.primer_end == other.primer_end and
                    self.strand == other.strand and
                    self.orientation == other.orientation and
                    self.tm == other.tm and
                    self.gc_percent == other.gc_percent and
                    self.self_any_th == other.self_any_th and
                    self.self_end_th == other.self_end_th and
                    self.hairpin_th == other.hairpin_th and
                    self.end_stability == other.end_stability
            )
        return False

    def __hash__(self):
        # Exclude the 'name' and 'pair_id' attributes from the comparison
        return hash((
            self.penalty,
            self.sequence,
            (self.coords.start, self.coords.end),
            self.primer_start,
            self.primer_end,
            self.strand,
            self.tm,
            self.gc_percent,
            self.self_any_th,
            self.self_end_th,
            self.hairpin_th,
            self.end_stability
        ))
