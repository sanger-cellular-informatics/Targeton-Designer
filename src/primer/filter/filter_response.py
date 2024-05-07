from dataclasses import dataclass
from typing import List

from primer.primer_pair import PrimerPair


@dataclass
class PrimerPairDiscarded:
    primer_pair: PrimerPair
    reason_discarded: str


@dataclass
class FilterResponse:
    primer_pairs_to_keep: List[PrimerPair]
    primer_pairs_to_discard: List[PrimerPairDiscarded]
