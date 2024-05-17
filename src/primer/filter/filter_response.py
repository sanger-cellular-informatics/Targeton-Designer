from dataclasses import dataclass
from typing import List

from primer.primer_pair import PrimerPair
from primer.primer_pair_discarded import PrimerPairDiscarded


@dataclass
class FilterResponse:
    primer_pairs_to_keep: List[PrimerPair]
    primer_pairs_to_discard: List[PrimerPairDiscarded]
