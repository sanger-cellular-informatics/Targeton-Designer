from typing import List

from primer.filter.filter import Filter
from primer.filter.filter_response import FilterResponse, PrimerPairDiscarded
from primer.primer_pair import PrimerPair


class DuplicatesFilter(Filter):
    key: str = 'duplicates'
    value_type: type = bool
    reason_discarded: str = "has duplicated with a higher stringency"

    def apply(self, pairs: List[PrimerPair]) -> FilterResponse:
        pairs_to_keep = []
        pairs_to_discard = []

        # TODO
        # for pair in pairs:
        #     if CONDITION:
        #         pairs_to_discard.append(PrimerPairDiscarded(pair, DuplicatesFilter.reason_discarded))
        #     else:
        #         pairs_to_keep.append(pair)

        return FilterResponse(pairs_to_keep, pairs_to_discard)
