from typing import List

from primer.filter.filter import Filter
from primer.filter.filter_response import FilterResponse
from primer.primer_pair import PrimerPair
from primer.primer_pair_discarded import PrimerPairDiscarded


class HAP1VariantFilter(Filter):
    key: str = 'HAP1_variant'
    value_type: type = bool
    reason_discarded: str = "contains background variants"

    def apply(self, pairs: List[PrimerPair]) -> FilterResponse:
        pairs_to_keep = []
        pairs_to_discard = []

        for pair in pairs:
            if pair.contain_hap_one_variant:
                pairs_to_discard.append(PrimerPairDiscarded(pair, HAP1VariantFilter.reason_discarded))
            else:
                pairs_to_keep.append(pair)

        return FilterResponse(pairs_to_keep, pairs_to_discard)
