from typing import List

from primer.filter.filter import Filter
from primer.filter.hap1_variant_filter import HAP1VariantFilter
from primer.filter.duplicates_filter import DuplicatesFilter
from primer.filter.filter_response import FilterResponse
from primer_designer import PrimerPair


class FilterManager:
    def __init__(self):
        self._filters_to_apply: List[Filter] = [HAP1VariantFilter(), DuplicatesFilter()]

    def apply_filters(self, data: List[PrimerPair]) -> FilterResponse:
        pairs_to_keep = data
        pairs_to_discard = []

        for _filter in self._filters_to_apply:
            filter_response = _filter.apply(pairs_to_keep)

            pairs_to_keep = filter_response.primer_pairs_to_keep
            pairs_to_discard += filter_response.primer_pairs_to_discard

        return FilterResponse(primer_pairs_to_keep=pairs_to_keep, primer_pairs_to_discard=pairs_to_discard)
