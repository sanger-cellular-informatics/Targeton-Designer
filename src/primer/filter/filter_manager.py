from collections import defaultdict
from inspect import getmembers, isclass, isabstract
from typing import List

import primer
from primer.filter.filter import Filter
from primer.filter.duplicates_filter import DuplicatesFilter
from primer.filter.hap1_variant_filter import HAP1VariantFilter
from primer.filter.filter_response import FilterResponse
from primer_designer import PrimerPair


class FilterManager:
    def __init__(self, filter_names: List[str]):
        self._filters_to_apply: List[Filter] = FilterManager.load_filters(filter_names)

    def apply_filters(self, data: List[PrimerPair]) -> FilterResponse:
        pairs_to_keep = data
        pairs_to_discard = []

        for _filter in self._filters_to_apply:
            filter_response = _filter.apply(pairs_to_keep)

            pairs_to_keep = filter_response.primer_pairs_to_keep
            pairs_to_discard += filter_response.primer_pairs_to_discard

        return FilterResponse(primer_pairs_to_keep=pairs_to_keep, primer_pairs_to_discard=pairs_to_discard)

    @staticmethod
    def load_filters(filter_names: List[str]) -> List[Filter]:
        classes = getmembers(primer.filter, lambda member: isclass(member) and not isabstract(member))
        filters = {cls.key: cls for _, cls in classes if issubclass(cls, Filter)}

        loaded_filters = []
        for name in set(filter_names):
            filter_cls = filters.get(name)
            if filter_cls:
                loaded_filters.append(filter_cls())
            else:
                raise ValueError(f"Filter name '{name}' does not exist")

        return loaded_filters
