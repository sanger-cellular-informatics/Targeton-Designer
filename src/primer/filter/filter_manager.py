from typing import List

from primer.filter.filter import Filter
from primer.filter.hap1_variant_filter import HAP1VariantFilter
from primer.filter.duplicates_filter import DuplicatesFilter
from primer.filter.filter_response import FilterResponse
from primer.primer_pair import PrimerPair

from custom_logger.custom_logger import CustomLogger

# Initialize logger
logger = CustomLogger(__name__)


class FilterManager:
    def __init__(self, filters: dict):
        
        self._filters_to_apply: List[Filter] = [DuplicatesFilter()]

        if not filters["hap1"]:
            logger.info("HAP1 filter is not applied. Use --filters hap1 in command arguments to apply HAP1.")

        if filters["hap1"]:
            logger.info("HAP1 filter is applied")
            self._filters_to_apply: List[Filter] = [DuplicatesFilter(), HAP1VariantFilter()]

    def apply_filters(self, data: List[PrimerPair]) -> FilterResponse:
        pairs_to_keep = data
        pairs_to_discard = []
        for _filter in self._filters_to_apply:
            filter_response = _filter.apply(pairs_to_keep)

            pairs_to_keep = filter_response.primer_pairs_to_keep
            pairs_to_discard += filter_response.primer_pairs_to_discard
     
        return FilterResponse(primer_pairs_to_keep=pairs_to_keep, primer_pairs_to_discard=pairs_to_discard)
