from typing import List
from itertools import chain

from primer.filter.filter import Filter
from primer.filter.hap1_variant_filter import HAP1VariantFilter
from primer.filter.duplicates_filter import DuplicatesFilter
from primer.filter.filter_response import FilterResponse
from primer.primer_pair import PrimerPair

from custom_logger.custom_logger import CustomLogger

# Initialize logger
logger = CustomLogger(__name__)


class FilterManager:
    def __init__(self, apply_filters: list):
      
        self.filters: List[Filter] = [DuplicatesFilter(), HAP1VariantFilter()]

        self._filters_to_apply: List[Filter] = []

        # Even if user forgot or do not add "duplicates" filter in config it will be applied.
        if not apply_filters and not "duplicates" in apply_filters:
            apply_filters.append("duplicates")

        for is_filter in self.filters:
            
            # Check if user added incorrect filter name.
            if not is_filter.key in apply_filters:
                incorrect_filter_names = [filter_name for filter_name in apply_filters if is_filter.key != filter_name]
                [logger.info(f"Incorrect filter name {incorrect_names}.") for incorrect_names in incorrect_filter_names]
                 
            
            if is_filter.key in apply_filters:
                self._filters_to_apply.append(is_filter)        
        

    def apply_filters(self, data: List[PrimerPair]) -> FilterResponse:
        primer_pairs = data

        pairs_to_keep = []
        pairs_to_discard = []

        for _filter in self._filters_to_apply:

            logger.info(f"Filter {_filter.key} is applied.")

            filter_response = _filter.apply(primer_pairs)

            pairs_to_keep.extend(filter_response.primer_pairs_to_keep)
            pairs_to_discard.extend(filter_response.primer_pairs_to_discard)
     
        return FilterResponse(primer_pairs_to_keep=pairs_to_keep, primer_pairs_to_discard=pairs_to_discard)
