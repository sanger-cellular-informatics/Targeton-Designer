from collections import defaultdict
from typing import List, Tuple

from primer.filter.filter import Filter
from primer.filter.filter_response import FilterResponse, PrimerPairDiscarded
from primer.primer_pair import PrimerPair


class DuplicatesFilter(Filter):
    key: str = 'duplicates'
    value_type: type = bool
    reason_discarded: str = "has duplicate with a higher stringency"

    def apply(self, pairs: List[PrimerPair]) -> FilterResponse:
        pairs_to_keep = []
        pairs_to_discard = []

        pair_duplicates_grouped = _group_duplicate_pairs(pairs)

        for duplicate_group in pair_duplicates_grouped:
            pair_max_stringency, others = _take_pair_with_max_stringency_and_others(duplicate_group)

            pairs_from_group_to_discard = [PrimerPairDiscarded(pair, DuplicatesFilter.reason_discarded) for pair in others]

            pairs_to_keep.append(pair_max_stringency)
            pairs_to_discard.extend(pairs_from_group_to_discard)

        return FilterResponse(pairs_to_keep, pairs_to_discard)


def _take_pair_with_max_stringency_and_others(pairs: List[PrimerPair]) -> Tuple[PrimerPair, List[PrimerPair]]:
    pair_max_stringency = _get_max_stringency_pair(pairs)
    other_pairs = [pair for pair in pairs if pair.stringency != pair_max_stringency.stringency]

    return pair_max_stringency, other_pairs


def _get_max_stringency_pair(duplicate_group: List[PrimerPair]) -> PrimerPair:
    return min(duplicate_group, key=lambda pair: pair.stringency)


def _group_duplicate_pairs(pairs: List[PrimerPair]) -> List[List[PrimerPair]]:
    groups_duplicates = defaultdict(list)
    for pair in pairs:
        groups_duplicates[pair].append(pair)
    return list(groups_duplicates.values())
