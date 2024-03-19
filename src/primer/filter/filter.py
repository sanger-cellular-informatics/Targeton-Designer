from abc import ABC, abstractmethod
from typing import List

from primer.filter.filter_response import FilterResponse
from primer.primer_pair import PrimerPair


class Filter(ABC):
    key: str
    value_type: type
    reason_discarded: str

    @abstractmethod
    def apply(self, mbs: List[PrimerPair]) -> FilterResponse:
        pass
