from dataclasses import dataclass


@dataclass
class RankingCriteria:
    name: str
    is_ascending: bool
    column: str

    @staticmethod
    def get_criterion_by_name(name: str):
        for rank_criteria_class in RankingCriteria.__subclasses__():
            if name == rank_criteria_class.name:
                return rank_criteria_class


class StringencyCriteria(RankingCriteria):
    name: str = 'stringency'
    is_ascending: bool = False
    column: str = 'stringency'

class ProductSizeCriteria(RankingCriteria):
    name: str = 'product_size'
    is_ascending: bool = False
    column: str = 'product_size'