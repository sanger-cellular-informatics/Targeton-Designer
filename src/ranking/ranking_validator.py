from venv import logger
from typing import List

from custom_logger.custom_logger import CustomLogger
from ranking.rank_criteria import RankingCriteria

# Initialize logger
logger = CustomLogger(__name__)

class RankCriteriaValidator:
    def __init__(self, config: dict):
        self._ranking_columns = config.get("ranking_criteria", [])
        self.created_criteria: List[RankingCriteria] = []
        self.is_validated: bool

        if not self._ranking_columns:
            logger.info("No Ranking criteria provided.")
            return


    def validated_criteria(self) -> List[RankingCriteria]:
        criteria = []

        for criterion_name in self._ranking_columns:
            rank_criterion = _get_rank_criterion(criterion_name)

            if rank_criterion in criteria:
                raise ValueError(f'Repeated ranking criteria: the given rank criterion "{criterion_name}" is repeated')

            criteria.append(rank_criterion)

        return criteria


def _get_rank_criterion(name: str) -> RankingCriteria:
    _check_criterion_name(name)

    return RankingCriteria.get_criterion_by_name(name)


def _check_criterion_name(name: str):
    if name not in [_criterion.name for _criterion in RankingCriteria.__subclasses__()]:
        raise ValueError(f'Invalid ranking criteria: the given rank criterion "{name}" is not valid')


#         for criteria in self.ranking_columns:
#             created = RankingCriteria.from_dict(criteria)
#             self.created_criteria.append(created)

#         validated, column_name, validated_rankers = validate_ranking_columns(self.created_criteria, self.csv_columns)

#         if not validated:
#             logger.warning(f"Invalid ranking column name {column_name}. Please Check configuration file.")
#             raise Exception(f"Invalid ranking column name {column_name}. Please Check configuration file.")
        
#         self.created_criteria = validated_rankers
    

# def validate_ranking_columns(created_criteria, csv_columns):

#     is_validated = True
#     column_name = ""
#     validated_rankers = []

#     for ranking_criteria in created_criteria:
#         if not ranking_criteria.column in csv_columns:
#             is_validated = not is_validated
#             column_name = ranking_criteria.column
#             return is_validated, column_name, validated_rankers
        
#         validated_rankers.append(ranking_criteria)

#     return is_validated, validated_rankers, validated_rankers