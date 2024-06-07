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
    