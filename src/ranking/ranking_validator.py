from typing import List

from .ranking_criteria import create_criteria


class RankCriteriaValidator:
    def __init__(self, config: Config):
        self.ranking_criteria = config.ranking_criteria

    def validated_criteria(self) -> List[RankCriteria]:
        criteria = []

        for criterion_name in self._ranking_order:
            rank_criterion = _get_rank_criterion(criterion_name)

            if rank_criterion in criteria:
                raise ValueError(f'Repeated ranking criteria: the given rank criterion "{criterion_name}" is repeated')

            criteria.append(rank_criterion)

        return criteria


def _get_rank_criterion(name: str) -> RankCriteria:
    _check_criterion_name(name)

    return RankCriteria.get_criterion_by_name(name)


def _check_criterion_name(name: str):
    if name not in [_criterion.name for _criterion in RankCriteria.__subclasses__()]:
        raise ValueError(f'Invalid ranking criteria: the given rank criterion "{name}" is not valid')
