from typing import List
import pandas as pd

from ranking.ranking_validator import RankCriteriaValidator, RankingCriteria


class Ranker:
    def __init__(self, config: dict):
        self.validated_rankers: List[RankingCriteria] = RankCriteriaValidator(config).created_criteria
    

    def rank(self):
        pass
