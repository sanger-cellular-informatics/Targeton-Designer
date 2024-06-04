from typing import List
import pandas as pd

from ranking.ranking_validator import RankingCriteria


class Ranker:
    def __init__(self, validated_rankers: List[RankingCriteria]):
        print(f"validated_rankers -->> {validated_rankers}")
