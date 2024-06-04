from dataclasses import dataclass
from venv import logger
from dataclasses_json import dataclass_json
from typing import List

from custom_logger.custom_logger import CustomLogger

# Initialize logger
logger = CustomLogger(__name__)


@dataclass_json
@dataclass
class RankingCriteria():
    name: str
    ascending: bool
    column: str

    def __init__(self, **criteria):
        self.name = criteria["name"]
        for key, value in criteria.items():
            setattr(self, key, value)
       
class RankCriteriaValidator:
    def __init__(self, config: dict):
        self.csv_columns = config["csv_column_order"]
        self.ranking_columns = config["ranking_criteria"]
        self.created_criteria: List[RankingCriteria] = []
        self.is_validated: bool

        for criteria in self.ranking_columns:
            created = RankingCriteria.from_dict(criteria)
            self.created_criteria.append(created)

        validated, column_name, validated_rankers = validate_ranking_columns(self.created_criteria, self.csv_columns)

        if not validated:
            logger.warning(f"Invalid ranking column name {column_name}. Please Check configuration file")
            raise Exception(f"Invalid ranking column name {column_name}. Please Check configuration file")
        
        self.created_criteria = validated_rankers
    

def validate_ranking_columns(created_criteria, csv_columns):

    is_validated = True
    column_name = ""
    validated_rankers = []

    for ranking_criteria in created_criteria:
        if ranking_criteria.column in csv_columns:
            validated_rankers.append(ranking_criteria)
        else:
            is_validated = not is_validated
            column_name = ranking_criteria.column
            return is_validated, column_name, validated_rankers
        
    return is_validated, validated_rankers, validated_rankers