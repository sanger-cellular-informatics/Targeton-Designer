from typing import List
import pandas as pd

from primer.write_primer_output import _get_primers_dataframe
from ranking.rank_criteria import RankingCriteria, ProductSizeCriteria, StringencyCriteria

from custom_logger.custom_logger import CustomLogger

# Initialize logger
logger = CustomLogger(__name__)

class Ranker:
    def __init__(self):
        self._ranking_order: List[RankingCriteria] = [StringencyCriteria, ProductSizeCriteria]

    def rank(self, primer_type:str, primer_pairs=[]) -> pd.DataFrame:

        # Primer pairs dataframe are grouped together.
        primers_df = _get_primers_dataframe(primer_pairs, primer_type)

        if primers_df.empty:
            logger.warning("No primer pairs to rank.")
            return primers_df

        if self._ranking_order:
            # Commented below line of code temporarily for future enhancement.
            # logger.info(f"Ranking is being applied by {', '.join([column.name for column in self._ranking_order])}")
            columns_to_sort = [column.column for column in self._ranking_order]
            is_ascending = [column.is_ascending for column in self._ranking_order]

            primers_df.sort_values(by=columns_to_sort, ascending=is_ascending, inplace=True)

        return primers_df
