from typing import List
import pandas as pd
import sys

from primer.write_primer_output import _get_primers_dataframe
from ranking.rank_criteria import RankingCriteria, ProductSizeCriteria, StringencyCriteria

from custom_logger.custom_logger import CustomLogger

# Initialize logger
logger = CustomLogger(__name__)


class Ranker:
    def __init__(self, ranking_config: dict):

        self._ranking_criteria: List[RankingCriteria] = [StringencyCriteria, ProductSizeCriteria]
        _ranking_criteria_names: List[str] = [criterion.name for criterion
                                              in self._ranking_criteria]

        ranking_retained: List[str] = [key for key, value in ranking_config.items()
                                       if value is True]

        if not ranking_config:
            msg: str = "No ranking criteria defined in config file under the ranking key"
            logger.warning(msg)
        else:
            for criterion in _ranking_criteria_names:
                if criterion not in ranking_config.keys():
                    msg: str = f"'{criterion}' missing in config file for ranking " \
                          "- Will not be used for ranking"
                    logger.warning(msg)

        incorrect_values = [key for key, value in ranking_config.items()
                            if not isinstance(value, bool)]
        if incorrect_values:
            msg: str = f"Wrong value(s) provided for '{', '.join(incorrect_values)}' in config " \
                  "file (only takes true or false). Unable to apply ranking - Exiting programme"
            logger.error(msg)
            sys.exit(1)

        incorrect_keys: List[str] = [key for key in ranking_config.keys()
                                     if key not in _ranking_criteria_names]
        if incorrect_keys:
            msg: str = "Invalid name(s) provided for ranking in config file: " \
                  f"'{', '.join(incorrect_keys)}'. The only valid names are: " \
                  f"{', '.join(_ranking_criteria_names)}. Unable to apply ranking - " \
                  "Exiting programme"
            logger.error(msg)
            sys.exit(1)

        self._ranking_order: List[RankingCriteria] = []
        for criterion in ranking_retained:
            criterion_index: int = _ranking_criteria_names.index(criterion)
            self._ranking_order.append(self._ranking_criteria[criterion_index])

    def rank(self, primer_type: str, primer_pairs=list) -> pd.DataFrame:

        # Primer pairs dataframe are grouped together.
        primers_df = _get_primers_dataframe(primer_pairs, primer_type)

        if primers_df.empty:
            logger.warning("No primer pairs to rank.")
            return primers_df

        if self._ranking_order:
            logger.info(f"Ranking is being applied by {', '.join([column.name for column in self._ranking_order])}")
            columns_to_sort = [column.column for column in self._ranking_order]
            is_ascending = [column.is_ascending for column in self._ranking_order]

            # kind = "stable" added to ensure stable sorting when sorting on one column
            primers_df.sort_values(by=columns_to_sort, ascending=is_ascending, inplace=True,
                                   kind="stable")

        else:
            logger.info("No ranking applied")

        return primers_df
