import logging
import unittest
from tests.utils.utils import CapturingStreamHandler


from ranking.rank_criteria import StringencyCriteria, ProductSizeCriteria
from ranking.ranking_validator import RankCriteriaValidator


class TestRankCriteriaValidator(unittest.TestCase):

    def setUp(self):

        self.handler = CapturingStreamHandler()
        self.logger = self.handler.get_logger(self.handler)

        self.mocked_config = {
              "ranking_criteria": ["stringency", "product_size"],
        }
    
    def tearDown(self):
            # Remove the handler after each test to reset logging
            logger = logging.getLogger()
            logger.removeHandler(self.handler)

    def test_validate_ranking_criteria(self):

        result = RankCriteriaValidator(self.mocked_config).validated_criteria()

        self.assertEqual(
            result, [StringencyCriteria, ProductSizeCriteria]
        )

    def test_validate_ranking_criteria_when_invalid_ranking_criterion(self):
        mocked_invalid_ranking_priority_order = {
              "ranking_criteria": ["invalid_column"],
        }

        with self.assertRaises(ValueError) as error:
            RankCriteriaValidator(mocked_invalid_ranking_priority_order).validated_criteria()

        self.assertEqual(
            str(error.exception), 'Invalid ranking criteria: the given rank criterion "invalid_column" is not valid'
        )
    
    def test_check_if_no_ranking_criteria_provided(self):

        RankCriteriaValidator({}).validated_criteria()

        logs = self.handler.buffer.getvalue().strip()

        self.assertTrue("No Ranking criteria provided." in logs)

    def test_validate_ranking_criteria_when_repeated_rank_criterion(self):
        mocked_ranking_priority_order_repeated = {
              "ranking_criteria": ["stringency", "stringency"],
        }

        with self.assertRaises(ValueError) as error:
            RankCriteriaValidator(mocked_ranking_priority_order_repeated).validated_criteria()

        self.assertEqual(
            str(error.exception), 'Repeated ranking criteria: the given rank criterion "stringency" is repeated'
        )

