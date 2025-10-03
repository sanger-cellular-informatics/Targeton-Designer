import logging
import unittest
from unittest.mock import patch

from tests.utils.utils import CapturingStreamHandler
from src.config.config import DesignerConfig


class TestDesignerConfig(unittest.TestCase):

    def setUp(self) -> None:

        # Create a custom stream handler to capture logs
        self.handler = CapturingStreamHandler()
        self.logger = self.handler.get_logger(self.handler)

        self.mock_config_with_no_filters_section = {
                    "stringency_vector": [1, 2, 3],
                    "csv_column_order": ["col1", "col2", "col3"],
                    "ranking": {},
                    "flanking_region": 0,
                    "exclusion_region": 0
                }

    def tearDown(self):
        # Remove the handler after each test to reset logging
        logger = logging.getLogger()
        logger.removeHandler(self.handler)

    @patch('src.config.config.DesignerConfig.read_config')
    def test_filters_exist_in_config_file(self, mock_config_with_no_filters_section):

        mock_config_with_no_filters_section.return_value = self.mock_config_with_no_filters_section.copy()
        expected_config = {}
        expected_config["filters"] = {"duplicates": True}
        config = DesignerConfig(self.mock_config_with_no_filters_section)

        self.assertEqual(config.filters, expected_config["filters"])

