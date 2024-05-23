import logging
import unittest
from unittest.mock import patch

from tests.utils.utils import CapturingStreamHandler
from src.config.config import DesignerConfig


class TestDesignerConfig(unittest.TestCase):

    def setUp(self) -> None:

        # Create a custom stream handler to capture logs
        self.handler = CapturingStreamHandler()
        # Get the logger and set the level to capture warnings (adjust if needed)
        logger = logging.getLogger()
        logger.setLevel(logging.DEBUG)
        logger.addHandler(self.handler)

        self.mock_config = {
            "stringency_vector": [1, 2, 3],
            "csv_column_order": ["col1", "col2", "col3"],
            "filters":{
                "hap1": False
            }
        }
    
    def tearDown(self):
        # Remove the handler after each test to reset logging
        logger = logging.getLogger()
        logger.removeHandler(self.handler)
    
    @patch('src.config.config.DesignerConfig.read_config')
    def test_hap1_filter_is_enabled(self, mock_read_config):
        mock_read_config.return_value = self.mock_config.copy()
        config = DesignerConfig(self.mock_config, "hap1")

        expected_params = {
            "stringency_vector": self.mock_config["stringency_vector"],
            "csv_column_order": self.mock_config["csv_column_order"],
            "filters": {
                "hap1": True
            }
        }
       
        self.assertEqual(config.params, expected_params)
    
    @patch('src.config.config.DesignerConfig.read_config')
    def test_hap1_filter_with_incorrect_name(self, mock_read_config):
        mock_read_config.return_value = self.mock_config.copy()
        incorrect_filter_name = "hap2"
        _ = DesignerConfig(self.mock_config, incorrect_filter_name)

        logs = self.handler.buffer.getvalue().strip()
       
        self.assertEqual(logs, f"Please check if you have used correct filter name. Eg.- use hap1 not {incorrect_filter_name}")
    
    @patch('src.config.config.DesignerConfig.read_config')
    def test_hap1_filter_is_not_used(self, mock_read_config):
        mock_read_config.return_value = self.mock_config.copy()
        config = DesignerConfig(self.mock_config)

        expected_params = {
            "stringency_vector": self.mock_config["stringency_vector"],
            "csv_column_order": self.mock_config["csv_column_order"],
            "filters": {
                "hap1": False
            }
        }

        self.assertEqual(config.params, expected_params)
