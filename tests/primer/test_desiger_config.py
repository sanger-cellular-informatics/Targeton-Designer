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
                "hap1": True
            }
        }

        self.mock_config_with_no_filters_section = {
            "stringency_vector": [1, 2, 3],
            "csv_column_order": ["col1", "col2", "col3"]
        }

        self.mock_config_with_incorrect_filter_name = {
            "stringency_vector": [1, 2, 3],
            "csv_column_order": ["col1", "col2", "col3"],
            "filters":{
                "hap2": True
            }
        }

    
    def tearDown(self):
        # Remove the handler after each test to reset logging
        logger = logging.getLogger()
        logger.removeHandler(self.handler)
    
    @patch('src.config.config.DesignerConfig.read_config')
    def test_filters_exist_in_config_file(self, mock_config_with_no_filters_section):
        mock_config_with_no_filters_section.return_value = self.mock_config_with_no_filters_section.copy()


        _ = DesignerConfig(self.mock_config_with_no_filters_section)

        logs = self.handler.buffer.getvalue().strip()
       
        self.assertEqual(logs, "No filters present in configuration file.")
    
    @patch('src.config.config.DesignerConfig.read_config')
    def test_hap1_filter_is_enabled(self, mock_read_config):
        mock_read_config.return_value = self.mock_config.copy()
        config = DesignerConfig(self.mock_config)

        expected_params = {
            "stringency_vector": self.mock_config["stringency_vector"],
            "csv_column_order": self.mock_config["csv_column_order"],
            "filters": {
                "hap1": True
            }
        }
       
        self.assertEqual(config.params, expected_params)
    
    @patch('src.config.config.DesignerConfig.read_config')
    def test_hap1_filter_with_incorrect_name_or_filter_not_in_list(self, mock_config_with_incorrect_filter_name):
        mock_config_with_incorrect_filter_name.return_value = self.mock_config_with_incorrect_filter_name.copy()

        _ = DesignerConfig(self.mock_config_with_incorrect_filter_name)

        logs = self.handler.buffer.getvalue().strip()
       
        self.assertEqual(logs, "Please check if you have used the correct filter name. Eg.- Use hap1")
