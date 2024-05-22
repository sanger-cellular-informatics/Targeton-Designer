import unittest
from unittest.mock import patch

from src.config.config import DesignerConfig


class TestDesignerConfig(unittest.TestCase):

    def setUp(self) -> None:

        self.mock_config = {
            "stringency_vector": [1, 2, 3],
            "csv_column_order": ["col1", "col2", "col3"],
            "filters":{
                "hap1": False
            }
        }
    
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
