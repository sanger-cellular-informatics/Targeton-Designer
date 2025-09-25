from unittest import TestCase

from unittest.mock import patch, call

from config.config import DesignerConfig
from utils.file_system import parse_json


class TestDesignerConfigClass(TestCase):
    def setUp(self):
        self.config_path = 'tests/config_files/test_user_designer.config.json'
        self.default_config_path = 'tests/config_files/test_default_designer.config.json'
        self.config_with_params_path = 'tests/config_files/test_user_designer.config_with_params.json'

    def test_stringency_is_set(self):
        expected = [1, 0.5, 0.1]

        designer_config = DesignerConfig(args={'conf': self.config_path})

        self.assertEqual(designer_config.stringency_vector, expected)

    def test_column_order_is_provided(self):
        expected = ["primer_type", "primer", "penalty", "stringency", "sequence", "primer_start", 
                    "primer_end", "tm", "gc_percent", "self_any_th", "self_end_th", "hairpin_th", 
                    "end_stability", "chromosome", "pre_targeton_start", "pre_targeton_end", 
                    "product_size", "targeton_id", "pair_uid"]

        designer_config = DesignerConfig(args={'conf': self.config_path})

        self.assertEqual(designer_config.csv_column_order, expected)

    def test_region_padding_errors(self):
        expected = ["primer_type", "primer", "penalty", "stringency", "sequence", "primer_start", 
                    "primer_end", "tm", "gc_percent", "self_any_th", "self_end_th", "hairpin_th", 
                    "end_stability", "chromosome", "pre_targeton_start", "pre_targeton_end", 
                    "product_size", "targeton_id", "pair_uid"]

        designer_config = DesignerConfig(args={'conf': self.config_path})

        self.assertEqual(designer_config.csv_column_order, expected)

    def test_read_config(self):
        expected = {'padding_region': 150,
                    'avoid_region': 5,
                    'stringency_vector': [1, 0.5, 0.1],
                    'csv_column_order': ["primer_type", "primer", "penalty", "stringency", "sequence", "primer_start", 
                                         "primer_end", "tm", "gc_percent", "self_any_th", "self_end_th", "hairpin_th", 
                                         "end_stability", "chromosome", "pre_targeton_start", "pre_targeton_end", 
                                         "product_size", "targeton_id", "pair_uid"]}

        result = DesignerConfig.read_config(
            default_config_file=self.default_config_path, config_file=self.config_path
        )
        self.assertEqual(result, expected)

    def test_no_config_file_found(self):
        incorrect_path = 'tests/config/111.config.json'

        with self.assertRaises(FileNotFoundError):
            DesignerConfig.read_config(incorrect_path)

    def test_use_default_config(self):
        expected = {'stringency_vector': [1, 0.1],
                    'csv_column_order': ["primer_type", "primer", "penalty", "stringency", "sequence", "primer_start", 
                                         "primer_end", "tm", "gc_percent", "self_any_th", "self_end_th", "hairpin_th", 
                                         "end_stability", "chromosome", "pre_targeton_start", "pre_targeton_end", 
                                         "product_size", "targeton_id", "pair_uid"]}

        result = DesignerConfig.read_config(default_config_file=self.default_config_path)
        self.assertEqual(result, expected)

    def test_no_default_config_file_found(self):
        incorrect_default_path = 'tests/config/111.config.json'

        with self.assertRaises(FileNotFoundError):
            DesignerConfig.read_config(
                default_config_file=incorrect_default_path, config_file=self.config_path
            )

    @patch('config.config.parse_json')
    def test_args_params_override_config_file_params(self, mock_parse_json):
        # Arrange
        primer3_params = {"key": "second_value"}

        def side_effect(*args, **kwargs):
            if mock_parse_json.call_count == 3:
                return primer3_params
            else:
                return parse_json(*args, **kwargs)

        mock_parse_json.side_effect = side_effect
        args = {
            "dir": "OUTPUT_DIR_FROM_ARGS",
            "fasta": "FASTA_DIR_FROM_ARGS",
            "primer3_params": "PRIMER3_PARAMS_FROM_ARGS",
            "conf": self.config_with_params_path
        }

        # Act
        config = DesignerConfig(args)

        # Assert
        self.assertEqual(config.prefix_output_dir, "OUTPUT_DIR_FROM_ARGS")
        self.assertEqual(config.fasta, "FASTA_DIR_FROM_ARGS")
        self.assertEqual(config.primer3_params, primer3_params)
        self.assertEqual(mock_parse_json.call_args_list[2], call("PRIMER3_PARAMS_FROM_ARGS"))

    @patch('config.config.parse_json')
    def test_get_params_from_config_file_when_no_args(self, mock_parse_json):
        # Arrange
        primer3_params = {"key": "second_value"}

        def side_effect(*args, **kwargs):
            if mock_parse_json.call_count == 3:
                return primer3_params
            else:
                return parse_json(*args, **kwargs)

        mock_parse_json.side_effect = side_effect
        args = {
            "conf": self.config_with_params_path
        }

        # Act
        config = DesignerConfig(args)

        # Assert
        json_config_expected = parse_json(self.config_with_params_path)
        self.assertEqual(config.prefix_output_dir, json_config_expected["dir"])
        self.assertEqual(config.fasta, json_config_expected["fasta"])
        self.assertEqual(config.primer3_params, primer3_params)
        self.assertEqual(mock_parse_json.call_args_list[2], call(json_config_expected["primer3_params"]))
