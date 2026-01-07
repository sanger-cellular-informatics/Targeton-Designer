from unittest import TestCase

from parameterized import parameterized
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

    @parameterized.expand([
        ("flanking_string", "STRING"),
        ("flanking_min_size_negative", -1),
    ])
    @patch.object(DesignerConfig, 'read_config')
    @patch('config.config_helpers.logger.error')
    def test_flanking_region_wrong_value_error(self, test_case, flanking_value,mock_logger_error, mock_read_config):
        mock_read_config.return_value = {
            'flanking_region': flanking_value,
            'exclusion_region': 5,
            'stringency_vector': [],
            'csv_column_order': [],
            'filters': {},
            'ranking': {}
        }

        with self.assertRaises(SystemExit):
            DesignerConfig(args={})

        expected_error_message = (
            f'Invalid config value for "flanking_region": '
            f'expected non-negative int, got {flanking_value!r} ({type(flanking_value).__name__})'
        )

        mock_logger_error.assert_called_once_with(expected_error_message)

    @patch.object(DesignerConfig, 'read_config')
    @patch('config.config.logger.info')
    def test_flanking_region_no_value(self, mock_logger_info, mock_read_config):
        mock_read_config.return_value = {
            'flanking_region': 0,
            'exclusion_region': 5,
            'stringency_vector': [],
            'csv_column_order': [],
            'filters': {},
            'ranking': {}
        }

        DesignerConfig(args={})

        expected_info_message = (
            "flanking_region set to 0, so primer placement will not be restricted the flanking_region, and "
            "exclusion_region will be ignored."
        )
        mock_logger_info.assert_called_once_with(expected_info_message)

    @parameterized.expand([
        ("exclusion_region_string", "STRING"),
        ("exclusion_region_min_size_negative", -1),
    ])
    @patch.object(DesignerConfig, 'read_config')
    @patch('config.config_helpers.logger.error')
    def test_exclusion_region_wrong_value_error(self, test_case, exclusion_region_value, mock_logger_error, mock_read_config):
        mock_read_config.return_value = {
            'flanking_region': 5,
            'exclusion_region': exclusion_region_value,
            'stringency_vector': [],
            'csv_column_order': [],
            'filters': {},
            'ranking': {}
        }

        with self.assertRaises(SystemExit):
            DesignerConfig(args={})

        expected_error_message = (
            f'Invalid config value for "exclusion_region": '
            f'expected non-negative int, got {exclusion_region_value!r} ({type(exclusion_region_value).__name__})'
        )

        mock_logger_error.assert_called_once_with(expected_error_message)

    def test_read_config(self):
        expected = {'flanking_region': 150,
                    'exclusion_region': 5,
                    'stringency_vector': [1, 0.5, 0.1],
                    'csv_column_order': ["primer_type", "primer", "penalty", "stringency", "sequence", "primer_start", 
                                         "primer_end", "tm", "gc_percent", "self_any_th", "self_end_th", "hairpin_th", 
                                         "end_stability", "chromosome", "pre_targeton_start", "pre_targeton_end", 
                                         "product_size", "targeton_id", "pair_uid"]}

        result = DesignerConfig.read_config(
            default_config_path=self.default_config_path, config_path=self.config_path
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

        result = DesignerConfig.read_config(default_config_path=self.default_config_path)
        self.assertEqual(result, expected)

    def test_no_default_config_file_found(self):
        incorrect_default_path = 'tests/config/111.config.json'

        with self.assertRaises(FileNotFoundError):
            DesignerConfig.read_config(
                default_config_path=incorrect_default_path, config_path=self.config_path
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


class TestIpcressOutputDesignerConfig(TestCase):
    def setUp(self):
        self.config_path = 'tests/config_files/test_user_designer.config.json'
        self.starting_config = {
            'flanking_region': 150,
            'exclusion_region': 5,
            'stringency_vector': [],
            'csv_column_order': [],
            'filters': {},
            'ranking': {}
        }

    @patch.object(DesignerConfig, 'read_config')
    def test_write_ipcress_file_false_if_ipcress_parameters_missing(self, mock_read_config):
        mock_read_config.return_value = self.starting_config

        designer_config = DesignerConfig(args={'conf': self.config_path})

        self.assertFalse(designer_config.ipcress_params_write_file)

    
    @patch.object(DesignerConfig, 'read_config')    
    def test_get_write_ipcress_file_from_ipcress_parameters_false(self, mock_read_config):
        ipcress_parameters = {
            "write_ipcress_file": False,
            "min_size": 5,
            "max_size": 300
        }

        config = self.starting_config.copy()
        config["ipcress_parameters"] = ipcress_parameters
        mock_read_config.return_value = config

        designer_config = DesignerConfig(args={'conf': self.config_path})

        self.assertFalse(designer_config.ipcress_params_write_file)


    @patch.object(DesignerConfig, 'read_config')
    def test_get_ipcress_parameters_from_config(self, mock_read_config):
        ipcress_parameters = {
            "write_ipcress_file": True,
            "min_size": 5,
            "max_size": 300
        }

        config = self.starting_config.copy()
        config["ipcress_parameters"] = ipcress_parameters
        mock_read_config.return_value = config

        designer_config = DesignerConfig(args={'conf': self.config_path})

        self.assertTrue(designer_config.ipcress_params_write_file)
        self.assertEqual(designer_config.ipcress_params_min_size, 5)
        self.assertEqual(designer_config.ipcress_params_max_size, 300)

    @patch.object(DesignerConfig, 'read_config')
    @patch("config.config.logger.error")
    def test_get_ipcress_parameters_from_config_when_wrong_write_ipcress_file_value(self, mock_logger_error, mock_read_config):
        ipcress_parameters = {
            "write_ipcress_file": 15,
            "min_size": 5,
            "max_size": 300
        }

        config = self.starting_config.copy()
        config["ipcress_parameters"] = ipcress_parameters
        mock_read_config.return_value = config

        with self.assertRaises(SystemExit):
            DesignerConfig(args={"conf": self.config_path})

        expected_error_message = (
            f"ipcress_params.write_ipcress_file must be a boolean (true/false), got int"
        )

        mock_logger_error.assert_called_once_with(expected_error_message)


    @parameterized.expand([
        ("min_size_unformatted", "min_size", "STRING"),
        ("min_size_negative", "min_size", -1),
        ("max_size_string", "max_size", "STRING"),
        ("max_size_negative", "max_size", -10),
    ])
    @patch.object(DesignerConfig, "read_config")
    @patch("config.config_helpers.logger.error")
    def test_get_ipcress_parameters_from_config_when_wrong_values(
            self, test_case_name, param_name, param_value, mock_logger_error, mock_read_config
    ):
        # Default correct parameters
        ipcress_parameters = {"write_ipcress_file": True, "min_size": 100, "max_size": 300 }
        # Override the parameter to an invalid value
        ipcress_parameters[param_name] = param_value

        config = self.starting_config.copy()
        config["ipcress_parameters"] = ipcress_parameters
        mock_read_config.return_value = config

        with self.assertRaises(SystemExit):
            DesignerConfig(args={"conf": self.config_path})

        expected_error_message = (
            f'Invalid config value for "ipcress_params.{param_name}": '
            f'expected non-negative int, got {param_value!r} ({type(param_value).__name__})'
        )

        mock_logger_error.assert_called_once_with(expected_error_message)
        mock_logger_error.reset_mock()
