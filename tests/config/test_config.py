import unittest
from unittest.mock import patch
from pyfakefs.fake_filesystem_unittest import TestCase

from config.config import DesignerConfig, Primer3ParamsConfig
from utils.arguments_parser import ParsedInputArguments


class TestDesignerConfigClass(TestCase):
    def setUp(self):
        self.config_path = 'tests/config/designer.config.json'
        self.default_config_path = 'tests/config/designer_default.config.json'
        self.designer_config = DesignerConfig(self.config_path)

    def test_stringency_is_set(self):
        expected = [1, 0.5, 0.1]

        self.assertEqual(self.designer_config.stringency_vector, expected)

    def test_read_config(self):
        expected = {'stringency_vector': [1, 0.5, 0.1]}

        result = self.designer_config.read_config(
            default_config_file=self.default_config_path, config_file=self.config_path
        )

        self.assertEqual(result, expected)

    def test_no_config_file_found(self):
        incorrect_path = 'tests/config/111.config.json'

        with self.assertRaises(FileNotFoundError):
            config = DesignerConfig(incorrect_path)
            config.read_config(incorrect_path)

    def test_use_default_config(self):
        expected = {'stringency_vector': [1, 0.1]}

        config = DesignerConfig()
        result = config.read_config(default_config_file=self.default_config_path)

        self.assertEqual(result, expected)

    def test_no_default_config_file_found(self):
        incorrect_default_path = 'tests/config/111.config.json'

        config = DesignerConfig(self.config_path)

        with self.assertRaises(FileNotFoundError):
            config.read_config(
                default_config_file=incorrect_default_path, config_file=self.config_path
            )


class TestPrimer3ParamsClass(TestCase):
    def setUp(self):
        self.config_path = 'tests/config/primer3_params_simple.config.json'
        self.default_config_path = 'tests/config/primer3_params_default.config.json'
        self.p3_params_config = Primer3ParamsConfig(self.config_path)

    def test_params_are_set(self):
        expected = {"PRIMER_TASK": "pick_cloning_primers"}

        self.assertEqual(self.p3_params_config.params, expected)

    def test_read_config(self):
        expected = {"PRIMER_TASK": "pick_cloning_primers"}

        result = self.p3_params_config.read_config(
            default_config_file=self.default_config_path, config_file=self.config_path
        )

        self.assertEqual(result, expected)

    def test_no_config_file_found(self):
        incorrect_path = 'tests/config/111.config.json'

        with self.assertRaises(FileNotFoundError):
            config = Primer3ParamsConfig(incorrect_path)
            config.read_config(incorrect_path)

    def test_use_default_config(self):
        expected = {"PRIMER_TASK": "generic"}

        config = Primer3ParamsConfig()
        result = config.read_config(default_config_file=self.default_config_path)

        self.assertEqual(result, expected)
