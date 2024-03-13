import unittest
from unittest.mock import patch
from pyfakefs.fake_filesystem_unittest import TestCase

from config.config import Config, prepare_config
from utils.arguments_parser import ParsedInputArguments


class TestConfigClass(TestCase):
    def test_stringency_is_set(self):
        config_path = 'tests/config/designer.config.json'
        expected = [1, 0.5, 0.1]

        config = Config(config_path)

        self.assertEqual(config.stringency_vector, expected)


class TestPrepareConfig(TestCase):
    def test_prepare_config(self):
        config_path = 'tests/config/designer.config.json'
        expected = {'stringency_vector': [1, 0.5, 0.1]}

        result = prepare_config(config_path)

        self.assertEqual(result, expected)

    def test_no_config_file_found(self):
        incorrect_path = 'tests/config/111.config.json'

        with self.assertRaises(FileNotFoundError):
            prepare_config(incorrect_path)

    def test_use_default_config(self):
        default_config_path = 'tests/config/designer_default.config.json'
        expected = {'stringency_vector': [1, 0.1]}

        result = prepare_config(None, default_config_path)

        self.assertEqual(result, expected)

    def test_no_default_config_file_found(self):
        incorrect_default_path = 'tests/config/111.config.json'
        config_path = 'tests/config/designer.config.json'

        with self.assertRaises(FileNotFoundError):
            prepare_config(config_path, incorrect_default_path)