import unittest
from pyfakefs.fake_filesystem_unittest import TestCase

from config.config import prepare_config


class TestPrepareConfig(TestCase):
    def test_prepare_config(self):
        config_path = 'tests/config/designer.config.json'
        expected = {'stringency_vector': [1, 0.5, 0.1]}

        result = prepare_config(config_path)

        self.assertEqual(result, expected)

    def test_use_default_config(self):
        default_config_path = 'tests/config/designer_default.config.json'
        expected = {'stringency_vector': [1, 0.1]}

        result = prepare_config(None, default_config_path)

        self.assertEqual(result, expected)
