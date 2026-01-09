from unittest import TestCase

from parameterized import parameterized
from config.ipcress_params import IpcressParameters


class TestIpcressParams(TestCase):

    def test_valid_parameters(self):
        params = {
            "write_ipcress_file": True,
            "min_size": 5,
            "max_size": 300
        }

        ipcress_params = IpcressParameters(params)

        self.assertTrue(ipcress_params.write_ipcress_file)
        self.assertEqual(ipcress_params.min_size, 5)
        self.assertEqual(ipcress_params.max_size, 300)

    def test_error_when_params_not_dict(self):
        with self.assertRaises(TypeError) as cm:
            IpcressParameters("not a dict")
        self.assertIn("ipcress_parameters must be a dict", str(cm.exception))

    @parameterized.expand([
        ("when_string", "STRING"),
        ("when_negative", -1),
        ("when_None", None),
    ])
    def test_write_ipcress_file_invalid(self, test_case, value):
        with self.assertRaises(ValueError) as cm:
            IpcressParameters({
                "write_ipcress_file": value,
                "min_size": 5,
                "max_size": 300
            })

        expected = f'Invalid config value for "write_ipcress_file": expected boolean, got {value} ({type(value).__name__})'
        self.assertIn(expected, str(cm.exception))

    @parameterized.expand([
        ("when_string", "STRING"),
        ("when_negative", -1),
        ("when_None", None),
    ])
    def test_min_size_invalid(self, test_case, min_size_value):
        with self.assertRaises(ValueError) as cm:
            IpcressParameters({
                "write_ipcress_file": True,
                "min_size": min_size_value,
                "max_size": 300
            })

        expected = f'Invalid config value for "min_size": expected non-negative int, got {min_size_value} ({type(min_size_value).__name__})'
        self.assertIn(expected, str(cm.exception))

    @parameterized.expand([
        ("when_string", "STRING"),
        ("when_negative", -1),
        ("when_None", None),
    ])
    def test_max_size_invalid(self, test_case, max_size_value):
        with self.assertRaises(ValueError) as cm:
            IpcressParameters({
                "write_ipcress_file": True,
                "min_size": 5,
                "max_size": max_size_value
            })

        expected = f'Invalid config value for "max_size": expected non-negative int, got {max_size_value} ({type(max_size_value).__name__})'
        self.assertIn(expected, str(cm.exception))
