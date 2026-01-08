from unittest import TestCase

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

    def test_params_not_dict(self):
        with self.assertRaises(TypeError) as cm:
            IpcressParameters("not a dict")
        self.assertIn("ipcress_parameters must be a dict", str(cm.exception))

    def test_write_ipcress_file_invalid(self):
        invalid_values = [None, 0, 1, "true", []]
        for value in invalid_values:
            with self.assertRaises(ValueError) as cm:
                IpcressParameters({
                    "write_ipcress_file": value,
                    "min_size": 5,
                    "max_size": 300
                })
            self.assertIn("ipcress_parameters.write_ipcress_file must be a boolean", str(cm.exception))

    def test_min_size_invalid(self):
        invalid_values = [-1, -100, 5.5, "10", None]
        for value in invalid_values:
            with self.assertRaises(ValueError) as cm:
                IpcressParameters({
                    "write_ipcress_file": True,
                    "min_size": value,
                    "max_size": 300
                })
            self.assertIn("ipcress_parameters.min_size must be a non-negative integer", str(cm.exception))

    def test_max_size_invalid(self):
        invalid_values = [-1, -100, 5.5, "10", None]
        for value in invalid_values:
            with self.assertRaises(ValueError) as cm:
                IpcressParameters({
                    "write_ipcress_file": True,
                    "min_size": 5,
                    "max_size": value
                })
            self.assertIn("ipcress_parameters.max_size must be a non-negative integer", str(cm.exception))
