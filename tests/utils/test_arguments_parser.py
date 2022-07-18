import unittest
import argparse

from unittest.mock import patch
from pyfakefs.fake_filesystem_unittest import TestCase

from src.utils.arguments_parser import positive_int, len_positive_int


class TestFileSystem(TestCase):
    def test_positive_int_positive_arg_success(self):
        # arrange
        test_arg = 1
        expected = 1

        # act
        actual = positive_int(test_arg)

        # assert
        self.assertEqual(actual, expected)

    def test_positive_int_negative_arg_fail(self):
        # arrange
        test_arg = -1
        expected = "Parameter must be above 0"

        # act
        with self.assertRaises(argparse.ArgumentTypeError) as exception_context:
            positive_int(test_arg)

        # assert
        self.assertEqual(str(exception_context.exception), expected)

    def test_positive_int_zero_fail(self):
        # arrange
        test_arg = 0
        expected = 'Parameter must be above 0'

        # act
        with self.assertRaises(argparse.ArgumentTypeError) as exception_context:
            positive_int(test_arg)

        # assert
        self.assertEqual(str(exception_context.exception), expected)

    def test_len_positive_int_positive_arg_success(self):
        # arrange
        test_arg = 1
        expected = 1

        # act
        actual = len_positive_int(test_arg)

        # assert
        self.assertEqual(actual, expected)

    def test_len_positive_int_zero_fail(self):
        # arrange
        test_arg = 0
        expected = 'Parameter must be above 0 and below 10000'

        # act
        with self.assertRaises(argparse.ArgumentTypeError) as exception_context:
            len_positive_int(test_arg)

        # assert
        self.assertEqual(str(exception_context.exception), expected)

    def test_len_positive_int_negative_arg_fail(self):
        # arrange
        test_arg = -1
        expected = "Parameter must be above 0 and below 10000"

        # act
        with self.assertRaises(argparse.ArgumentTypeError) as exception_context:
            len_positive_int(test_arg)

        # assert
        self.assertEqual(str(exception_context.exception), expected)

    def test_len_positive_int_large_arg_fail(self):
        # arrange
        test_arg = 10001
        expected = "Parameter must be above 0 and below 10000"

        # act
        with self.assertRaises(argparse.ArgumentTypeError) as exception_context:
            len_positive_int(test_arg)

        # assert
        self.assertEqual(str(exception_context.exception), expected)


if __name__ == '__main__':
    unittest.main()

