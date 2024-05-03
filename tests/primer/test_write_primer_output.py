import unittest
from io import StringIO
import sys
import pandas as pd
from pyfakefs.fake_filesystem_unittest import TestCase

from primer.write_primer_output import _reorder_columns

class TestWritePrimerOutputFiles(TestCase):
        def test_reorder_columns_when_duplicate_column_names(self):

            # Arrange list indicating column order (contains duplicates) and dataframe for reordering
            column_names = ['Name', 'Age', 'Name']
            data = {
                'Name': ['Juan', 'Maria', 'Pedro'],
                'Age': [25, 30, 35]
            }
            df = pd.DataFrame(data)

            # Act (reorder dataframe according to list)
            ordered_df = _reorder_columns(column_names, df)

            # Assertion
            pd.testing.assert_frame_equal(ordered_df, df)


        def test_reorder_columns_when_empty_column_names(self):
            expected_stdout = StringIO()
            sys.stdout = expected_stdout

            # Arrange list indicating column order (empty list) and dataframe for reordering
            column_names = []
            data = {
                'Name': ['Juan', 'Maria', 'Pedro'],
                'Age': [25, 30, 35]
            }
            df = pd.DataFrame(data)

            # Act (reorder dataframe according to list)
            ordered_df = _reorder_columns(column_names, df)

            # Assertion
            std_result = expected_stdout.getvalue().strip()

            pd.testing.assert_frame_equal(ordered_df, df)
            self.assertEqual(std_result, 
                                "Warning: empty csv_column_order list provided in config file, returning dataframe with default column order")


        def test_reorder_columns_when_all_column_names_wrong(self):
            expected_stdout = StringIO()
            sys.stdout = expected_stdout

            # Arrange list indicating column order (only wrong names inputed) and dataframe for reordering
            column_names = ['WRONG']
            data = {
                'Name': ['Juan', 'Maria', 'Pedro'],
                'Age': [25, 30, 35]
            }
            df = pd.DataFrame(data)

            # Act (reorder dataframe according to list)
            with self.assertRaises(ValueError) as value_error:
                _reorder_columns(column_names, df)

            # Assertion
            std_result = expected_stdout.getvalue().strip()

            self.assertEqual(str(value_error.exception), "All column names in config file are wrong")
            self.assertEqual(std_result, "Warning: 'WRONG' specified in config file not is not a column name")


        def test_reorder_columns_when_some_column_names_wrong(self):
            expected_stdout = StringIO()
            sys.stdout = expected_stdout

            # Arrange list indicating column order (some wrong names inputed) and dataframe for reordering
            column_names = ['WRONG', 'Name', 'Age']
            data = {
                'Name': ['Juan', 'Maria', 'Pedro'],
                'Age': [25, 30, 35]
            }
            df = pd.DataFrame(data)

            # Act (reorder dataframe according to list)
            ordered_df = _reorder_columns(column_names, df)

            # Assertion
            std_result = expected_stdout.getvalue().strip()

            pd.testing.assert_frame_equal(ordered_df, df[['Name', 'Age']])
            self.assertEqual(std_result, "Warning: 'WRONG' specified in config file not is not a column name")


        def test_reorder_columns_when_some_column_names_missing(self):
            expected_stdout = StringIO()
            sys.stdout = expected_stdout

            # Arrange list indicating column order (columns missing) and dataframe for reordering
            column_names = ['Name']
            data = {
                'Name': ['Juan', 'Maria', 'Pedro'],
                'Age': [25, 30, 35]
            }
            df = pd.DataFrame(data)

            # Act (reorder dataframe according to list)
            ordered_df = _reorder_columns(column_names, df)

            # Assertion
            std_result = expected_stdout.getvalue().strip()

            pd.testing.assert_frame_equal(ordered_df, df[['Name']])
            self.assertEqual(std_result, "'Age' column discarded as it is not in config file")


        def test_reorder_columns_success(self):

            # Arrange list indicating column order and dataframe for reordering
            column_names = ['Age', 'Name']
            data = {
                'Name': ['Juan', 'Maria', 'Pedro'],
                'Age': [25, 30, 35]
            }
            df = pd.DataFrame(data)

            # Act (reorder dataframe according to list)
            ordered_df = _reorder_columns(column_names, df)

            # Assertion
            pd.testing.assert_frame_equal(ordered_df, df[['Age', 'Name']])


if __name__ == '__main__':
    unittest.main()