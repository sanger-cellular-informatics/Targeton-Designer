import logging
import subprocess
import unittest
from io import StringIO
import pandas as pd
from pyfakefs.fake_filesystem_unittest import TestCase


from primer.write_primer_output import _reorder_columns

class CapturingStreamHandler(logging.StreamHandler):
  """Custom StreamHandler to capture logs in memory."""
  def __init__(self):
    super().__init__()
    self.buffer = StringIO()
    self.stream = self.buffer

class TestWritePrimerOutputFiles(TestCase):
        
        def setUp(self):
            # Create a custom stream handler to capture logs
            self.handler = CapturingStreamHandler()
            # Get the logger and set the level to capture warnings (adjust if needed)
            logger = logging.getLogger()
            logger.setLevel(logging.DEBUG)
            logger.addHandler(self.handler)
            subprocess.call(['rm', '-rf', 'logs/*'])

        def tearDown(self):
            # Remove the handler after each test to reset logging
            logger = logging.getLogger()
            logger.removeHandler(self.handler)
    
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

            # Arrange list indicating column order (empty list) and dataframe for reordering
            column_names = []
            data = {
                'Name': ['Juan', 'Maria', 'Pedro'],
                'Age': [25, 30, 35]
            }
            df = pd.DataFrame(data)

            ordered_df = _reorder_columns(column_names, df)

            pd.testing.assert_frame_equal(ordered_df, df)
            logs = self.handler.buffer.getvalue().strip()
            self.assertEqual(logs, "Warning: empty csv_column_order list provided in config file, returning dataframe with default column order")


        def test_reorder_columns_when_all_column_names_wrong(self):

            # Arrange list indicating column order (only wrong names inputed) and dataframe for reordering
            column_names = ['WRONG']
            data = {
                'Name': ['Juan', 'Maria', 'Pedro'],
                'Age': [25, 30, 35]
            }
            df = pd.DataFrame(data)

            with self.assertRaises(ValueError) as value_error:
                _reorder_columns(column_names, df)
                self.assertEqual(str(value_error.exception), "All column names in config file are wrong")
            
            logs = self.handler.buffer.getvalue().strip()
            self.assertEqual(logs, "Warning: 'WRONG' specified in config file not is not a column name")
                



        def test_reorder_columns_when_some_column_names_wrong(self):
      
            # Arrange list indicating column order (some wrong names inputed) and dataframe for reordering
            column_names = ['WRONG', 'Name', 'Age']
            data = {
                'Name': ['Juan', 'Maria', 'Pedro'],
                'Age': [25, 30, 35]
            }
            df = pd.DataFrame(data)

            # Act (reorder dataframe according to list)
            ordered_df = _reorder_columns(column_names, df)

            pd.testing.assert_frame_equal(ordered_df, df[['Name', 'Age']])
            logs = self.handler.buffer.getvalue().strip()
            self.assertEqual(logs, "Warning: 'WRONG' specified in config file not is not a column name")


        def test_reorder_columns_when_some_column_names_missing(self):
    
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
            pd.testing.assert_frame_equal(ordered_df, df[['Name']])

            logs = self.handler.buffer.getvalue().strip()
            self.assertEqual(logs, "'Age' column discarded as it is not in config file")


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