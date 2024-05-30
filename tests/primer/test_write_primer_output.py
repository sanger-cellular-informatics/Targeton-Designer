import logging
import unittest
from io import StringIO
import pandas as pd
from pyfakefs.fake_filesystem_unittest import TestCase

from tests.utils.utils import CapturingStreamHandler

from collections import defaultdict
from primer.primer_pair import PrimerPair
from primer.designed_primer import DesignedPrimer, Interval
from primer.write_primer_output import _reorder_columns, _add_primer_pair


class TestWritePrimerOutputFiles(TestCase):
        
        def setUp(self):
            # Create a custom stream handler to capture logs
            self.handler = CapturingStreamHandler()
            self.logger = self.handler.get_logger(self.handler)

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
            logs = self.handler.buffer.getvalue().strip()
            self.assertEqual(logs, "'Name' duplicated in config file, only first instance retained")
        
    
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
            self.assertEqual(logs, "Empty csv_column_order list provided in config file, returning dataframe with default column order")


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
            self.assertEqual(logs, "'WRONG' specified in config file not is not a column name")


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
            self.assertEqual(logs, "'WRONG' specified in config file not is not a column name")


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


class TestDataFrameBuild(TestCase):
    def test__add_primer_pair_when_dict_empty(self):
        test_dict = defaultdict(list)
        primer_type = "LibAmp"

        pair = PrimerPair(
                pair_id="ABCD_0_str1",
                chromosome="1",
                pre_targeton_start=42930996,
                pre_targeton_end=42931206,
                product_size=129,
                stringency=1,
                targeton_id="ABCD",
                uid="uid0")

        pair.forward = DesignedPrimer(
                name="ABCD_LibAmpF_0",
                penalty=0.16,
                pair_id="ABCD_0_str1",
                sequence="AGAAAACGCTGGTGTGGTGA",
                coords=Interval(start=169, end=20),
                primer_start=42931146,
                primer_end=42931166,
                strand="+",
                tm=60.0,
                gc_percent=50.0,
                self_any_th=0.0,
                self_end_th=0.0,
                hairpin_th=0.0,
                end_stability=4.02
            )

        pair.reverse = DesignedPrimer(
                name="ABCD_LibAmpR_0",
                penalty=0.30,
                pair_id="ABCD_0_str1",
                sequence="TGTTGCTCTTTTCCCAGGCT",
                coords=Interval(start=41, end=20),
                primer_start=42930996,
                primer_end=42931016,
                strand="-",
                tm=59.8,
                gc_percent=50.1,
                self_any_th=0.1,
                self_end_th=0.2,
                hairpin_th=0.2,
                end_stability=4.58
            )

        expected_dict = defaultdict(list)
        expected_dict = {
                'primer_type': ['LibAmp', 'LibAmp'],
                'primer': ['ABCD_LibAmpF_0', 'ABCD_LibAmpR_0'],
                'penalty': [0.16, 0.3],
                'sequence': ['AGAAAACGCTGGTGTGGTGA', 'TGTTGCTCTTTTCCCAGGCT'],
                'primer_start': [42931146, 42930996],
                'primer_end': [42931166, 42931016],
                'tm': [60.0, 59.8],
                'gc_percent': [50.0, 50.1],
                'self_any_th': [0.0, 0.1],
                'self_end_th': [0.0, 0.2],
                'hairpin_th': [0.0, 0.2],
                'end_stability': [4.02, 4.58],
                'pair_uid': ['uid0', 'uid0'],
                'stringency': [1, 1],
                'chromosome': ['1', '1'],
                'pre_targeton_start': [42930996, 42930996],
                'pre_targeton_end': [42931206, 42931206],
                'product_size': [129, 129],
                'targeton_id': ['ABCD', 'ABCD']
                }

        self.maxDiff = None
        self.assertIsNone(_add_primer_pair(test_dict, pair, primer_type))
        self.assertDictEqual(test_dict, expected_dict)

    def test__add_primer_pair_when_dict_not_empty(self):
        test_dict = defaultdict(list)
        primer_type = "LibAmp"

        test_dict = {
                'primer_type': ['LibAmp', 'LibAmp'],
                'primer': ['ABCD_LibAmpF_0', 'ABCD_LibAmpR_0'],
                'penalty': [0.16, 0.3],
                'sequence': ['AGAAAACGCTGGTGTGGTGA', 'TGTTGCTCTTTTCCCAGGCT'],
                'primer_start': [42931146, 42930996],
                'primer_end': [42931166, 42931016],
                'tm': [60.0, 59.8],
                'gc_percent': [50.0, 50.1],
                'self_any_th': [0.0, 0.1],
                'self_end_th': [0.0, 0.2],
                'hairpin_th': [0.0, 0.2],
                'end_stability': [4.02, 4.58],
                'pair_uid': ['uid0', 'uid0'],
                'stringency': [1, 1],
                'chromosome': ['1', '1'],
                'pre_targeton_start': [42930996, 42930996],
                'pre_targeton_end': [42931206, 42931206],
                'product_size': [129, 129],
                'targeton_id': ['ABCD', 'ABCD']
            }

        pair = PrimerPair(
                pair_id="ABCD_1_str1",
                chromosome="1",
                pre_targeton_start=42930996,
                pre_targeton_end=42931206,
                product_size=122,
                stringency=1,
                targeton_id="ABCD",
                uid="uid1")

        pair.forward = DesignedPrimer(
                name="ABCD_LibAmpF_1",
                penalty=0.17,
                pair_id="ABCD_1_str1",
                sequence="GCTGGTGTGGTGATATGCCT",
                coords=Interval(start=169, end=20),
                primer_start=42931139,
                primer_end=42931159,
                strand="+",
                tm=60.1,
                gc_percent=50.2,
                self_any_th=0.2,
                self_end_th=0.3,
                hairpin_th=0.3,
                end_stability=4.75
            )

        pair.reverse = DesignedPrimer(
                name="ABCD_LibAmpR_1",
                penalty=0.31,
                pair_id="ABCD_1_str1",
                sequence="TGTTGCTCTTTTCCCAGGCT",
                coords=Interval(start=41, end=20),
                primer_start=42930996,
                primer_end=42931016,
                strand="-",
                tm=59.9,
                gc_percent=50.3,
                self_any_th=0.3,
                self_end_th=0.4,
                hairpin_th=0.4,
                end_stability=4.58
            )

        expected_dict = defaultdict(list)
        expected_dict = {
                'primer_type': ['LibAmp', 'LibAmp', 'LibAmp', 'LibAmp'],
                'primer': ['ABCD_LibAmpF_0', 'ABCD_LibAmpR_0',
                           'ABCD_LibAmpF_1', 'ABCD_LibAmpR_1'],
                'penalty': [0.16, 0.3, 0.17, 0.31],
                'sequence': ['AGAAAACGCTGGTGTGGTGA', 'TGTTGCTCTTTTCCCAGGCT',
                             'GCTGGTGTGGTGATATGCCT', 'TGTTGCTCTTTTCCCAGGCT'],
                'primer_start': [42931146, 42930996, 42931139, 42930996],
                'primer_end': [42931166, 42931016, 42931159, 42931016],
                'tm': [60.0, 59.8, 60.1, 59.9],
                'gc_percent': [50.0, 50.1, 50.2, 50.3],
                'self_any_th': [0.0, 0.1, 0.2, 0.3],
                'self_end_th': [0.0, 0.2, 0.3, 0.4],
                'hairpin_th': [0.0, 0.2, 0.3, 0.4],
                'end_stability': [4.02, 4.58, 4.75, 4.58],
                'pair_uid': ['uid0', 'uid0', 'uid1', 'uid1'],
                'stringency': [1, 1, 1, 1],
                'chromosome': ['1', '1', '1', '1'],
                'pre_targeton_start': [42930996, 42930996, 42930996, 42930996],
                'pre_targeton_end': [42931206, 42931206, 42931206, 42931206],
                'product_size': [129, 129, 122, 122],
                'targeton_id': ['ABCD', 'ABCD', 'ABCD', 'ABCD']
            }

        self.maxDiff = None
        self.assertIsNone(_add_primer_pair(test_dict, pair, primer_type))
        self.assertDictEqual(test_dict, expected_dict)


if __name__ == '__main__':
    unittest.main()
