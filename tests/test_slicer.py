import unittest
import pandas as pd
from io import StringIO

from designer.slicer import import_bed_file, add_slice_coordinates

class TestSlicer(unittest.TestCase):

    def test_import_bed_file(self):
        expected = pd.DataFrame({
            'chr': ['chr1'],
            'start': [123],
            'end': [456]
        })
        fake_bed = StringIO('chr1\t123\t456')
        actual = import_bed_file(fake_bed)
        pd.testing.assert_frame_equal(expected, actual)

    def test_add_slice_coordinates(self):
        expected_slices = {
            50: 260,
            55: 265,
            60: 270,
            65: 275,
            70: 280,
            75: 285,
            80: 290,
            85: 295,
            90: 300
        }
        expected = pd.DataFrame({
            'chr': ['chr1'],
            'start': [100],
            'end': [250],
            'slices': [expected_slices]
        })
        data = pd.DataFrame({
            'chr': ['chr1'],
            'start': [100],
            'end': [250]
        })
        actual = add_slice_coordinates(data, 50, 50, 210, 5)
        pd.testing.assert_frame_equal(expected, actual)

if __name__ == '__main__':
    unittest.main()
