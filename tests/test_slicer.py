import unittest
import pandas as pd
from io import StringIO

from designer.slicer import import_bed_file

class TestSlicer(unittest.TestCase):
    @classmethod
    def test_import_bed_file(self):
        expected = pd.DataFrame({'chr': ['chr1'], 'start': [123], 'end': [456]})
        fake_bed = StringIO('chr1\t123\t456')
        actual = import_bed_file(fake_bed)
        pd.testing.assert_frame_equal(expected, actual)

if __name__ == '__main__':
    unittest.main()
