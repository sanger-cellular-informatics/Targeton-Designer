import unittest
import pybedtools

from designer.slicer import get_slice_coordinates, get_slice_sequences

class TestSlicer(unittest.TestCase):

    def test_get_slice_coordinates(self):
        expected = [
            ('chr1', 50, 260),
            ('chr1', 55, 265),
            ('chr1', 60, 270),
            ('chr1', 65, 275),
            ('chr1', 70, 280),
            ('chr1', 75, 285),
            ('chr1', 80, 290),
            ('chr1', 85, 295),
            ('chr1', 90, 300)
        ]
        bed = pybedtools.BedTool('chr1\t100\t250', from_string=True)
        params = {
            'flank_5': 50,
            'flank_3': 50,
            'length': 210,
            'offset': 5
        }
        actual = get_slice_coordinates(bed, params)
        self.assertEqual(expected, actual)

    def test_get_slice_sequences(self):
        expected = {
            'chr1:5-10': 'AGTCT',
            'chr1:15-20': 'ATTTT'
        }
        bed = pybedtools.BedTool('chr1\t5\t10\nchr1\t15\t20', from_string=True)
        actual = get_slice_sequences(bed, pybedtools.example_filename('test.fa'))
        self.assertEqual(expected, actual)

if __name__ == '__main__':
    unittest.main()
