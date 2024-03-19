from unittest import TestCase
from unittest.mock import Mock

from primer.primer_pair import PrimerPair
from primer.filter.rename_primers import rename_primers


class TestRenamePrimers(TestCase):
    def setUp(self):
        pair1 = Mock(spec=PrimerPair)
        pair1.forward={'primer': 'primer1f'}
        pair1.reverse={'primer': 'primer1r'}

        pair2 = Mock(spec=PrimerPair)
        pair2.forward = {'primer': 'primer2f'}
        pair2.reverse = {'primer': 'primer2r'}

        self.primer_pairs = [pair1, pair2]

    def test_rename_primers(self):
        expected_primer_names1 = ('slice_name_LibAmpF_0', 'slice_name_LibAmpR_0')
        expected_primer_names2 = ('slice_name_LibAmpF_1', 'slice_name_LibAmpR_1')

        result = rename_primers(self.primer_pairs)

        self.assertEqual(
            (result[0].forward['primer'], result[0].reverse['primer']),
            expected_primer_names1
        )
        self.assertEqual(
            (result[1].forward['primer'], result[1].reverse['primer']),
            expected_primer_names2
        )
