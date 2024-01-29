from unittest.mock import patch

from pyfakefs.fake_filesystem_unittest import TestCase

from post_primer_pairs import filter_primer_pairs, parse_primer_json, post_primer_data  # NOQA


class TestPostPrimerPairs(TestCase):
    def setUp(self):
        self.setUpPyfakefs()
        json_contents = '''[
            {
                "left": {
                    "seq": "CTGTTCTGACAGTAGAAAGGCA"
                },
                "pair": "exon1_2_LibAmp_0",
                "right": {
                    "seq": "AAGAATTTTCCCCAATGGTTGCT"
                },
                "score": "1.0",
                "targeton": "exon1"
            },
            {
                "left": {
                    "seq": "CTGTTCTGACAGTAGAAAGGCA"
                },
                "pair": "exon1_2_LibAmp_1",
                "right": {
                    "seq": "AAGAATTTTCCCCAATGGTTGC"
                },
                "score": "2.0",
                "targeton": "exon1"
            },
            {
                "left": {
                    "seq": "CTGTTCTGACAGTAGAAAGGCAT"
                },
                "pair": "exon1_2_LibAmp_2",
                "right": {
                    "seq": "AAGAATTTTCCCCAATGGTTGCT"
                },
                "score": "0.0",
                "targeton": "exon1"
            },
            {
                "left": {
                    "seq": "CTGTTCTGACAGTAGAAAGGCAT"
                },
                "pair": "exon1_2_LibAmp_3",
                "right": {
                    "seq": "AAGAATTTTCCCCAATGGTTGC"
                },
                "score": "0.0",
                "targeton": "exon1"
            },
            {
                "left": {
                    "seq": "CTGTTCTGACAGTAGAAAGGCA"
                },
                "pair": "exon2_2_LibAmp_0",
                "right": {
                    "seq": "AAGAATTTTCCCCAATGGTTGCT"
                },
                "score": "1.0",
                "targeton": "exon2"
            }
        ]'''
        self.fs.create_file('/test.json', contents=json_contents)
        self.primer_pair_dict = {
            'exon1': [
                {
                    "left": {"seq": "CTGTTCTGACAGTAGAAAGGCA"},
                    "pair": "exon1_2_LibAmp_0",
                    "right": {"seq": "AAGAATTTTCCCCAATGGTTGCT"},
                    "score": "1.0",
                    "targeton": "exon1",
                },
                {
                    "left": {"seq": "CTGTTCTGACAGTAGAAAGGCA"},
                    "pair": "exon1_2_LibAmp_1",
                    "right": {"seq": "AAGAATTTTCCCCAATGGTTGC"},
                    "score": "2.0",
                    "targeton": "exon1",
                },
                {
                    "left": {"seq": "CTGTTCTGACAGTAGAAAGGCAT"},
                    "pair": "exon1_2_LibAmp_2",
                    "right": {"seq": "AAGAATTTTCCCCAATGGTTGCT"},
                    "score": "0.0",
                    "targeton": "exon1",
                },
                {
                    "left": {"seq": "CTGTTCTGACAGTAGAAAGGCAT"},
                    "pair": "exon1_2_LibAmp_3",
                    "right": {"seq": "AAGAATTTTCCCCAATGGTTGC"},
                    "score": "0.0",
                    "targeton": "exon1",
                },
            ],
            'exon2': [
                {
                    "left": {"seq": "CTGTTCTGACAGTAGAAAGGCA"},
                    "pair": "exon2_2_LibAmp_0",
                    "right": {"seq": "AAGAATTTTCCCCAATGGTTGCT"},
                    "score": "1.0",
                    "targeton": "exon2",
                },
            ],
        }
        self.primer_pair_list = [
            {
                "left": {"seq": "CTGTTCTGACAGTAGAAAGGCAT"},
                "pair": "exon1_2_LibAmp_2",
                "right": {"seq": "AAGAATTTTCCCCAATGGTTGCT"},
                "score": "0.0",
                "targeton": "exon1",
            },
            {
                "left": {"seq": "CTGTTCTGACAGTAGAAAGGCAT"},
                "pair": "exon1_2_LibAmp_3",
                "right": {"seq": "AAGAATTTTCCCCAATGGTTGC"},
                "score": "0.0",
                "targeton": "exon1",
            },
            {
                "left": {"seq": "CTGTTCTGACAGTAGAAAGGCA"},
                "pair": "exon1_2_LibAmp_0",
                "right": {"seq": "AAGAATTTTCCCCAATGGTTGCT"},
                "score": "1.0",
                "targeton": "exon1",
            },
            {
                "left": {"seq": "CTGTTCTGACAGTAGAAAGGCA"},
                "pair": "exon2_2_LibAmp_0",
                "right": {"seq": "AAGAATTTTCCCCAATGGTTGCT"},
                "score": "1.0",
                "targeton": "exon2",
            },
        ]

    def test_parse_primer_json_outputs_expected_dict(self):
        # arrange
        expected = self.primer_pair_dict

        # act
        actual = parse_primer_json('/test.json')

        # assert
        self.assertEqual(expected, actual)

    @patch('builtins.print')
    def test_filter_primer_pairs_outputs_top_primer_pairs(self, mock_print):
        # arrange
        expected = self.primer_pair_list

        # act
        actual = filter_primer_pairs(self.primer_pair_dict)

        # assert
        self.assertEqual(expected, actual)
        mock_print.assert_called_once_with('Only 1 primer pair(s) for targeton: exon2')

    @patch('builtins.print')
    @patch('post_primer_pairs.requests.post')
    def test_post_primer_data_success(self, mock_post, mock_print):
        # arrange
        mock_post.return_value.status_code = 201

        # act
        post_primer_data(self.primer_pair_list)

        # assert
        mock_post.assert_called_once_with(
            'https://sge-service.link:8081/libamp',
            json=self.primer_pair_list,
            headers={'Content-Type': 'application/json'},
        )
        mock_print.assert_called_once_with('Successfully posted primers!')

    @patch('builtins.print')
    @patch('post_primer_pairs.requests.post')
    def test_post_primer_data_fail(self, mock_post, mock_print):
        # arrange
        mock_post.return_value.status_code = 401
        mock_post.return_value.reason = 'Error'

        # act
        post_primer_data(self.primer_pair_list)

        # assert
        mock_post.assert_called_once_with(
            'https://sge-service.link:8081/libamp',
            json=self.primer_pair_list,
            headers={'Content-Type': 'application/json'},
        )
        mock_print.assert_called_once_with('Issue with post request: 401 Error')
