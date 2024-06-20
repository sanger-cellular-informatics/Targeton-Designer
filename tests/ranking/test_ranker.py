import logging
from io import StringIO
from unittest import TestCase
from unittest.mock import Mock
import pandas as pd

from pandas.testing import assert_frame_equal
from pysam import index
from ranking.ranker import Ranker
from primer.primer_pair import PrimerPair
from primer.designed_primer import DesignedPrimer, Interval
from tests.utils.utils import CapturingStreamHandler

class RankerTest(TestCase):

    def setUp(self):

        # Create a custom stream handler to capture logs
        self.handler = CapturingStreamHandler()
        self.logger = self.handler.get_logger(self.handler)

        # Mock config
        # mocked_ranking_config_all_true also tests for order
        # mocked_ranking_wrong_key also tests for missing key
        self.mocked_ranking_config_all_true = {'product_size': True, 'stringency': True}
        self.mocked_ranking_config_empty = {}
        self.mocked_ranking_wrong_key = {'product_sizes': True, 'stringency': True}
        self.mocked_ranking_wrong_value = {'product_size': 'True', 'stringency': True}
        self.mocked_ranking_config_all_false = {'product_size': False, 'stringency': False}
        self.mocked_ranking_config_one_true = {'product_size': True, 'stringency': False}

        # Mock primers and pairs
        self.mocked_primer_forward = DesignedPrimer(
            name="primer_fr",
            penalty=0.5,
            pair_id="pair_id",
            sequence="ATCGATCG",
            coords=Interval(start=199, end=18),
            primer_start=11540,
            primer_end=11545,
            strand="+",
            tm=60.0,
            gc_percent=50.0,
            self_any_th=30.0,
            self_end_th=10.0,
            hairpin_th=20.0,
            end_stability=25.0
        )

        self.mocked_primer_reverse = DesignedPrimer(
            name="primer_rv",
            penalty=0.5,
            pair_id="pair_id",
            sequence="CTCGATCG",
            coords=Interval(start=299, end=18),
            primer_start=11640,
            primer_end=11645,
            strand="-",
            tm=60.0,
            gc_percent=50.0,
            self_any_th=30.0,
            self_end_th=10.0,
            hairpin_th=20.0,
            end_stability=25.0
        )

        self.mocked_pair_1 = PrimerPair(
            pair_id="pair1",
            chromosome="ch1",
            pre_targeton_start=11540,
            pre_targeton_end=11645,
            product_size=3,
            stringency=0.1,
            targeton_id="targeton_id",
            uid="uid1")

        self.mocked_pair_1.forward = self.mocked_primer_forward
        self.mocked_pair_1.reverse = self.mocked_primer_reverse

        self.mocked_pair_2 = PrimerPair(
            pair_id="pair2",
            chromosome="ch1",
            pre_targeton_start=11540,
            pre_targeton_end=11645,
            product_size=1,
            stringency=0.2,
            targeton_id="targeton_id",
            uid="uid2")

        self.mocked_pair_2.forward = self.mocked_primer_forward
        self.mocked_pair_2.reverse = self.mocked_primer_reverse

        self.mocked_pair_3 = PrimerPair(
            pair_id="pair3",
            chromosome="ch1",
            pre_targeton_start=11540,
            pre_targeton_end=11645,
            product_size=3,
            stringency=0.3,
            targeton_id="targeton_id",
            uid="uid3")

        self.mocked_pair_3.forward = self.mocked_primer_forward
        self.mocked_pair_3.reverse = self.mocked_primer_reverse

        self.mocked_primer_pairs = [self.mocked_pair_1, self.mocked_pair_2, self.mocked_pair_3]

        self.mocked_pair_with_same_stringency_1 = PrimerPair(
            pair_id="pair4",
            chromosome="ch1",
            pre_targeton_start=11540,
            pre_targeton_end=11545,
            product_size=6,
            stringency=0.4,
            targeton_id="targeton_id",
            uid="uid4")

        self.mocked_pair_with_same_stringency_1.forward = self.mocked_primer_forward
        self.mocked_pair_with_same_stringency_1.reverse = self.mocked_primer_reverse

        self.mocked_pair_with_same_stringency_2 = PrimerPair(
            pair_id="pair5",
            chromosome="ch1",
            pre_targeton_start=11540,
            pre_targeton_end=11545,
            product_size=7,
            stringency=0.4,
            targeton_id="targeton_id",
            uid="uid5")

        self.mocked_pair_with_same_stringency_2.forward = self.mocked_primer_forward
        self.mocked_pair_with_same_stringency_2.reverse = self.mocked_primer_reverse

        self.mocked_primer_data_with_stringency = [self.mocked_pair_with_same_stringency_1, self.mocked_pair_with_same_stringency_2]

        self.mocked_pair_4 = PrimerPair(
            pair_id="pair4",
            chromosome="ch1",
            pre_targeton_start=11540,
            pre_targeton_end=11545,
            product_size=1,
            stringency=0.1,
            targeton_id="targeton_id",
            uid="uid6")

        self.mocked_pair_4.forward = self.mocked_primer_forward
        self.mocked_pair_4.reverse = self.mocked_primer_reverse

        self.mocked_pair_5 = PrimerPair(
            pair_id="pair5",
            chromosome="ch1",
            pre_targeton_start=11540,
            pre_targeton_end=11545,
            product_size=1,
            stringency=0.1,
            targeton_id="targeton_id",
            uid="uid7")

        self.mocked_pair_5.forward = self.mocked_primer_forward
        self.mocked_pair_5.reverse = self.mocked_primer_reverse

        self.mocked_same_stringency_same_product_size = [self.mocked_pair_4 , self.mocked_pair_5]

        self.expected_ranked_primers = """
primer_type,primer,penalty,sequence,primer_start,primer_end,tm,gc_percent,self_any_th,self_end_th,hairpin_th,end_stability,pair_uid,stringency,chromosome,pre_targeton_start,pre_targeton_end,product_size,targeton_id
mockType,primer_fr,0.5,ATCGATCG,11540,11545,60.0,50.0,30.0,10.0,20.0,25.0,uid3,0.3,ch1,11540,11645,3,targeton_id
mockType,primer_rv,0.5,CTCGATCG,11640,11645,60.0,50.0,30.0,10.0,20.0,25.0,uid3,0.3,ch1,11540,11645,3,targeton_id
mockType,primer_fr,0.5,ATCGATCG,11540,11545,60.0,50.0,30.0,10.0,20.0,25.0,uid1,0.1,ch1,11540,11645,3,targeton_id
mockType,primer_rv,0.5,CTCGATCG,11640,11645,60.0,50.0,30.0,10.0,20.0,25.0,uid1,0.1,ch1,11540,11645,3,targeton_id
mockType,primer_fr,0.5,ATCGATCG,11540,11545,60.0,50.0,30.0,10.0,20.0,25.0,uid2,0.2,ch1,11540,11645,1,targeton_id
mockType,primer_rv,0.5,CTCGATCG,11640,11645,60.0,50.0,30.0,10.0,20.0,25.0,uid2,0.2,ch1,11540,11645,1,targeton_id
"""

        self.expected_same_stringency_primer_pairs = """
primer_type,primer,penalty,sequence,primer_start,primer_end,tm,gc_percent,self_any_th,self_end_th,hairpin_th,end_stability,pair_uid,stringency,chromosome,pre_targeton_start,pre_targeton_end,product_size,targeton_id
mockType,primer_fr,0.5,ATCGATCG,11540,11545,60.0,50.0,30.0,10.0,20.0,25.0,uid5,0.4,ch1,11540,11545,7,targeton_id
mockType,primer_rv,0.5,CTCGATCG,11640,11645,60.0,50.0,30.0,10.0,20.0,25.0,uid5,0.4,ch1,11540,11545,7,targeton_id
mockType,primer_fr,0.5,ATCGATCG,11540,11545,60.0,50.0,30.0,10.0,20.0,25.0,uid4,0.4,ch1,11540,11545,6,targeton_id
mockType,primer_rv,0.5,CTCGATCG,11640,11645,60.0,50.0,30.0,10.0,20.0,25.0,uid4,0.4,ch1,11540,11545,6,targeton_id
"""

        self.expected_same_stringency_same_product_size = """
primer_type,primer,penalty,sequence,primer_start,primer_end,tm,gc_percent,self_any_th,self_end_th,hairpin_th,end_stability,pair_uid,stringency,chromosome,pre_targeton_start,pre_targeton_end,product_size,targeton_id
mockType,primer_fr,0.5,ATCGATCG,11540,11545,60.0,50.0,30.0,10.0,20.0,25.0,uid6,0.1,ch1,11540,11545,1,targeton_id
mockType,primer_rv,0.5,CTCGATCG,11640,11645,60.0,50.0,30.0,10.0,20.0,25.0,uid6,0.1,ch1,11540,11545,1,targeton_id
mockType,primer_fr,0.5,ATCGATCG,11540,11545,60.0,50.0,30.0,10.0,20.0,25.0,uid7,0.1,ch1,11540,11545,1,targeton_id
mockType,primer_rv,0.5,CTCGATCG,11640,11645,60.0,50.0,30.0,10.0,20.0,25.0,uid7,0.1,ch1,11540,11545,1,targeton_id
"""
        self.expected_unranked = """
primer_type,primer,penalty,sequence,primer_start,primer_end,tm,gc_percent,self_any_th,self_end_th,hairpin_th,end_stability,pair_uid,stringency,chromosome,pre_targeton_start,pre_targeton_end,product_size,targeton_id
mockType,primer_fr,0.5,ATCGATCG,11540,11545,60.0,50.0,30.0,10.0,20.0,25.0,uid1,0.1,ch1,11540,11645,3,targeton_id
mockType,primer_rv,0.5,CTCGATCG,11640,11645,60.0,50.0,30.0,10.0,20.0,25.0,uid1,0.1,ch1,11540,11645,3,targeton_id
mockType,primer_fr,0.5,ATCGATCG,11540,11545,60.0,50.0,30.0,10.0,20.0,25.0,uid2,0.2,ch1,11540,11645,1,targeton_id
mockType,primer_rv,0.5,CTCGATCG,11640,11645,60.0,50.0,30.0,10.0,20.0,25.0,uid2,0.2,ch1,11540,11645,1,targeton_id
mockType,primer_fr,0.5,ATCGATCG,11540,11545,60.0,50.0,30.0,10.0,20.0,25.0,uid3,0.3,ch1,11540,11645,3,targeton_id
mockType,primer_rv,0.5,CTCGATCG,11640,11645,60.0,50.0,30.0,10.0,20.0,25.0,uid3,0.3,ch1,11540,11645,3,targeton_id
"""
        self.expected_product_size_true = """
primer_type,primer,penalty,sequence,primer_start,primer_end,tm,gc_percent,self_any_th,self_end_th,hairpin_th,end_stability,pair_uid,stringency,chromosome,pre_targeton_start,pre_targeton_end,product_size,targeton_id
mockType,primer_fr,0.5,ATCGATCG,11540,11545,60.0,50.0,30.0,10.0,20.0,25.0,uid1,0.1,ch1,11540,11645,3,targeton_id
mockType,primer_rv,0.5,CTCGATCG,11640,11645,60.0,50.0,30.0,10.0,20.0,25.0,uid1,0.1,ch1,11540,11645,3,targeton_id
mockType,primer_fr,0.5,ATCGATCG,11540,11545,60.0,50.0,30.0,10.0,20.0,25.0,uid3,0.3,ch1,11540,11645,3,targeton_id
mockType,primer_rv,0.5,CTCGATCG,11640,11645,60.0,50.0,30.0,10.0,20.0,25.0,uid3,0.3,ch1,11540,11645,3,targeton_id
mockType,primer_fr,0.5,ATCGATCG,11540,11545,60.0,50.0,30.0,10.0,20.0,25.0,uid2,0.2,ch1,11540,11645,1,targeton_id
mockType,primer_rv,0.5,CTCGATCG,11640,11645,60.0,50.0,30.0,10.0,20.0,25.0,uid2,0.2,ch1,11540,11645,1,targeton_id
"""


    def tearDown(self):
        # Remove the handler after each test to reset logging
        logger = logging.getLogger()
        logger.removeHandler(self.handler)


    def test_ranker_all_true(self):

        expected_ranked_df = _create_dataframe(self.expected_ranked_primers)

        ranked_df = Ranker(self.mocked_ranking_config_all_true).rank("mockType", self.mocked_primer_pairs)

        # index being dropped to make comparison as ranking does not reset the index
        assert_frame_equal(ranked_df.reset_index(drop=True), expected_ranked_df)


    def test_ranker_with_same_stringency_rank_by_product_size_all_true(self):
        expected_ranked_with_same_stringency_df = _create_dataframe(self.expected_same_stringency_primer_pairs)

        ranked_df = Ranker(self.mocked_ranking_config_all_true).rank("mockType",  self.mocked_primer_data_with_stringency)

        # index being dropped to make comparison as ranking does not reset the index
        assert_frame_equal(ranked_df.reset_index(drop=True), expected_ranked_with_same_stringency_df)


    def test_ranker_with_same_stringency_same_product_size_all_true(self):
        expected_ranked_with_same_stringency_same_product_size_df = _create_dataframe(self.expected_same_stringency_same_product_size)

        unranked_df = Ranker(self.mocked_ranking_config_all_true).rank("mockType", self.mocked_same_stringency_same_product_size)
        # designer_config.params['ranking']

        # index being dropped to make comparison as ranking does not reset the index
        assert_frame_equal(unranked_df.reset_index(drop=True), expected_ranked_with_same_stringency_same_product_size_df)


    def test_empty_df_when_no_primer_pairs_provided(self):

        ranked_df = Ranker(self.mocked_ranking_config_all_true).rank("mockType", [])

        self.assertTrue(ranked_df.empty)

    def test_ranker_when_config_empty(self):

        # Arrange
        expected_unranked_df = _create_dataframe(self.expected_unranked)
        expected_warning = "No ranking criteria defined in config file under the ranking key"\
                           "\nNo ranking applied"

        # Act
        result = Ranker(self.mocked_ranking_config_empty).rank("mockType",
                                                                        self.mocked_primer_pairs)
        logs = self.handler.buffer.getvalue().strip()

        # Assert
        self.assertEqual(logs, expected_warning)
        assert_frame_equal(result.reset_index(drop = True), expected_unranked_df)


    def test_ranker_when_wrong_key(self):

        # Arrange
        expected_error = "'product_size' missing in config file for ranking - Will not be " \
                         "used for ranking\nInvalid name(s) provided for ranking in config " \
                         "file: 'product_sizes'. The only valid names are: stringency, " \
                         "product_size. Unable to apply ranking - Exiting programme"

        # Act
        with self.assertRaises(SystemExit):
            Ranker(self.mocked_ranking_wrong_key).rank("mockType", self.mocked_primer_pairs)
        logs = self.handler.buffer.getvalue().strip()

        # Assert
        self.assertEqual(logs, expected_error)


    def test_ranker_when_wrong_value(self):

        # Arrange
        expected_error = "Wrong value(s) provided for 'product_size' in config file (only takes " \
                         "true or false). Unable to apply ranking - Exiting programme"

        # Act
        with self.assertRaises(SystemExit):
            Ranker(self.mocked_ranking_wrong_value).rank("mockType", self.mocked_primer_pairs)
        logs = self.handler.buffer.getvalue().strip()

        # Assert
        self.assertEqual(logs, expected_error)


    def test_ranker_when_all_false(self):

        # Arrange
        expected_unranked_df = _create_dataframe(self.expected_unranked)

        # Act
        result = Ranker(self.mocked_ranking_config_all_false).rank("mockType",
                                                                   self.mocked_primer_pairs)

        # Assert
        assert_frame_equal(result.reset_index(drop = True), expected_unranked_df)


    def test_ranker_when_one_true(self):

        # Arrange
        expected_ranked_df = _create_dataframe(self.expected_product_size_true)

        # Act
        result = Ranker(self.mocked_ranking_config_one_true).rank("mockType",
                                                                  self.mocked_primer_pairs)

        # Assert
        assert_frame_equal(result.reset_index(drop = True), expected_ranked_df)


def _create_dataframe(data: str) -> pd.DataFrame:
    data_io = StringIO(data)
    df = pd.read_csv(data_io, delim_whitespace=False)
    df["chromosome"] = df["chromosome"].astype(str)
    return df
