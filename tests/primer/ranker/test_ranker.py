import logging
from unittest import TestCase

from primer.ranker.ranker import Ranker
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

    def tearDown(self):
        # Remove the handler after each test to reset logging
        logger = logging.getLogger()
        logger.removeHandler(self.handler)


    def test_ranker_all_true(self):

        ranked_pairs = Ranker(self.mocked_ranking_config_all_true).rank(
            self.mocked_primer_pairs,
        )

        _assert_pairs_by_id_stringency_product_size(
            self,
            ranked_pairs,
            [
                ("pair1", 0.1, 3),
                ("pair3", 0.3, 3),
                ("pair2", 0.2, 1),
            ],
        )


    def test_ranker_with_same_stringency_rank_by_product_size_all_true(self):
        ranked_pairs = Ranker(self.mocked_ranking_config_all_true).rank(
            self.mocked_primer_data_with_stringency,
        )

        _assert_pairs_by_id_stringency_product_size(
            self,
            ranked_pairs,
            [
                ("pair5", 0.4, 7),
                ("pair4", 0.4, 6),
            ],
        )


    def test_ranker_with_same_stringency_same_product_size_all_true(self):
        ranked_pairs = Ranker(self.mocked_ranking_config_all_true).rank(
            self.mocked_same_stringency_same_product_size,
        )

        _assert_pairs_by_id_stringency_product_size(
            self,
            ranked_pairs,
            [
                ("pair4", 0.1, 1),
                ("pair5", 0.1, 1),
            ],
        )


    def test_empty_df_when_no_primer_pairs_provided(self):

        ranked_pairs = Ranker(self.mocked_ranking_config_all_true).rank([])

        self.assertEqual(ranked_pairs, [])

    def test_ranker_when_config_empty(self):

        # Arrange
        expected_warning = "No ranking criteria defined in config file under the ranking key"\
                           "\nNo ranking applied"

        # Act
        result = Ranker(self.mocked_ranking_config_empty).rank(self.mocked_primer_pairs)
        logs = self.handler.buffer.getvalue().strip()

        # Assert
        self.assertEqual(logs, expected_warning)
        _assert_pairs_by_id_stringency_product_size(
            self,
            result,
            [
                ("pair1", 0.1, 3),
                ("pair2", 0.2, 1),
                ("pair3", 0.3, 3),
            ],
        )


    def test_ranker_when_wrong_key(self):

        # Arrange
        expected_error = "'product_size' missing in config file for ranking - Will not be " \
                         "used for ranking\nInvalid name(s) provided for ranking in config " \
                         "file: 'product_sizes'. The only valid names are: stringency, " \
                         "product_size. Unable to apply ranking - Exiting programme"

        # Act
        with self.assertRaises(SystemExit):
            Ranker(self.mocked_ranking_wrong_key).rank(self.mocked_primer_pairs)
        logs = self.handler.buffer.getvalue().strip()

        # Assert
        self.assertEqual(logs, expected_error)


    def test_ranker_when_wrong_value(self):

        # Arrange
        expected_error = "Wrong value(s) provided for 'product_size' in config file (only takes " \
                         "true or false). Unable to apply ranking - Exiting programme"

        # Act
        with self.assertRaises(SystemExit):
            Ranker(self.mocked_ranking_wrong_value).rank(self.mocked_primer_pairs)
        logs = self.handler.buffer.getvalue().strip()

        # Assert
        self.assertEqual(logs, expected_error)


    def test_ranker_when_all_false(self):

        # Act
        result = Ranker(self.mocked_ranking_config_all_false).rank(self.mocked_primer_pairs)

        # Assert
        _assert_pairs_by_id_stringency_product_size(
            self,
            result,
            [
                ("pair1", 0.1, 3),
                ("pair2", 0.2, 1),
                ("pair3", 0.3, 3),
            ],
        )

    def test_ranker_when_one_true(self):
        # Act
        result = Ranker(self.mocked_ranking_config_one_true).rank(self.mocked_primer_pairs)

        # Assert
        _assert_pairs_by_id_stringency_product_size(
            self,
            result,
            [
                ("pair1", 0.1, 3),
                ("pair3", 0.3, 3),
                ("pair2", 0.2, 1),
            ],
        )
        

def _assert_pairs_by_id_stringency_product_size(
    test_case: TestCase,
    pairs: list,
    expected: list,
) -> None:
    test_case.assertEqual(
        [(pair.id, pair.stringency, pair.product_size) for pair in pairs],
        expected,
    )
