from io import StringIO
from unittest import TestCase
from unittest.mock import Mock
import pandas as pd
from ranking.ranker import Ranker
from primer.primer_pair import PrimerPair
from primer.designed_primer import DesignedPrimer, Interval

class RankerTest(TestCase):

    def setUp(self):

        self.mocked_primer = DesignedPrimer(
            name="primer_with_variant",
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

        self.mocked_pair_1 = PrimerPair(
            pair_id="pair1",
            chromosome="1",
            pre_targeton_start=11540,
            pre_targeton_end=11545,
            product_size=1,
            stringency=0.1,
            targeton_id="targeton_id",
            uid="uid")
        
        self.mocked_pair_1.forward = self.mocked_primer
        self.mocked_pair_1.reverse = self.mocked_primer

        self.mocked_pair_2 = PrimerPair(
            pair_id="pair2",
            chromosome="1",
            pre_targeton_start=11540,
            pre_targeton_end=11545,
            product_size=2,
            stringency=0.2,
            targeton_id="targeton_id",
            uid="uid")
        
        self.mocked_pair_2.forward = self.mocked_primer
        self.mocked_pair_2.reverse = self.mocked_primer

        self.mocked_pair_3 = PrimerPair(
            pair_id="pair3",
            chromosome="1",
            pre_targeton_start=11540,
            pre_targeton_end=11545,
            product_size=3,
            stringency=0.3,
            targeton_id="targeton_id",
            uid="uid")

        self.mocked_pair_3.forward = self.mocked_primer
        self.mocked_pair_3.reverse = self.mocked_primer

        self.mocked_primer_pairs = [self.mocked_pair_1, self.mocked_pair_2, self.mocked_pair_3]

        self.expected_ranked_primers = """primer_type,primer,penalty,sequence,primer_start,primer_end,tm,gc_percent,self_any_th,self_end_th,hairpin_th,end_stability,pair_uid,stringency,chromosome,pre_targeton_start,pre_targeton_end,product_size,targeton_id
mockType,primer_with_variant,0.5,ATCGATCG,11540,11545,60.0,50.0,30.0,10.0,20.0,25.0,uid,0.3,1,11540,11545,3,targeton_id
mockType,primer_with_variant,0.5,ATCGATCG,11540,11545,60.0,50.0,30.0,10.0,20.0,25.0,uid,0.3,1,11540,11545,3,targeton_id
mockType,primer_with_variant,0.5,ATCGATCG,11540,11545,60.0,50.0,30.0,10.0,20.0,25.0,uid,0.2,1,11540,11545,2,targeton_id
mockType,primer_with_variant,0.5,ATCGATCG,11540,11545,60.0,50.0,30.0,10.0,20.0,25.0,uid,0.2,1,11540,11545,2,targeton_id
mockType,primer_with_variant,0.5,ATCGATCG,11540,11545,60.0,50.0,30.0,10.0,20.0,25.0,uid,0.1,1,11540,11545,1,targeton_id
mockType,primer_with_variant,0.5,ATCGATCG,11540,11545,60.0,50.0,30.0,10.0,20.0,25.0,uid,0.1,1,11540,11545,1,targeton_id"""


    def test_ranker(self):

        expected_ranked_df = _create_dataframe(self.expected_ranked_primers)

        ranked_df = Ranker().rank("mockType", self.mocked_primer_pairs)

        ranked_df.reset_index(drop=True).equals(expected_ranked_df.astype(object).reset_index(drop=True))
    
    def test_empty_df_when_no_primer_pairs_provided(self):

        ranked_df = Ranker().rank("mockType", [])

        self.assertEqual(len(ranked_df), 0)
        

def _config(ranking_priority_order: list) -> Mock:
    config = Mock()
    config.ranking_priority_order = ranking_priority_order
    return config

def _create_dataframe(data: str) -> pd.DataFrame:
    data_io = StringIO(data)
    df = pd.read_csv(data_io, delim_whitespace=False)
    df["chromosome"] = df["chromosome"].astype(int) 
    return df
