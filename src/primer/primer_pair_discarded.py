from primer.primer_pair import PrimerPair


class PrimerPairDiscarded(PrimerPair):
    def __init__(self, primer_pair: PrimerPair, reason_discarded: str):
        super().__init__(pair_id=primer_pair.id,
                         chromosome=primer_pair.chromosome,
                         pre_targeton_start=primer_pair.pre_targeton_start,
                         pre_targeton_end=primer_pair.pre_targeton_end,
                         product_size=primer_pair.product_size,
                         stringency=primer_pair.stringency,
                         targeton_id=primer_pair.targeton_id,
                         uid=primer_pair.uid,
                         )
        self.forward = primer_pair.forward
        self.reverse = primer_pair.reverse
        self.reason_discarded = reason_discarded
