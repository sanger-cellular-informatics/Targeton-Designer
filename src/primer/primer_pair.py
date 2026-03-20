from utils.get_data.hap1 import contain_variant


class PrimerPair:
    def __init__(self, pair_id: str,
                       chromosome: str,
                       pre_targeton_start: int,
                       pre_targeton_end: int,
                       product_size: int,
                       stringency: float,
                       targeton_id: str,
                       uid: str):
        self.id = pair_id
        self.uid = uid
        self.chromosome = chromosome
        self.pre_targeton_start = pre_targeton_start
        self.pre_targeton_end = pre_targeton_end
        self.product_size = product_size
        self.stringency = stringency
        self.targeton_id = targeton_id
        self.reverse = None
        self.forward = None


    def __repr__(self):
        return (f"PrimerPair(pair_id='{self.id}', "
                f"uid='{self.uid}', "
                f"chromosome='{self.chromosome}', "
                f"pre_targeton_start='{self.pre_targeton_start}', "
                f"pre_targeton_end='{self.pre_targeton_end}', "
                f"product_size='{self.product_size}', "
                f"stringency='{self.stringency}',"
                f"targeton_id='{self.targeton_id}', "
                f"forward={self.forward}, "
                f"reverse={self.reverse})"
                )

    def __eq__(self, other):
        if isinstance(other, PrimerPair):
            return (
                    self.chromosome == other.chromosome and
                    self.forward == other.forward and
                    self.reverse == other.reverse
            )
        return False

    def __hash__(self):
        return hash((self.chromosome, self.forward, self.reverse))

    @property
    def contain_hap_one_variant(self) -> bool:
        forward_start, forward_end = self.forward.primer_start, self.forward.primer_end
        reverse_start, reverse_end = self.reverse.primer_start, self.reverse.primer_end

        return (contain_variant(self.chromosome, forward_start, forward_end) or
                contain_variant(self.chromosome, reverse_start, reverse_end))
