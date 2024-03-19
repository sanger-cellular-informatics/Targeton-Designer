from typing import List
from primer.primer_pair import PrimerPair


def rename_primers(pairs: List[PrimerPair], primer_type: str = 'LibAmp'):
    for index, pair in enumerate(pairs):
        slice_name = 'slice_name'

        pair.forward['primer'] = slice_name + '_' + primer_type + 'F_' + str(index)
        pair.reverse['primer'] = slice_name + '_' + primer_type + 'R_' + str(index)

    return pairs
