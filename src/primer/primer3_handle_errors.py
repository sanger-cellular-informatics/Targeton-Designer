import primer3

from typing import List

from utils.exceptions import Primer3Error
from primer.primer_pair import PrimerPair
from custom_logger.custom_logger import CustomLogger

# Initialize logger
logger = CustomLogger(__name__)


def format_no_primer_pairs_message(stringency: int,
                                   primer_explain_flag: int,
                                   designs: dict) -> str:
    msg = 'Stringency level ' + str(stringency) + " -- "
    if primer_explain_flag:
        msg += _get_primer3_explain(designs)
    else:
        msg += 'No primer pairs returned; add PRIMER_EXPLAIN_FLAG == 1 to config file for more details'
    return msg


def handle_primer3_errors(primer_explain: List[str], primer_pairs: List[PrimerPair]) -> None:
    message = '\n'.join([msg for msg in primer_explain])

    if not primer_pairs:
        message = 'NO PRIMER PAIRS BUILT BY PRIMER3:\n' + message
        logger.exception(Primer3Error(message))
        raise Primer3Error(message)

    else:
        message = 'Warning: No primer pairs built by Primer3 with the following stringencies:\n' + message
        logger.warning(message)


def _get_primer3_explain(designs: dict) -> str:
    keys = ["PRIMER_LEFT_EXPLAIN",
            "PRIMER_RIGHT_EXPLAIN",
            "PRIMER_PAIR_EXPLAIN"]
    msg = dict((key, designs[key]) for key in keys)
    msg_formatted = '; '.join([f"{key}: {value}" for key, value in msg.items()])
    return msg_formatted
