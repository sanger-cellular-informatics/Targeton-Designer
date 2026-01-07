from custom_logger.custom_logger import CustomLogger

import sys

# Initialize logger
logger = CustomLogger(__name__)

def check_non_negative_integer(name, value):
    if not isinstance(value, int) or value < 0:
        logger.error(
            f'Invalid config value for "{name}": expected non-negative int, '
            f'got {value!r} ({type(value).__name__})'
        )
        sys.exit(1)