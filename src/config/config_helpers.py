from custom_logger.custom_logger import CustomLogger

import sys

# Initialize logger
logger = CustomLogger(__name__)

def check_non_negative_integer(name, value):
    if not isinstance(value, int) or value < 0:
        raise ValueError(
             f'Invalid config value for "{name}": expected non-negative int, '
            f'got {value} ({type(value).__name__})'
        )

def check_boolean(name, value):
    if not isinstance(value, bool) or value < 0:
        raise ValueError(
             f'Invalid config value for "{name}": expected boolean, '
            f'got {value} ({type(value).__name__})'
        )
