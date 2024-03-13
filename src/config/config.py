from pathlib import Path

from utils.arguments_parser import ParsedInputArguments
from utils.file_system import parse_json


DEFAULT_CONFIG_FILE = Path(__file__).parent / '../../config/designer.config.json'

class Config:
    def __init__(self, input_args: ParsedInputArguments):
        args = input_args.arguments
        config = prepare_config(args.get('conf'))

        self.stringency_vector = config['stringency_vector']

def prepare_config(
        config_file: str = None,
        default_config_file: str = DEFAULT_CONFIG_FILE
) -> dict:
    default_config = parse_json(default_config_file)
    keys = default_config.keys()

    print('CONFIG', config_file)

    if config_file is None or config_file == default_config_file:
        return default_config
    else:
        config = parse_json(config_file)

        for field in keys:
            config.setdefault(field, default_config[field])

    return config

