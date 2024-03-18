from pathlib import Path
from abc import ABC, abstractmethod

from utils.file_system import parse_json
from utils.exceptions import InvalidConfigError


class Config(ABC):
    def __init__(self):
        pass

    @abstractmethod
    def read_config(self):
        pass
        

class DesignerConfig(Config):
    def __init__(self, config_file: str):
        self._default_config_file = \
            Path(__file__).parent / '../../config/designer.config.json'

        config = self.read_config(self._default_config_file, config_file)
        self.stringency_vector = config['stringency_vector']

    @staticmethod
    def read_config(
            default_config_file: str,
            config_file: str = '',
    ) -> dict:
        default_config = parse_json(default_config_file)
        keys = default_config.keys()

        if config_file == '' or config_file == default_config_file:
            return default_config
        else:
            config = parse_json(config_file)

            for field in keys:
                config.setdefault(field, default_config[field])

        return config


class Primer3ParamsConfig(Config):
    def __init__(self, config_file: str):
        self._default_config_file = \
            Path(__file__).parent / 'primer' / 'primer3.config.json'

        self.params = self.read_config(self._default_config_file, config_file)

    @staticmethod
    def read_config(
            default_config_file: str,
            config_file: str = '',
    ) -> dict:
        file = config_file if config_file != '' else default_config_file

        try:
            config_data = parse_json(file)
        except Exception:
            raise InvalidConfigError(
                'Primer3 config file is not a correct JSON')

        return config_data





