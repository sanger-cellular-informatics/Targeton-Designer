from abc import ABC, abstractmethod
from utils.file_system import parse_json

from custom_logger.custom_logger import CustomLogger

# Initialize logger
logger = CustomLogger(__name__)


class Config(ABC):
    def __init__(self):
        pass

    @abstractmethod
    def read_config(self):
        pass


class DesignerConfig(Config):
    def __init__(self, config_file: str = None):

        self._default_config_file = 'config/designer.config.json'

        config = self.read_config(self._default_config_file, config_file)

        # Check if filters exist in configuration.
        if not config.get("filters"):
            # If duplicates is False by mistake, it will be enabled to run as a default filter.
            # Because we want to run duplicates filter even there's no filter added in configuration.
            config["filters"] = {"duplicates": True}

        if not config.get("ranking"):
            config["ranking"] = {}

        self.params = {'stringency_vector': config['stringency_vector'],
                       'csv_column_order': config['csv_column_order'],
                       'filters': config['filters'],
                       'ranking': config['ranking']}


    @staticmethod
    def read_config(
            default_config_file: str,
            config_file: str = None,
    ) -> dict:
        default_config = parse_json(default_config_file)
        keys = default_config.keys()

        if config_file is None or config_file == default_config_file:
            return default_config
        else:
            config = parse_json(config_file)

            for field in keys:
                config.setdefault(field, default_config[field])

        return config


class Primer3ParamsConfig(Config):
    def __init__(self, config_file: str = None):
        self._default_config_file = 'src/primer/primer3.config.json'
        self.params = self.read_config(self._default_config_file, config_file)

    @staticmethod
    def read_config(
            default_config_file: str,
            config_file: str = None,
    ) -> dict:

        file = default_config_file
        if config_file is not None:
            file = config_file

        config_data = parse_json(file)

        return config_data
