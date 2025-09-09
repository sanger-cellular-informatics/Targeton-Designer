from utils.file_system import parse_json

from custom_logger.custom_logger import CustomLogger

# Initialize logger
logger = CustomLogger(__name__)


class DesignerConfig:
    def __init__(self, args: dict):
        default_config_file = 'config/default_designer.config.json'

        config = DesignerConfig.read_config(default_config_file, args.get('conf', None))

        # Check if filters exist in configuration.
        if not config.get("filters"):
            # If duplicates is False by mistake, it will be enabled to run as a default filter.
            # Because we want to run duplicates filter even there's no filter added in configuration.
            config["filters"] = {"duplicates": True}

        if not config.get("ranking"):
            config["ranking"] = {}

        self.region_padding = config['region_padding']
        self.stringency_vector = config['stringency_vector']
        self.csv_column_order = config['csv_column_order']
        self.filters = config['filters']
        self.ranking = config['ranking']

        self.prefix_output_dir = args.get('dir', None) or config.get('dir', None)
        self.fasta = args.get('fasta', None) or config.get('fasta', None)

        primer3_params_path = (args.get('primer3_params', None) or config.get('primer3_params', None)
                               or 'config/default_primer3.config.json')
        self.primer3_params = parse_json(primer3_params_path)

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
