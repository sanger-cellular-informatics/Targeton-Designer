from utils.file_system import parse_json

import sys
from custom_logger.custom_logger import CustomLogger

# Initialize logger
logger = CustomLogger(__name__)


class DesignerConfig:
    def __init__(self, args: dict):
        default_config_path = 'config/default_designer.config.json'

        config = DesignerConfig.read_config(default_config_path, args.get('conf', None))

        # Check if filters exist in configuration.
        if not config.get("filters"):
            # If duplicates is False by mistake, it will be enabled to run as a default filter.
            # Because we want to run duplicates filter even there's no filter added in configuration.
            config["filters"] = {"duplicates": True}

        # Check ranking
        if not config.get("ranking"):
            config["ranking"] = {}

        # Check primer region restriction parameters
        self.flanking_region = config['flanking_region'] or 0
        self.exclusion_region = config['exclusion_region'] or 0
        if not isinstance(self.flanking_region, int) or self.flanking_region < 0:
            logger.error("flanking_region must be a non-negative integer")
            sys.exit(1)
        if not self.flanking_region:
            logger.info(("flanking_region set to 0, so primer placement will not be restricted the flanking_region, "
                        "and exclusion_region will be ignored."))
        if not isinstance(self.exclusion_region, int) or self.exclusion_region < 0:
            logger.error("exclusion_region must be a non-negative integer")
            sys.exit(1)

        # Check Ipcress output file parameters
        ipcress_params = config.get("ipcress_parameters") or {}
        self.ipcress_params_write_file = False

        if ipcress_params.get("write_ipcress_file"):
            if not isinstance(ipcress_params.get("min_size"), int) or ipcress_params.get("min_size") < 0:
                logger.error("Designer config error: ipcress_parameters.min_size must be a non-negative integer")
            if not isinstance(ipcress_params.get("max_size"), int) or ipcress_params.get("max_size") < 0:
                logger.error("Designer config error: ipcress_parameters.max_size must be a non-negative integer")
  
            self.ipcress_params_write_file = ipcress_params["write_ipcress_file"] 
            self.ipcress_params_min_size = ipcress_params.get("min_size")
            self.ipcress_params_max_size = ipcress_params.get("max_size")

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
            default_config_path: str,
            config_path: str = None,
    ) -> dict:
        default_config = parse_json(default_config_path)
        keys = default_config.keys()

        if config_path is None or config_path == default_config_path:
            return default_config
        else:
            config = parse_json(config_path)

            for field in keys:
                config.setdefault(field, default_config[field])

        return config
