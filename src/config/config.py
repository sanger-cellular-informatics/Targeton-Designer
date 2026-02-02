from math import log
from config.config_helpers import check_non_negative_integer
from config.ipcress_params import IpcressParameters
from utils.file_system import parse_json

import sys
from custom_logger.custom_logger import CustomLogger
import os
from pathlib import Path

# Initialize logger
logger = CustomLogger(__name__)


BASE_DIR = Path(__file__).resolve().parent.parent.parent
DEFAULT_CONFIG_PATH = BASE_DIR / "config" / "default_designer.config.json"
DEFAULT_PRIMER3_CONFIG_PATH = BASE_DIR / "config" / "default_primer3.config.json"

class DesignerConfig:
    def __init__(self, args: dict):


        designer_config_path = args.get('conf')

        logger.info(f"Base directory for configuration: {BASE_DIR}")
        logger.info(f"Using configuration file: {designer_config_path or DEFAULT_CONFIG_PATH}")

        config = DesignerConfig.read_config(DEFAULT_CONFIG_PATH, designer_config_path)

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

        check_non_negative_integer(name="flanking_region", value=self.flanking_region)
        if not self.flanking_region:
            logger.info(("flanking_region set to 0, so primer placement will not be restricted the flanking_region, "
                         "and exclusion_region will be ignored."))
        check_non_negative_integer(name="exclusion_region", value=self.exclusion_region)

        ipcress_block = config.get("ipcress_parameters")
        if ipcress_block is not None:
            self.ipcress_params = IpcressParameters(ipcress_block)
        else:
            self.ipcress_params = None

        self.stringency_vector = config['stringency_vector']
        self.csv_column_order = config['csv_column_order']
        self.filters = config['filters']
        self.ranking = config['ranking']

        self.prefix_output_dir = args.get('dir', None) or config.get('dir', None)
        self.fasta = args.get('fasta', None) or config.get('fasta', None)

        # Use provided primer3_params path from args, otherwise use default primer3 config path
        primer3_params_path = (args.get('primer3_params', None) or config.get('primer3_params', None)
                               or DEFAULT_PRIMER3_CONFIG_PATH)

        logger.info(f"Using primer3 parameters file: {primer3_params_path}")

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
