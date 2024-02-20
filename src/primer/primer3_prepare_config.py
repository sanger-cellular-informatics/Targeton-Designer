import copy

def prepare_config(base_config: dict, stringency: str) -> dict:
    result_config = copy.deepcopy(base_config)

    if result_config.get('PRIMER_MASK_FAILURE_RATE') is not None:
        result_config['PRIMER_MASK_FAILURE_RATE'] = stringency

    return result_config