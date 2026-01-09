from config.config_helpers import check_non_negative_integer, check_boolean


class IpcressParameters:
    def __init__(self, params: dict):
        if not isinstance(params, dict):
            raise TypeError(
                f"ipcress_parameters must be a dict, got {type(params).__name__}"
            )

        self.write_ipcress_file = self._get_bool(params, "write_ipcress_file")
        self.min_size = self._get_non_negative_int(params, "min_size", )
        self.max_size = self._get_non_negative_int(params, "max_size", )

    def _get_bool(self, params: dict, key: str) -> bool:
        value = params.get(key)
        check_boolean(key, value)
        return value

    def _get_non_negative_int(self, params: dict, key: str) -> int:
        value = params.get(key)
        check_non_negative_integer(key, value)
        return value
