class IpcressParameters:
    def __init__(self, params: dict):
        if not isinstance(params, dict):
            raise TypeError(
                f"ipcress_parameters must be a dict, got {type(params).__name__}"
            )

        self.write_ipcress_file = self._get_bool(
            params, "write_ipcress_file"
        )
        self.min_size = self._get_non_negative_int(
            params, "min_size"
        )
        self.max_size = self._get_non_negative_int(
            params, "max_size"
        )

    @staticmethod
    def _get_bool(params: dict, key: str) -> bool:
        value = params.get(key)
        if not isinstance(value, bool):
            raise ValueError(
                f"ipcress_parameters.{key} must be a boolean (true/false), "
                f"got {type(value).__name__}"
            )
        return value

    @staticmethod
    def _get_non_negative_int(params: dict, key: str) -> int:
        value = params.get(key)
        if not isinstance(value, int) or value < 0:
            raise ValueError(
                f"ipcress_parameters.{key} must be a non-negative integer, "
                f"got {value}"
            )
        return value
