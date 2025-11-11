# config_manager.py
from pyCATHY.version_config import CONFIG_MAP


class VersionConfigManager:
    """
    Generic version-dependent configuration manager that uses modular dictionaries.

    Example:
        cfg = VersionConfigManager("prepro", "update_prepo_inputs", version="withIrr")
        params = cfg.get("terrain_parameter")
    """

    def __init__(self, module: str, function: str, version: str = "default"):
        self.module = module
        self.function = function
        self.version = version

        # Find the dictionary for this module and function
        module_dict = CONFIG_MAP.get(module, {})
        function_dict = module_dict.get(function, {})

        # Select the version configuration (or default)
        self.version_config = function_dict.get(version, function_dict.get("default", {}))

    def get(self, key: str, fallback=None):
        """Return a specific parameter or a fallback."""
        return self.version_config.get(key, fallback)

    def all(self):
        """Return the full configuration dictionary for this function/version."""
        return self.version_config

    def list_versions(self):
        """Return all available versions for this function."""
        module_dict = CONFIG_MAP.get(self.module, {})
        function_dict = module_dict.get(self.function, {})
        return list(function_dict.keys())

    def __repr__(self):
        return f"<VersionConfigManager module='{self.module}' function='{self.function}' version='{self.version}'>"
