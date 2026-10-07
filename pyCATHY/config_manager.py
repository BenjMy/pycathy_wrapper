"""
config_manager.py

Small helper that resolves version-specific configuration dictionaries
for pyCATHY, backed by version_config.py's CONFIG_MAP.

This module did not exist yet in the installed package -- cathy_tools.py
(and the experimental cathy_tools_NEWBRANCH.py) import
`VersionConfigManager` from here, but nothing had actually implemented
it. This is that implementation, matching the interface already in use
(`.all()`, `.get(key, default)`) so no other file needs to change.

Expected location: pyCATHY/config_manager.py, alongside
pyCATHY/version_config.py (see the commented-out import in
cathy_tools_NEWBRANCH.py: "# from pyCATHY.version_config import ...").
"""

from pyCATHY.version_config import CONFIG_MAP


class VersionConfigManager:
    """
    Resolve a version-specific configuration dict from CONFIG_MAP.

    Parameters
    ----------
    section : str
        Top-level CONFIG_MAP key, e.g. "cathy" or "prepro".
    key : str
        Second-level key identifying which config table to use, e.g.
        "update_cathyH", "header", "soil_format", "update_veg_map".
    version : str
        The CATHY object's self.version (e.g. "SCF_variable",
        "withIrr", "1.0.0", or anything else). If this exact version
        name isn't a key in the table, falls back to "default" -- so an
        unrecognized or unset version silently gets default behavior
        instead of raising, which matters for backward compatibility:
        existing scripts that never pass a special version keep working
        unchanged.

    Raises
    ------
    KeyError
        If `section` or `key` themselves are not found in CONFIG_MAP --
        that indicates a real configuration/programming error (a typo in
        the section/key name), not just an unconfigured version, so it
        is not silently swallowed the way a missing version is.

    Examples
    --------
    >>> cfg = VersionConfigManager("cathy", "soil_format", "SCF_variable")
    >>> cfg.get("scf_per_veg", False)
    True
    >>> cfg = VersionConfigManager("cathy", "soil_format", "some_other_version")
    >>> cfg.get("scf_per_veg", False)
    False
    """

    def __init__(self, section, key, version):
        try:
            section_map = CONFIG_MAP[section]
        except KeyError:
            raise KeyError(
                f"VersionConfigManager: unknown section '{section}' -- "
                f"available sections are {list(CONFIG_MAP.keys())}"
            )
        try:
            table = section_map[key]
        except KeyError:
            raise KeyError(
                f"VersionConfigManager: unknown key '{key}' for section "
                f"'{section}' -- available keys are {list(section_map.keys())}"
            )

        self.section = section
        self.key = key
        self.version = version
        self._cfg = table.get(version, table.get("default", {}))

    def all(self):
        """Return the full resolved configuration dict for this version."""
        return self._cfg

    def get(self, name, default=None):
        """Return a single value from the resolved config, with a fallback."""
        return self._cfg.get(name, default)

    def __repr__(self):
        return (
            f"VersionConfigManager(section={self.section!r}, "
            f"key={self.key!r}, version={self.version!r}, "
            f"resolved={self._cfg!r})"
        )
