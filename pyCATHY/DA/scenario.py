"""Validation of a perturbation scenario (``scenario.toml``), before any ensemble is built.

``perturbate.perturbate`` ignores what it does not know: a misspelled parameter name
creates no perturbation, a bad bound leaves a parameter unbounded, a list one entry too
short fails halfway through. ``validate_scenario`` reports every such problem at once.

Only definite mistakes are errors; things that are legal but suspicious are warnings.
"""
import difflib
import re
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

# parameters perturbate.perturbate knows how to build
KNOWN_PARAMETERS = (
    "ic", "WTPOSITION", "atmbc", "ZROOT", "PCREF", "PCWLT", "OMGC", "PZ",
    "Ks", "porosity", "n_VG", "alpha_VG", "thetar_VG", "SATCELL_VG",
)
# named in the scenario format but not built: perturbate() would silently ignore them
NOT_IMPLEMENTED = {
    "VGP": "its perturbation is commented out in perturbate.py, nothing would be built",
    "alpha_VG": (
        "it can be sampled, but writing it into the member soil files fails with 'Unknown "
        "parameter' (no soil column is mapped to it, a conversion to VGPSATCELL = -1/alpha "
        "would be needed); perturb SATCELL_VG instead"
    ),
}
# Feddes (vegetation) parameters: the rows of the Feddes table are the vegetation classes and
# pert_control_name is ignored. ZROOT and PCREF take one entry per row when the table has
# several rows; PCWLT, OMGC and PZ take ONE scalar nominal applied (and sampled independently)
# to every row, a list fails inside the sampler. per_type is a single value for all of them.
FEDDES_LIST = ("ZROOT", "PCREF")
FEDDES_SCALAR = ("PCWLT", "OMGC", "PZ")
# controls each parameter really supports (measured with perturbate() on a stub simulation,
# see tests/test_scenario_layouts.py); the others fail with an obscure error inside numpy
SUPPORTED_CONTROLS = {
    "ic": (None, "layers", "zone"),
    "WTPOSITION": (None,),
    "porosity": (None, "layers", "zone"),
}
LIST_KEYS = ("per_type", "per_name", "per_nom", "per_mean", "per_sigma", "sampling_type", "transf_type")
CONTROLS = (None, "None", "layers", "zone", "root_map")
SAMPLING = ("normal", "lognormal", "uniform")
PER_TYPES = (None, "None", "additive", "multiplicative")
TRANSFORMS = (None, "None", "log", "log-ratio", "hyperbolic")
STATE_UPDATE = "St. var."
NODATA = -9999


class ScenarioError(ValueError):
    """Raised with every error found in the scenario, one per line."""


def _as_none(x):
    return None if (x is None or (isinstance(x, str) and x.lower() == "none")) else x


def _suggest(name: str, choices: Sequence[str]) -> str:
    by_lower = {str(c).lower(): str(c) for c in choices}  # case-insensitive: Zroot -> ZROOT
    close = difflib.get_close_matches(str(name).lower(), list(by_lower), n=1)
    return f" (did you mean '{by_lower[close[0]]}'?)" if close else ""


def _classes(raster) -> Optional[np.ndarray]:
    """Sorted class values of a zone or vegetation raster, nodata excluded."""
    if raster is None:
        return None
    values = np.unique(np.asarray(raster))
    return values[values != NODATA]


def validate_scenario(
    scenario: Dict[str, Any],
    nlayers: Optional[int] = None,
    zone_raster=None,
    veg_raster=None,
    raise_on_error: bool = True,
    n_feddes_rows: Optional[int] = None,
) -> Tuple[List[str], List[str]]:
    """Check a scenario dictionary. Returns ``(errors, warnings)``.

    Parameters
    ----------
    scenario : dict
        The scenario as read from the TOML file (``scenario[name]``).
    nlayers : int, optional
        Number of soil layers, to check the length of per-layer lists.
    zone_raster, veg_raster : array-like, optional
        Zone map and vegetation map, to check per-zone and per-class lists and the
        numbering of the classes (they must be 1..K without gaps).
    raise_on_error : bool
        Raise ``ScenarioError`` listing all errors (default) instead of returning them.
    n_feddes_rows : int, optional
        Number of rows of the Feddes table (vegetation classes), to check the layout of
        ``ZROOT`` and ``PCREF``.
    """
    errors: List[str] = []
    warnings: List[str] = []

    # ---- 1. keys and list lengths ------------------------------------------------
    for key in ("per_name", "per_nom", "per_mean", "per_sigma"):
        if key not in scenario:
            errors.append(f"missing key '{key}'")
    if errors:
        return _finish(errors, warnings, raise_on_error)

    n = len(scenario["per_name"])
    for key in LIST_KEYS + ("per_bounds", "pert_control_name"):
        if key in scenario and len(scenario[key]) != n:
            errors.append(
                f"'{key}' has {len(scenario[key])} entries but per_name has {n}: "
                "all lists are aligned with per_name"
            )

    if errors:  # misaligned lists: the per-parameter checks below would index out of range
        return _finish(errors, warnings, raise_on_error)

    names = list(scenario["per_name"])
    for dup in {x for x in names if names.count(x) > 1}:
        errors.append(f"parameter '{dup}' appears more than once in per_name")

    # ---- 2. per parameter --------------------------------------------------------
    n_zone_classes = None if zone_raster is None else len(_classes(zone_raster))
    n_veg_classes = None if veg_raster is None else len(_classes(veg_raster))

    for i, name in enumerate(names):
        if name not in KNOWN_PARAMETERS:
            errors.append(
                f"unknown parameter '{name}'{_suggest(name, KNOWN_PARAMETERS)}: it would be "
                "silently ignored"
            )
        if name in NOT_IMPLEMENTED:
            errors.append(f"{name}: {NOT_IMPLEMENTED[name]}")
            continue
        control = _as_none(scenario.get("pert_control_name", [None] * n)[i])
        if control not in CONTROLS:
            errors.append(f"{name}: unknown pert_control_name '{control}'; use {CONTROLS[2:]}")

        nom = scenario["per_nom"][i]
        listed = isinstance(nom, list)
        feddes = name in FEDDES_LIST + FEDDES_SCALAR

        if feddes:
            if control is not None:
                warnings.append(
                    f"{name}: pert_control_name '{control}' is ignored for Feddes parameters "
                    "(their rows are the vegetation classes)"
                )
            if "per_type" in scenario and isinstance(scenario["per_type"][i], list):
                errors.append(f"{name}: per_type is a single value for Feddes parameters, not a list")
            if name in FEDDES_SCALAR and listed:
                errors.append(
                    f"{name}: per_nom, per_mean and per_sigma must be scalars (one nominal applied to "
                    "every vegetation class, sampled independently per class); a list fails in the sampler"
                )
            if name in FEDDES_LIST and n_feddes_rows is not None:
                if n_feddes_rows > 1 and not (listed and len(nom) == n_feddes_rows):
                    errors.append(
                        f"{name}: the Feddes table has {n_feddes_rows} rows (vegetation classes): "
                        f"per_nom, per_mean and per_sigma need a list with {n_feddes_rows} entries"
                    )
                elif n_feddes_rows == 1 and listed:
                    errors.append(f"{name}: the Feddes table has 1 row: per_nom, per_mean and per_sigma must be scalars")
        else:
            supported = SUPPORTED_CONTROLS.get(name)
            if supported is not None and control not in supported:
                errors.append(
                    f"{name}: pert_control_name '{control}' is not supported for this parameter "
                    f"(use {[c for c in supported]}); it fails inside the sampler"
                )
            elif control is not None and not listed:
                errors.append(
                    f"{name}: pert_control_name '{control}' needs one entry per "
                    f"{'layer' if control == 'layers' else 'zone' if control == 'zone' else 'class'} "
                    "(a list) in per_nom, per_mean and per_sigma"
                )
            elif control is None and listed:
                errors.append(
                    f"{name}: per_nom is a list but pert_control_name is not set: a list needs "
                    "'layers', 'zone' or 'root_map' (otherwise it fails inside numpy)"
                )
        if listed and not feddes:
            expected = {"layers": nlayers, "zone": n_zone_classes, "root_map": n_veg_classes}.get(control)
            for key in ("per_nom", "per_mean", "per_sigma", "per_type"):
                if key in scenario:
                    v = scenario[key][i]
                    if not isinstance(v, list):
                        if key != "per_type":
                            errors.append(f"{name}: '{key}' must be a list like per_nom")
                        continue
                    if len(v) != len(nom):
                        errors.append(f"{name}: '{key}' has {len(v)} entries, per_nom has {len(nom)}")
            if expected is not None and len(nom) != expected:
                what = {"layers": "layers", "zone": "zones", "root_map": "vegetation classes"}[control]
                errors.append(f"{name}: {len(nom)} entries in the lists but the model has {expected} {what}")

        # sigma
        sig = scenario["per_sigma"][i]
        for v in (sig if isinstance(sig, list) else [sig]):
            if not isinstance(v, (int, float)) or v < 0:
                errors.append(f"{name}: per_sigma must be a number >= 0, got {v!r}")
        # sampling / type / transform
        for key, allowed in (("sampling_type", SAMPLING), ("per_type", PER_TYPES), ("transf_type", TRANSFORMS)):
            if key in scenario:
                v = scenario[key][i]
                for item in (v if isinstance(v, list) else [v]):
                    if item not in allowed:
                        errors.append(f"{name}: {key} '{item}' not in {allowed}")
        # bounds (same parser as the perturbation itself)
        if "per_bounds" in scenario and n == len(scenario["per_bounds"]):
            from pyCATHY.DA.perturbate import parse_bounds

            spec = scenario["per_bounds"][i]
            specs = spec if (isinstance(spec, list) and not (len(spec) == 2 and all(isinstance(x, (int, float)) for x in spec))) else [spec]
            for s in specs:
                try:
                    lo, hi = parse_bounds(s, where=f"per_bounds for {name}")
                except ValueError as exc:
                    errors.append(str(exc))
                    continue
                mean = scenario["per_mean"][i]
                if not isinstance(mean, list):
                    if (lo is not None and mean < lo) or (hi is not None and mean > hi):
                        errors.append(f"{name}: per_mean {mean} is outside the bounds [{lo}, {hi}]")
                    if (lo is not None and hi is not None and isinstance(sig, (int, float)) and sig > 0
                            and (hi - lo) < 2 * sig):
                        warnings.append(f"{name}: bounds [{lo}, {hi}] are narrower than 2 sd ({sig}): "
                                        "the truncated prior will be far from normal")

    # ---- 3. class numbering (B20, B21) ------------------------------------------
    for label, raster, uses in (("zone", zone_raster, "zone"), ("vegetation", veg_raster, "root_map")):
        controls = [_as_none(c) for c in scenario.get("pert_control_name", [])]
        if raster is not None and uses in controls:
            vals = _classes(raster)
            if len(vals) and not np.array_equal(vals, np.arange(1, len(vals) + 1)):
                errors.append(
                    f"{label} map classes are {vals.tolist()}: they must be numbered 1..{len(vals)} "
                    "without gaps (entry k of the lists is written to class k + 1)"
                )

    # ---- 4. which parameters are updated ----------------------------------------
    for upd in scenario.get("listUpdateParm", []):
        if upd == STATE_UPDATE:
            continue
        root = re.split(r"(\d+)$", str(upd))[0]
        if root not in names:
            errors.append(
                f"listUpdateParm: '{upd}' is not perturbed (per_name: {names}); a parameter must be "
                f"perturbed to be updated{_suggest(root, names)}"
            )
            continue
        m = re.search(r"(\d+)$", str(upd))
        nom = scenario["per_nom"][names.index(root)]
        if m and isinstance(nom, list) and int(m.group(1)) >= len(nom):
            errors.append(f"listUpdateParm: '{upd}' but {root} has {len(nom)} entries (0 to {len(nom) - 1})")
    if scenario.get("listUpdateParm") and scenario["listUpdateParm"].count(STATE_UPDATE) > 1:
        warnings.append("listUpdateParm lists 'St. var.' more than once")

    return _finish(errors, warnings, raise_on_error)


def _finish(errors, warnings, raise_on_error):
    if errors and raise_on_error:
        raise ScenarioError("invalid scenario:\n  - " + "\n  - ".join(errors))
    return errors, warnings


def check_update_parameters(list_pert, list_update):
    """Every name of ``listUpdateParm`` (except ``"St. var."``) must be a perturbed key.

    ``list_pert`` is the output of ``perturbate.perturbate``: its ``type_parm`` are the
    keys the update looks up. A name that is not among them used to be dropped without
    a message, so the parameter was perturbed but never updated.
    """
    built = [e["type_parm"] for e in list_pert]
    missing = [u for u in list_update if STATE_UPDATE not in str(u) and u not in built]
    if missing:
        hints = "; ".join(f"'{u}'{_suggest(u, built)}" for u in missing)
        raise ScenarioError(
            f"listUpdateParm names that were not built: {hints}. Perturbed keys: {built}"
        )


def select_zone_parameters(list_update, zone_index, state_label=STATE_UPDATE):
    """Entries of ``listUpdateParm`` updated by the analysis of zone number ``zone_index``.

    Pairing is **by name**: the parameter ``<name><k>`` belongs to zone (or vegetation
    class) ``k``, 0-based, the same index as in the perturbation keys (``ZROOT0`` is the
    first class). The state entries (``"St. var."``) are kept for every zone.

    Legacy fallback, used only when a parameter name has no numeric suffix (e.g. a scalar
    ``Ks``): the parameter at position ``zone_index + 1`` after a leading state entry, as
    before. That positional rule silently pairs ``["St. var.", "ZROOT1", "ZROOT0"]`` the
    wrong way round, which is why names are preferred.

    A zone without a parameter of its own is analysed for the states only.
    """
    states = [u for u in list_update if state_label in str(u)]
    params = [u for u in list_update if state_label not in str(u)]
    suffixes = [re.search(r"(\d+)$", str(u)) for u in params]
    if params and all(m is not None for m in suffixes):
        return states + [u for u, m in zip(params, suffixes) if int(m.group(1)) == zone_index]
    # legacy positional pairing
    if len(list_update) > 1:
        return [list_update[0], list_update[zone_index + 1]]
    return list(list_update[:1])


# ---------------------------------------------------------------------------------
# Resolved configuration: what each perturbed parameter touches
# ---------------------------------------------------------------------------------
FEDDES_PARAMETERS = ("ZROOT", "PCREF", "PCWLT", "OMGC", "PZ")
CONDITIONS = ("ic", "WTPOSITION", "atmbc")


def describe_ensemble(list_pert, list_update=(), soil_index=None, zone_raster=None, veg_raster=None):
    """Table of what every perturbed parameter touches, to catch a layout mistake before a run.

    Parameters
    ----------
    list_pert : list of dict
        Output of ``perturbate.perturbate`` (one entry per layer, zone, class or scalar).
    list_update : sequence of str
        ``listUpdateParm`` (``"St. var."`` is ignored here).
    soil_index : pandas.MultiIndex, optional
        Index ``(zone, layer)`` of the soil table (``read_inputs('soil')[0].index``).
    zone_raster, veg_raster : array-like, optional
        Zone map and vegetation map, needed for the ``root_map`` control.

    Returns
    -------
    pandas.DataFrame
        One row per entry: key, parameter, control, target, prior, bounds, transform and
        whether it is updated. ``df.attrs["warnings"]`` lists writes that override each
        other. The soil rows come from the very function that writes the member files
        (``param_tables._assign_soil_value``), so the table cannot disagree with the run.
    """
    import pandas as pd

    from pyCATHY.DA.param_tables import _SOIL_PARM_COLUMNS, _assign_soil_value

    rows, written = [], {}
    for entry in list_pert:
        key = entry["type_parm"]
        m = re.match(r"^(.*?)(\d+)$", key)
        name, index = (m.group(1), int(m.group(2))) if m else (key, None)
        control = _as_none(entry.get("pert_control_name"))

        if name in _SOIL_PARM_COLUMNS:
            cols = _SOIL_PARM_COLUMNS[name]
            if soil_index is None:
                target = "soil table (index not given)"
            else:
                try:
                    table = pd.DataFrame(np.nan, index=soil_index, columns=cols)
                    veg = None if veg_raster is None else np.asarray(veg_raster)
                    zone = None if zone_raster is None else np.asarray(zone_raster)
                    if control == "root_map" and (veg is None or zone is None):
                        raise ValueError("root_map needs the zone and vegetation maps")
                    _assign_soil_value(table, control, index, cols, 1.0, veg, zone)
                    hit = table.dropna(how="all").index
                    target = "soil rows (zone, layer): " + ", ".join(map(str, hit.tolist()))
                    for r in hit:
                        written.setdefault((tuple(cols), r), []).append(key)
                except Exception as exc:  # the same write would fail in the run
                    target = f"ERROR: {type(exc).__name__}: {exc}"
            target += f"  [columns {', '.join(cols)}]"
        elif name in FEDDES_PARAMETERS:
            target = (f"Feddes table row {index} (vegetation class order, 0-based)"
                      if index is not None else "Feddes table, all rows")
        elif name in CONDITIONS:
            target = "initial / boundary condition (all nodes)"
        else:
            target = "other"

        rows.append({
            "key": key, "parameter": name, "control": control or "global", "target": target,
            "mean": entry.get("mean"), "sd": entry.get("sd"),
            "min": entry.get("clip_min"), "max": entry.get("clip_max"),
            "transform": _as_none(entry.get("transf_type")),
            "updated": key in set(list_update),
        })

    df = pd.DataFrame(rows)
    messages = []
    for (cols, row), keys in written.items():
        if len(keys) > 1:
            messages.append(
                f"{', '.join(keys)} all write {'/'.join(cols)} of row {row}: the last one written wins"
            )
    df.attrs["warnings"] = sorted(set(messages))
    return df
